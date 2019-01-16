#' Load packages and prepare data
#'
#' Data has 7 variables. They are
#' Time independent covariates: HEMOG, AIDS, TRT, SEX
#' Time dependent covariates: CD4, KSC, OP
#' ID = 1:467

library(tidyverse)
library(rocTree)
library(survival)
library(ggplot2)
library(rpart)
library(party)
library(partykit)

load("dat.RData")

## Required names
N <- dim(dat)[1]
dat.col <- colnames(dat)
tnames <- dat.col[which(substr(dat.col,1,2) == "T2")]
snames <- substr(tnames,3,length(tnames))
oi.names <- snames[5:25]
oit.names <- tnames[5:25]
OP <- dat[,oit.names] 
CD4 <-  dat[,names(dat)[substr(names(dat), 1, 3) == "CD4"]]
KSC <- dat[,names(dat)[substr(names(dat), 1, 6) == "KSCORE"]]
AGE <- dat$AGE
AIDS <- dat$PREVOI
TRT <- dat$RANDGRP
SEX <- dat$GENDER
Y <- dat$T2DEATH
E <- dat$DEATH
HEMOG <- dat$HEMOBL

## Sorting data
OP <- OP[order(Y),]
CD4 <- CD4[order(Y),]
KSC <- KSC[order(Y),]
AGE <- AGE[order(Y)]
HEMOG <- HEMOG[order(Y)]
AIDs <- AIDS[order(Y)]
TRT <- TRT[order(Y)]
SEX <- SEX[order(Y)]
E <- E[order(Y)]
Y <- Y[order(Y)]

## Put everything in a data frame
DF <- NULL

## pre-specified observation times for time dependent covariates
tKSC <- c(c(0, 14), 30.5 * c(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20))
tCD4 <- c(0, 30.5 * c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20))

DF <- NULL
for (i in 1:length(Y)) {
    tOP <- sort(unique(as.numeric(OP[i,])))
    tOP <- tOP[!is.na(tOP)]
    tOP <- tOP[tOP < Y[i]]
    tmp <- NULL
    if (length(tOP) == 0) {
        tmp <- cbind(i, Y[1:i], 0, AGE[i], HEMOG[i], AIDS[i], TRT[i], SEX[i], 
                     approx(tKSC, KSC[i,], Y[1:i], method = "constant",
                            yleft = KSC[i,1],
                            yright = KSC[i, max(which(!is.na(KSC[i,])))])$y,
                     approx(tCD4, CD4[i,], Y[1:i], method = "constant",
                            yleft = CD4[i,1],
                            yright = CD4[i, max(which(!is.na(CD4[i,])))])$y,
                     0)
    } else {
        tmp <- cbind(i, Y[1:i], 0, AGE[i], HEMOG[i], AIDS[i], TRT[i], SEX[i], 
                     approx(tKSC, KSC[i,], Y[1:i], method = "constant",
                            yleft = KSC[i,1],
                            yright = KSC[i, max(which(!is.na(KSC[i,])))])$y,
                     approx(tCD4, CD4[i,], Y[1:i], method = "constant",
                            yleft = CD4[i,1],
                            yright = CD4[i, max(which(!is.na(CD4[i,])))])$y,
                     rowSums(sapply(tOP, function(x) (Y[1:i] >= x))))
    }
    DF <- rbind(DF, tmp)
}

DF <- data.frame(DF)
names(DF) <- c("ID", "Y", "Status", "AGE", "HEMOG", "AIDS", "TRT", "SEX", "KSC", "CD4", "OP")
DF$Y <- DF$Y / 365.25
DF$Status[cumsum(1:length(Y))] <- E
head(DF)

## add noise to Y to break ties
tmp <- matrix(1:length(Y), length(Y), length(Y))
tmp <- tmp[upper.tri(tmp, diag = TRUE)] * 1e-4
DF$Y <- DF$Y + tmp
head(DF, 20)
length(unique(DF$Y)) == length(unique(DF$ID))

fm <- Surv(Y, Status) ~ HEMOG + AIDS + TRT + SEX + KSC + CD4 + OP

#' ------------------------------------------------------------------------------------------
#' Fitting rocTree
#' ------------------------------------------------------------------------------------------
set.seed(1)
system.time(fit1 <- rocTree(fm, data = DF, id = ID, 
                           control = list(disc = c(0, 1, 1, 1, 0, 0, 1), tau = 1.5,
                                          minsp = 20, minsp2 = 5, CV = TRUE,
                                          parallel = T, parCluster = 8)))
## adding (t) to emphasize they are time-dependent covariates.
fit1$vNames[5] <- "KSC(t)"
fit1$vNames[7] <- "OP(t)"

fit1
## Root                    
##  ¦--2) KSC <= 0.39615*  
##  °--3) KSC > 0.39615
##      ¦--6) OP <= 0.0000*
##      °--7) OP > 0.0000* 

plot(fit1, control = list(shape = "rect"))
## plot(fit1, control = list(savePlot = TRUE, shape = "rect"))

plotTreeHaz(fit1)
plotTreeHaz(fit1) + theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
ggsave(filename = "haz-fit1-hn04.pdf")

plotTreeHaz(fit1, control = list(ghN = .2))
plotTreeHaz(fit1, control = list(ghN = .3))
plotTreeHaz(fit1, control = list(ghN = .4))

#' ------------------------------------------------------------------------------------------
#' Random Forest
#' ------------------------------------------------------------------------------------------
set.seed(1)
system.time(fit3 <- rocForest(fm, data = DF, id = ID,
                              control = list(disc = c(0, 1, 1, 1, 0, 0, 1), ghN = .5,
                                             tau = 1.5, minsp = 3, minsp2 = 1,
                                             parallel = T, parCluster = 16)))
## save(fit3, file = "fit3.RData")
load("fit3.RData")

#############################################################################################
#' Testing set with increasing/decreasing KS
#############################################################################################

dat0 <- tibble(Y = sort(unique(DF$Y)), HEMOG = 12.5, AIDS = 1, TRT = 2, SEX = 1, CD4 = 27, OP = 1)
dim(dat0) # 467
quantile(DF$KSC, .39615)
ecdf(DF$KSC)(75)
ecdf(DF$KSC)(80)

#' KSC ranges from 10 to 100
## datA <- dat0 %>% mutate(ID = 1, KSC = seq(40, 60, length = 467))
datA <- dat0 %>% mutate(ID = 1, KSC = seq(60, 90, length = 467))
datB <- dat0 %>% mutate(ID = 2, KSC = seq(80, 40, length = 467))
datC <- dat0 %>% mutate(ID = 3, KSC = 80)
datD <- dat0 %>% mutate(ID = 4, KSC = seq(80, 40, length = 467)[2])

predA <- predB <- predC <- predD <- NULL
system.time(predA <- predict(fit3, datA, type = "hazard"))
system.time(predB <- predict(fit3, datB, type = "hazard"))
system.time(predC <- predict(fit3, datC, type = "hazard"))

## system.time(predA <- predict(fit3, datA, type = "cumHaz"))
## system.time(predB <- predict(fit3, datB, type = "cumHaz"))
## system.time(predC <- predict(fit3, datC, type = "cumHaz"))

datgg <- rbind(predA$pred[[1]], predB$pred[[1]], predC$pred[[1]])
datgg$patient <- rep(LETTERS[1:3], each = dim(predA$pred[[1]])[1])

gg <- ggplot(datgg, aes(x = Time, y = haz, group = patient)) +
    geom_line(aes(linetype = patient, color = patient), lwd = I(1.1)) +
    xlab("Time") + ylab("Hazard")

gg + theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 


## ggsave(filename = "pred-fit3-hn05.pdf")

#############################################################################################
#' Testing set with all covariate fixed
#############################################################################################
dat0 <- tibble(Y = sort(unique(DF$Y)), HEMOG = 12.5, AIDS = 1, TRT = 2, SEX = 1, CD4 = 27, OP = 1)
dat0 <- tibble(Y = sort(unique(DF$Y)), HEMOG = 12.5, AIDS = 1, TRT = 2, SEX = 1, CD4 = 27, OP = 0)

datA <- dat0 %>% mutate(ID = 1, KSC = 90)
datB <- dat0 %>% mutate(ID = 2, KSC = 80)
datC <- dat0 %>% mutate(ID = 3, KSC = 70)
datD <- dat0 %>% mutate(ID = 4, KSC = 60)
datE <- dat0 %>% mutate(ID = 5, KSC = 50)

predA <- predB <- predC <- predD <- predE <- NULL
system.time(predA <- predict(fit3, datA, type = "hazard"))
system.time(predB <- predict(fit3, datB, type = "hazard"))
system.time(predC <- predict(fit3, datC, type = "hazard"))
system.time(predD <- predict(fit3, datD, type = "hazard"))
system.time(predE <- predict(fit3, datE, type = "hazard"))

datgg <- rbind(predA$pred[[1]], predB$pred[[1]], predC$pred[[1]], predD$pred[[1]], predE$pred[[1]])
datgg$patient <- rep(c("KSC = 90", "KSC = 80", "KSC = 70", "KSC = 60", "KSC = 50"),
                     each = dim(predA$pred[[1]])[1])

gg <- ggplot(datgg, aes(x = Time, y = haz, group = patient)) +
    stat_smooth(aes(linetype = patient, color = patient), lwd = I(1.1), se = FALSE, n = 10) +
    xlab("Time") + ylab("Hazard")

gg + theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank()) + xlim(0, 1.55)

## ggsave(filename = "pred-hn05-sm10-op1.pdf")
## ggsave(filename = "pred-hn05-sm10-op0.pdf")

###########################################################################################
#' Testing sets with fixed KS(t)
#'
#' CD4, KSC, OP are time varying covariates
#' The order of the covariates can be seem from `fm`
#'
#' dat0 has all time independent variable
#' 
#' Fix CD4 and OP at .5
###########################################################################################

dat0 <- tibble(Y = sort(unique(DF$Y)), HEMOG = 12.5, AIDS = 1, TRT = 2, SEX = 1, CD4 = 27, OP = 1)
dat0 <- tibble(Y = sort(unique(DF$Y)),
               HEMOG = apply(fit3$xlist0[[1]], 1,
                             function(e) quantile(e, na.rm = TRUE, type = 1, prob = .5)),
               AIDS = apply(fit3$xlist0[[2]], 1,
                            function(e) quantile(e, na.rm = TRUE, type = 1, prob = .5)),
               TRT = apply(fit3$xlist0[[3]], 1,
                           function(e) quantile(e, na.rm = TRUE, type = 1, prob = .5)),
               SEX = apply(fit3$xlist0[[4]], 1,
                           function(e) quantile(e, na.rm = TRUE, type = 1, prob = .5)),
               CD4 = apply(fit3$xlist0[[6]], 1,
                           function(e) quantile(e, na.rm = TRUE, type = 1, prob = .5)),
               OP = apply(fit3$xlist0[[7]], 1,
                          function(e) quantile(e, na.rm = TRUE, type = 1, prob = .5)))

datA <- dat0 %>%
    mutate(ID = 1, KSC = quantile(fit3$xlist0[[5]][1,], type = 1, prob = .1))
datB <- dat0 %>%
    mutate(ID = 2, KSC = quantile(fit3$xlist0[[5]][1,], type = 1, prob = .3))
datC <- dat0 %>%
    mutate(ID = 3, KSC = quantile(fit3$xlist0[[5]][1,], type = 1, prob = .5))

## Use predict to get the W matrix
predA <- predB <- predC <- predD <- NULL
system.time(predA <- predict(fit3, datA, type = "hazard"))
system.time(predB <- predict(fit3, datB, type = "hazard"))
system.time(predC <- predict(fit3, datC, type = "hazard"))


matk <- sapply(fit3$Y0, function(z)
    fit3$E0 * rocTree:::K3(z, fit3$Y0, fit3$ctrl$ghN) / fit3$ctrl$ghN)

datgg <- data.frame(Time = rep(fit3$Y0, 3),
                    haz = c(colSums(matk * predA$Wi[[1]][,1]),
                            colSums(matk * predB$Wi[[1]][,1]),
                            colSums(matk * predC$Wi[[1]][,1])))
datgg$patient <- rep(c("KSC(t) = 0.1", "KSC(t) = 0.3", "KSC(t) = 0.5"), 
                     each = dim(predA$pred[[1]])[1])

gg <- ggplot(datgg, aes(x = Time, y = haz, group = patient)) +
    geom_line(aes(linetype = patient, color = patient), lwd = I(1.1)) +
    xlab("Time") + ylab("Hazard")

gg + theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank()) + xlim(0, 1.5)
ggsave("constant_haz_fixKSt.pdf")

###########################################################################################
#' Fitting rpart and ctree
#'
#' Ploting with DiagrammeR
###########################################################################################

fm <- Surv(Y, Status) ~ HEMOG + AIDS + TRT + SEX + KSC + CD4 + OP
dat0 <- DF[cumsum(with(DF, unlist(lapply(split(ID, ID), length), use.names = FALSE))),]
rownames(dat0) <- NULL

fit.rpart <- rpart(fm, data = dat0)
fit.ctree <- ctree(fm, data = dat0)

dat0 <- dat0 %>% mutate(HEMOG = ecdf(HEMOG)(HEMOG), KSC = ecdf(KSC)(KSC), CD4 = ecdf(CD4)(CD4))
fit.rpart0 <- rpart(fm, data = dat0)
fit.ctree0 <- ctree(fm, data = dat0)

fit.rpart
fit.rpart0
fit.ctree
fit.ctree0

str(fit.rpart)

library(data.tree)
library(DiagrammeR)

## rpart
root <- Node$new("Root", type = "root", decision = "", nd = 1)
Node2 <- root$AddChild("2) KSC >= 75", type = "interior", nd = 2)
Node3 <- root$AddChild("3) KSC < 75", type = "interior", nd = 3)
Node4 <- Node2$AddChild("4) HEMOG >= 11.25", type = "interior", nd = 4)
Node5 <- Node2$AddChild("5) HEMOG < 11.25", type = "terminal", nd = 5)
Node8 <- Node4$AddChild("8) KSC >= 92.5", type = "terminal", nd = 8)
Node9 <- Node4$AddChild("9) KSC < 92.5", type = "terminal", nd = 9)
Node6 <- Node3$AddChild("6) HEMOG >= 10.35", type = "interior", nd = 6)
Node7 <- Node3$AddChild("7) HEMOG < 10.35", type = "terminal", nd = 7)
Node12 <- Node6$AddChild("12) KSC >= 42.5", type = "interior", nd = 12)
Node13 <- Node6$AddChild("13) KSC < 4.5", type = "interior", nd = 13)
Node24 <- Node12$AddChild("24) CD4 >= 68", type = "terminal", nd = 24)
Node25 <- Node12$AddChild("25) CD4 < 68", type = "interior", nd = 25)
Node50 <- Node25$AddChild("50) HEMOG >= 13.25", type = "terminal", nd = 50)
Node51 <- Node25$AddChild("51) HEMOG < 13.25", type = "terminal", nd = 51)
Node26 <- Node13$AddChild("26) OP >= 0.5", type = "terminal", nd = 26)
Node27 <- Node13$AddChild("27) OP < 0.5", type = "terminal", nd = 27)
SetNodeStyle(root, fontname = 'helvetica', shape = "rectangle")

graph <- ToDiagrammeRGraph(root, direction = "climb", pruneFun = NULL)
render_graph(graph, output = "graph")
export_graph(graph, file_name = "rpart-data.pdf", file_type = "pdf")


nodeOnly = FALSE
savePlot = FALSE
file_name = "pic.pdf"
file_type = "pdf"
export_graph(graph, file_name = control$file_name, file_type = control$file_type)



## rpart
root <- Node$new("Root", type = "root", decision = "", nd = 1)
Node2 <- root$AddChild("2) KSC <= 60", type = "interior", nd = 2)
Node3 <- root$AddChild("3) KSC > 60", type = "interior", nd = 3)
Node4 <- Node2$AddChild("4) HEMOG <= 13.3", type = "interior", nd = 4)
Node5 <- Node2$AddChild("5) HEMOG > 13.3", type = "terminal", nd = 5)
Node6 <- Node3$AddChild("6) HEMOG <= 11", type = "interior", nd = 6)
Node7 <- Node3$AddChild("7) HEMOG > 11", type = "terminal", nd = 7)
Node8 <- Node4$AddChild("8) OP <= 0", type = "interior", nd = 8)
Node9 <- Node4$AddChild("9) OP > 0", type = "terminal", nd = 9)
Node12 <- Node6$AddChild("8) KSC <= 70", type = "terminal", nd = 12)
Node13 <- Node6$AddChild("9) KSC > 70", type = "terminal", nd = 13)
SetNodeStyle(root, fontname = 'helvetica', shape = "rectangle")

graph <- ToDiagrammeRGraph(root, direction = "climb", pruneFun = NULL)
render_graph(graph, output = "graph")
export_graph(graph, file_name = "ctree-data.pdf", file_type = "pdf")
