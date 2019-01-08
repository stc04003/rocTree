#' Load packages and prepare data
#'
#' Data has 7 variables. They are
#' Time independent covariates: HEMOG, AIDS, TRT, SEX
#' Time dependent covariates: CD4, KSC, OP
#' ID = 1:467

library(rocTree)
library(survival)
library(ggplot2)

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

#' ------------------------------------------------------------------------------------------
#' Testing set
#' ------------------------------------------------------------------------------------------
library(tidyverse)

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


#' ------------------------------------------------------------------------------------------
#' Forest plot
#' ------------------------------------------------------------------------------------------

system.time(forest1 <- rocForest(fm, data = DF, id = ID,
                                 control = list(disc = c(0, 1, 1, 1, 0, 0, 1), ghN = .5,
                                                tau = 1.5, minsp = 3, minsp2 = 1,
                                                parallel = T, parCluster  = 6)))

system.time(forest2 <- rocForest(fm, data = DF, id = ID,
                                 control = list(disc = c(0, 1, 1, 1, 0, 0, 1), ghN = 1,
                                                tau = 1.5, minsp = 3, minsp2 = 1,
                                                parallel = T, parCluster  = 6)))

system.time(forest3 <- rocForest(fm, data = DF, id = ID,
                                 control = list(disc = c(0, 1, 1, 1, 0, 0, 1), ghN = 1,
                                                tau = 1.5, minsp = 10, minsp2 = 3,
                                                parallel = T, parCluster  = 6)))

system.time(forest4 <- rocForest(fm, data = DF, id = ID,
                                 control = list(disc = c(0, 1, 1, 1, 0, 0, 1), ghN = 1,
                                                tau = 1.5, minsp = 3, minsp2 = 1,
                                                fsz = function(x) round(x / 3), 
                                                parallel = T, parCluster  = 6)))

system.time(forest5 <- rocForest(fm, data = DF, id = ID,
                                 control = list(disc = c(0, 1, 1, 1, 0, 0, 1), ghN = 1,
                                                tau = 1.5, minsp = 3, minsp2 = 1,
                                                fsz = function(x) round(0.8 * x),
                                                parallel = T, parCluster  = 6)))


gg <- function(forest) {
    dat0 <- tibble(Y = sort(unique(DF$Y)), HEMOG = 12.5, AIDS = 1, TRT = 2, SEX = 1, CD4 = 27, OP = 1)
    ## datA <- dat0 %>% mutate(ID = 1, KSC = seq(60, 90, length = 467))
    ## datB <- dat0 %>% mutate(ID = 2, KSC = seq(80, 40, length = 467))
    ## datC <- dat0 %>% mutate(ID = 3, KSC = 80)
    ## datD <- dat0 %>% mutate(ID = 4, KSC = seq(80, 40, length = 467)[2])
    datA <- dat0 %>% mutate(ID = 1, KSC = 90)
    datB <- dat0 %>% mutate(ID = 2, KSC = 70)
    datC <- dat0 %>% mutate(ID = 3, KSC = 50)
    datD <- dat0 %>% mutate(ID = 4, KSC = 30)
    predA <- predB <- predC <- predD <- NULL
    system.time(predA <- predict(forest, datA, type = "hazard"))
    system.time(predB <- predict(forest, datB, type = "hazard"))
    system.time(predC <- predict(forest, datC, type = "hazard"))
    system.time(predD <- predict(forest, datD, type = "hazard"))
    datgg <- rbind(predA$pred[[1]], predB$pred[[1]], predC$pred[[1]], predD$pred[[1]])
    datgg$patient <- rep(LETTERS[1:4], each = dim(predA$pred[[1]])[1])
    gg <- ggplot(datgg, aes(x = Time, y = haz, group = patient)) +
        geom_line(aes(linetype = patient, color = patient), lwd = I(1.1)) +
        xlab("Time") + ylab("Hazard")
    gg + theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()) 
}

gg(fit3)

gg(forest1)
gg(forest2)
gg(forest3)
gg(forest4)
gg(forest5)

head(predA$pred[[1]])
head(predB$pred[[1]])
head(predC$pred[[1]])

debug(rocTree:::predict.rocForest)
rocTree:::predict.rocForest(fit3, datB, type = "hazard")
rocTree:::predict.rocForest(fit3, datB, type = "hazard0")
rocTree:::predict.rocForest(fit3, datC, type = "hazard")


plot(Y0, colSums(matk * W[[1]] / 1000), 's')
invisible(sapply(1:467, function(z) lines(Y0, colSums(matk * W[[1]][,z]) / 1000, 's', col = "gray")))
lines(Y0, colSums(matk * W[[1]] / 1000), 's')
lines(Y0, colSums(matk * W[[1]][,467] / 1000), 's')
lines(Y0, colSums(matk * rowMeans(W[[1]]) / 1000), 's', col = 2)

plot(Y0, colSums(matk * W[[1]][,1] / W.rs[[1]][,1]), 's', ylim = c(0, 0.21))
invisible(sapply(1:467, function(z)
    lines(Y0, colSums(matk * W[[1]][,z] / W.rs[[1]][,z]), 's', col = "gray")))

plot(Y0, colSums(matk * W[[1]] / W.rs[[1]]), 's')
plot(Y0, colSums(matk * W[[1]]), 's')


set.seed(1)
dat <- simu(200, 0, 1.2)
dat.test <- simuTest(dat)
fit <- rocForest(Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10,
                 data = dat, id = id,
                 control =  list(ghN = .5, tau = 1.5, minsp = 3, minsp2 = 1))

haz.pred <- predict(fit, dat.test, type = "hazard")
surv.pred <- predict(fit, dat.test, type = "survival")
cumH.pred <- predict(fit, dat.test, type = "cumHaz")


with(haz.pred$pred[[1]], plot(Time, haz, 's'))

with(surv.pred$pred[[1]], plot(Time, Surv, 's'))
with(surv.pred$pred[[1]], lines(Time, trueSurv(dat.test)(Time), 's', col = 2))

with(cumH.pred$pred[[1]], plot(Time, cumHaz, 's'))
with(cumH.pred$pred[[1]], lines(Time, trueHaz(dat.test)(Time), 's', col = 2))


##############################################################################

set.seed(11)
dat <- simu(10, 0.25, 2.1)
dat.test1 <- simuTest(dat)
dat.test2 <- simuTest(dat)
dat.test3 <- simuTest(dat)
dat.test4 <- simuTest(dat)

foo <- rocForest(Surv(Y, death) ~ z1 + z2, data = dat, id = id,
                 control = list(parallel = T, parCluster = 16,
                                ghN = .5, minsp = 3, minsp2 = 1))

foo.pred1 <- predict(foo, dat.test1, type = "hazard")
foo.pred2 <- predict(foo, dat.test2, type = "hazard")
foo.pred3 <- predict(foo, dat.test3, type = "hazard")
foo.pred4 <- predict(foo, dat.test4, type = "hazard")

with(foo.pred1$pred[[1]], plot(Time, haz, 's', ylim = c(0.8, 1.8)))
with(foo.pred2$pred[[1]], lines(Time, haz, 's', col = 2))
with(foo.pred3$pred[[1]], lines(Time, haz, 's', col = 2))
with(foo.pred4$pred[[1]], lines(Time, haz, 's', col = 2))


undebug(rocTree:::predict.rocForest)
debug(rocTree:::predict.rocForest)
rocTree:::predict.rocForest(predict(foo, dat.test1, type = "hazard"))


rr <- 7
matrix(oneW(ndInd, xlist, object$forest[[rr]])[[1]], 10)
matrix(oneV(ndInd, xlist, object$forest[[rr]])[[1]], 10)
dim(unique(t(oneV(ndInd, xlist, object$forest[[rr]])[[1]])))
rbind(diag(matrix(oneW(ndInd, xlist, object$forest[[rr]])[[1]], 10)),
      oneV(ndInd, xlist, object$forest[[rr]])[[1]][,1])


table(t(sapply(1:500, function(x) dim(unique(t(oneV(ndInd, xlist, object$forest[[x]])[[1]]))))))


matrix(lapply(1:dim(xlist[[1]])[2], function(z) giveW(ndInd[, z], idB2, ndInd2, ndTerm, szL2))[[1]],10)

lapply(1:dim(xlist[[1]])[2], function(z) giveV(ndInd[, z], idB2, ndInd2, ndTerm, szL2))
diag(matrix(giveW(rep(2, 10), idB2, ndInd2, ndTerm, szL2), 10))
diag(matrix(giveW(rep(3, 10), idB2, ndInd2, ndTerm, szL2), 10))
diag(matrix(giveW(ndInd[,1], idB2, ndInd2, ndTerm, szL2), 10))

diag(lapply(1:dim(xlist[[1]])[2], function(z) giveV(ndInd[, z], idB2, ndInd2, ndTerm, szL2))[[1]])


giveND <- function(xx) {
    ndInd <- matrix(1, 467, 1)
    for (i in 1:dim(Frame)[1]) {
        if (Frame$terminal[i] == 0) {
            ndInd[ndInd == Frame$nd[i] & xx[[Frame$p[i]]] <= Frame$cut[i]] <- Frame$nd[i] * 2
            ndInd[ndInd == Frame$nd[i] & xx[[Frame$p[i]]] > Frame$cut[i]] <- Frame$nd[i] * 2 + 1
        }
    }
    return(ndInd)
}


head(ndInd[,1], 10)
xlist0 <- lapply(xlist, function(e) matrix(rep(e[4], 467), ncol = 1))
table(giveND(xlist))
table(giveND(xlist0))


sapply(1:467, function(g) unique(giveND(lapply(xlist, function(e) matrix(rep(e[g], 467), ncol = 1)))))
table(sapply(1:467, function(g) unique(giveND(lapply(xlist, function(e) matrix(rep(e[g], 467), ncol = 1))))))
