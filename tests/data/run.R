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

#' Fitting rocTree
set.seed(1)
system.time(fit <- rocTree(fm, data = DF, id = ID,
                           control = list(disc = c(0, 1, 1, 1, 0, 0, 1),
                                          tau = 1.5, minsp = 30, minsp2 = 5)))
## 25 secs
fit

## Root                                        
##  ¦--2) KSC <= 0.39615                       
##  ¦   ¦--4) KSC <= 0.18844                   
##  ¦   ¦   ¦--8) HEMOG <= 0.14775             
##  ¦   ¦   ¦   ¦--16) HEMOG <= 0.053533*      
##  ¦   ¦   ¦   °--17) HEMOG > 0.053533*       
##  ¦   ¦   °--9) HEMOG > 0.14775              
##  ¦   ¦       ¦--18) HEMOG <= 0.45610        
##  ¦   ¦       ¦   ¦--36) AIDS <= 0.0000*     
##  ¦   ¦       ¦   °--37) AIDS > 0.0000*      
##  ¦   ¦       °--19) HEMOG > 0.45610         
##  ¦   ¦           ¦--38) AIDS <= 0.0000*     
##  ¦   ¦           °--39) AIDS > 0.0000*      
##  ¦   °--5) KSC > 0.18844                    
##  ¦       ¦--10) HEMOG <= 0.24411*           
##  ¦       °--11) HEMOG > 0.24411*            
##  °--3) KSC > 0.39615                        
##      ¦--6) OP <= 0.0000                     
##      ¦   ¦--12) CD4 <= 0.46467              
##      ¦   ¦   ¦--24) KSC <= 0.72163          
##      ¦   ¦   ¦   ¦--48) TRT <= 1.0000*      
##      ¦   ¦   ¦   °--49) TRT > 1.0000*       
##      ¦   ¦   °--25) KSC > 0.72163           
##      ¦   ¦       ¦--50) HEMOG <= 0.47966*   
##      ¦   ¦       °--51) HEMOG > 0.47966*    
##      ¦   °--13) CD4 > 0.46467               
##      ¦       ¦--26) HEMOG <= 0.55246        
##      ¦       ¦   ¦--52) KSC <= 0.72163*     
##      ¦       ¦   °--53) KSC > 0.72163       
##      ¦       ¦       ¦--106) CD4 <= 0.78587*
##      ¦       ¦       °--107) CD4 > 0.78587* 
##      ¦       °--27) HEMOG > 0.55246         
##      ¦           ¦--54) HEMOG <= 0.68737*   
##      ¦           °--55) HEMOG > 0.68737     
##      ¦               ¦--110) KSC <= 0.72163*
##      ¦               °--111) KSC > 0.72163* 
##      °--7) OP > 0.0000                      
##          ¦--14) TRT <= 1.0000*              
##          °--15) TRT > 1.0000*     
                               
set.seed(1)
system.time(fit2 <- rocTree(fm, data = DF, id = ID, 
                            control = list(disc = c(0, 1, 1, 1, 0, 0, 1), tau = 1.5,
                                           minsp = 30, minsp2 = 5, CV = TRUE, parallel = TRUE)))
## 29.5 secs
fit2

## ROC-guided survival tree
## node), split
##    * denotes terminal node
## Root                    
##  ¦--2) KSC <= 0.39615*  
##  °--3) KSC > 0.39615    
##      ¦--6) OP <= 0.0000*
##      °--7) OP > 0.0000* 

plot(fit, control = list(savePlot = TRUE))
plot(fit2, control = list(savePlot = TRUE))

plotTreeHaz(fit2)

plotTreeHaz(fit2) + theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
ggsave(filename = "haz.pdf")


set.seed(1)
system.time(fit3 <- rocForest(fm, data = DF, id = ID,
                              control = list(disc = c(0, 1, 1, 1, 0, 0, 1),
                                             tau = 1.5, minsp = 30, minsp2 = 5, parallel = TRUE)))
## 31.9 secs

debug(predict)
system.time(pred.fit3 <- predict(fit3))

