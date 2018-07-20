#######################################################################
## Load package
#######################################################################

## general
library(parallel)

## tree
library(survival)
library(rocTree)
library(rpart)
library(party)
library(partykit)

## forest
library(randomForestSRC)
library(grf)
library(ranger)

#######################################################################
## Load function
#######################################################################

do.Forest <- function(n, cen, sce = 1.1) {
    ## data preparation
    dat <- simu(n, cen, sce)
    dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    n3 <- 500
    dat3 <- lapply(1:n3, function(x) cbind(id = x, simuTest(dat)))
    dat.test <- do.call(rbind, dat3)
    dat0.test <- dat.test[cumsum(with(dat.test, unlist(lapply(split(id, id), length), use.names =FALSE))),]
    rownames(dat0.test) <- NULL
    dat0.test$Y <- dat0.test$id <- dat0.test$death <- NA
    ctrl <- list(CV = TRUE, tau = quantile(dat0$Y, .95), minsp = 3, minsp2 = 1)
    tt <- seq(0, ctrl$tau, length = 100)
    ## Define formula
    if (sce == 1.5) fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5
    if (sce %in% c(1.6, 1.7, 1.8, 1.9))
        fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
    if (!(sce %in% 15:19/10)) fm <- Surv(Y, death) ~ z1 + z2
    ## Fitting & predicting/testing    
    fit <- rocForest(fm, data = dat, id = id, control = ctrl)
    fit.dcon <- rocForest(fm, data = dat, id = id, control = c(ctrl, splitBy = "dCON"))
    fit.ranger <- ranger(fm, data = dat0)
    pred <- predict(fit, dat.test)
    pred.dcon <- predict(fit.dcon, dat.test)
    pred.ranger <- predict(fit.ranger, dat0.test)
    fit.rf <- rfsrc(fm, data = dat0)
    if (sce == 1.5)
        pred.rf <- exp(-predict(fit.rf, newdata = subset(dat0.test, select = c(z1, z2, z3, z4, z5)))$chf)
    if (sce %in% c(1.6, 1.7, 1.8, 1.9))
        pred.rf <- exp(-predict(fit.rf, newdata = subset(dat0.test, select = c(z1, z2, z3, z4, z5, z6, z7, z8, z9, 10)))$chf)
    if (!(sce %in% 15:19/10)) 
        pred.rf <- exp(-predict(fit.rf, newdata = subset(dat0.test, select = c(z1, z2)))$chf)
    ## Get errors
    err.ranger <- err.rfsrc <- err <- err.dcon <- matrix(NA, length(tt), n3)
    for (i in 1:n3) {
        dat.tmp <- dat3[[i]]
        attr(dat.tmp, "prepBy") <- attr(dat, "prepBy")
        attr(dat.tmp, "scenario") <- attr(dat, "scenario")
        truth <- trueSurv(dat.tmp)(tt)
        err[,i] <- with(pred$pred[[i]], stepfun(Time , c(1, Surv)))(tt)
        err.dcon[,i] <- with(pred.dcon$pred[[i]], stepfun(Time , c(1, Surv)))(tt)
        err.rfsrc[,i] <- stepfun(sort(unique(subset(dat0, dat0$death > 0)$Y)), c(1, pred.rf[i,]))(tt)
        err.ranger[,i] <- with(pred.ranger, stepfun(unique.death.times, c(1, survival[i,])))(tt)
        ## absolute error
        err.rfsrc[,i] <- abs(err.rfsrc[,i]- truth)
        err[,i] <- abs(err[,i] - truth)
        err.dcon[,i] <- abs(err.dcon[,i] - truth)
        err.ranger[,i] <- abs(err.ranger[,i] - truth)
    }
    rm(fit)
    rm(fit.dcon)
    rm(fit.ranger)
    c(mean(err), mean(err.dcon), mean(err.ranger), mean(err.rfsrc))
}

## check package version:
set.seed(1)
do.Forest(100, 0, 1.3) ## Sce IV in the paper
## [1] 0.06866184 0.06513320 0.13121325 0.14945714

set.seed(1)
do.Forest(100, 0, 1.4) ## Sce V in the paper
## [1] 0.06450766 0.06452592 0.13401107 0.15137860

set.seed(1)
do.Forest(100, 0, 1.9) ## Sce I in the paper
## [1] 0.1263728 0.1301165 0.1703579 0.1843670

set.seed(1)
do.Forest(100, 0, 1.6) ## Sce II in the paper
## [1] 0.06225787 0.05847085 0.11120561 0.12835964
