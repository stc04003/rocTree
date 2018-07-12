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


sceCtrl <- function(cen, sce) {
    ## Pre-determined control list
    ## tau is set at the 95th percentiles of Y
    ctrl <- list(CV = TRUE)
    if (sce %in% c(1.5, 1.1, 1.7)) {
        if (cen == 0) ctrl <- c(ctrl, tau = .8)
        if (cen == .25) ctrl <- c(ctrl, tau = .6)
        if (cen == .50) ctrl <- c(ctrl, tau = .5)
    }
    if (sce == 1.2) {
        if (cen == 0) ctrl <- c(ctrl, tau = 2)
        if (cen == .25) ctrl <- c(ctrl, tau = 1.7)
        if (cen == .50) ctrl <- c(ctrl, tau = 1)
    }
    if (sce == 1.3) {
        if (cen == 0) ctrl <- c(ctrl, tau = 3)
        if (cen == .25) ctrl <- c(ctrl, tau = 2.4)
        if (cen == .50) ctrl <- c(ctrl, tau = 1.7)
    }
    if (sce == 1.4) {
        if (cen == 0) ctrl <- c(ctrl, tau = 2.5)
        if (cen == .25) ctrl <- c(ctrl, tau = 2)
        if (cen == .50) ctrl <- c(ctrl, tau = 1.2)
    }
    if (sce == 1.6) {
        if (cen == 0) ctrl <- c(ctrl, tau = 1.41)
        if (cen == .25) ctrl <- c(ctrl, tau = 1.28)
        if (cen == .50) ctrl <- c(ctrl, tau = 1.03)
    }    
    if (sce == 1.8) {
        if (cen == 0) ctrl <- c(ctrl, tau = 30)
        if (cen == .25) ctrl <- c(ctrl, tau = 5.4)
        if (cen == .50) ctrl <- c(ctrl, tau = 1.68)
    }
    if (sce == 1.9) {
        if (cen == 0) ctrl <- c(ctrl, tau = 2.26)
        if (cen == .25) ctrl <- c(ctrl, tau = 1.86)
        if (cen == .50) ctrl <- c(ctrl, tau = 1.41)
    }
    if (sce == 2.1) {
        if (cen == 0) ctrl <- c(ctrl, tau = .88)
        if (cen == .25) ctrl <- c(ctrl, tau = .76)
        if (cen == .50) ctrl <- c(ctrl, tau = .64)
    }
    if (sce == 2.2) {
        if (cen == 0) ctrl <- c(ctrl, tau = 1.30)
        if (cen == .25) ctrl <- c(ctrl, tau = 1.10)
        if (cen == .50) ctrl <- c(ctrl, tau = .82)
    }
    if (sce == 2.3) {
        if (cen == 0) ctrl <- c(ctrl, tau = .55)
        if (cen == .25) ctrl <- c(ctrl, tau = .49)
        if (cen == .50) ctrl <- c(ctrl, tau = .42)
    }
    if (sce == 3.1) {
        if (cen == 0) ctrl <- c(ctrl, tau = 15.95)
        if (cen == .25) ctrl <- c(ctrl, tau = 9.37)
        if (cen == .50) ctrl <- c(ctrl, tau = 4.87)
    }
    if (sce == 3.2) {
        if (cen == 0) ctrl <- c(ctrl, tau = 15.80)
        if (cen == .25) ctrl <- c(ctrl, tau = 9.60)
        if (cen == .50) ctrl <- c(ctrl, tau = 4.94)
    }
    if (sce == 3.3) {
        if (cen == 0) ctrl <- c(ctrl, tau = 1.48)
        if (cen == .25) ctrl <- c(ctrl, tau = 1.29)
        if (cen == .50) ctrl <- c(ctrl, tau = .99)
    }
    return(ctrl)
}

do.Forest <- function(n, cen, sce = 1.1) {
    ctrl <- sceCtrl(cen, sce)
    ## data preparation
    dat <- simu(n, cen, sce)
    dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    n3 <- 1000
    dat3 <- lapply(1:n3, function(x) cbind(id = x, simuTest(dat)))
    dat.test <- do.call(rbind, dat3)
    dat0.test <- dat.test[cumsum(with(dat.test, unlist(lapply(split(id, id), length), use.names =FALSE))),]
    rownames(dat0.test) <- NULL
    dat0.test$Y <- dat0.test$id <- dat0.test$death <- NA
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
    ## removing grf
    ## if (sce == 1.5)
    ##     fit.grf <- regression_forest(subset(dat0, select = c(z1, z2, z3, z4, z5)),
    ##                                  dat0$Y, mtry = 2, honesty = TRUE)
    ## if (sce != 1.5)
    ##     fit.grf <- regression_forest(subset(dat0, select = c(z1, z2)), dat0$Y, mtry = 2, honesty = TRUE)
    ## w <- get_sample_weights(fit.grf, newdata = dat0.test)
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
        ## sw <- survfit(Surv(dat0$Y, dat0$death) ~ 1, weights = w[i,])
        ## err.grf[,i] <- stepfun(sw$time, c(1, sw$surv))(tt)
        err.rfsrc[,i] <- stepfun(sort(unique(subset(dat0, dat0$death > 0)$Y)), c(1, pred.rf[i,]))(tt)
        err.ranger[,i] <- with(pred.ranger, stepfun(unique.death.times, c(1, survival[i,])))(tt)
        ## absolute error
        err.rfsrc[,i] <- abs(err.rfsrc[,i]- truth)
        err[,i] <- abs(err[,i] - truth)
        err.dcon[,i] <- abs(err.dcon[,i] - truth)
        err.ranger[,i] <- abs(err.ranger[,i] - truth)
    }
    c(mean(err), mean(err.dcon), mean(err.ranger), mean(err.rfsrc))
}


cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.Forest"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(randomForestSRC)))
invisible(clusterEvalQ(NULL, library(ranger)))

sim1.1.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 1.1))
sim1.1.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 1.1))
sim1.1.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 1.1))

sim1.2.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 1.2))
sim1.2.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 1.2))
sim1.2.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 1.2))

sim1.3.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 1.3))
sim1.3.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 1.3))
sim1.3.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 1.3))

sim1.4.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 1.4))
sim1.4.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 1.4))
sim1.4.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 1.4))

sim1.5.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 1.5))
sim1.5.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 1.5))
sim1.5.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 1.5))

sim1.9.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(100, 0, 1.9), error = function(e) rep(NA, 4)))
sim1.9.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(100, .25, 1.9), error = function(e) rep(NA, 4)))
sim1.9.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(100, .5, 1.9), error = function(e) rep(NA, 4)))

sim1.8.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(100, 0, 1.8), error = function(e) rep(NA, 4)))
sim1.8.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(100, .25, 1.8), error = function(e) rep(NA, 4)))
sim1.8.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(100, .5, 1.8), error = function(e) rep(NA, 4)))

sim1.7.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(100, 0, 1.7), error = function(e) rep(NA, 4)))
sim1.7.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(100, .25, 1.7), error = function(e) rep(NA, 4)))
sim1.7.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(100, .5, 1.7), error = function(e) rep(NA, 4)))

sim1.9.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(200, 0, 1.9), error = function(e) rep(NA, 4)))
sim1.9.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(200, 0, 1.9), error = function(e) rep(NA, 4)))
sim1.9.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.Forest(200, 0, 1.9), error = function(e) rep(NA, 4)))

stopCluster(cl)

rowMeans(sim1.1.100.00) ## 0.09511860 0.09513851 0.19322780 0.12996838
rowMeans(sim1.1.100.25) ## 0.1041594 0.1041606 0.2258248 0.1296001
rowMeans(sim1.1.100.50) ## 0.1204257 0.1203956 0.2315092 0.1345418

rowMeans(sim1.1.200.00) ## 0.06823050 0.06850718 0.21070689 0.12894672
rowMeans(sim1.1.200.25) ## 0.08239985 0.08239461 0.24581324 0.12815634
rowMeans(sim1.1.200.50) ## 0.09634210 0.09635808 0.25589998 0.13260609

sim1.100 <- list(sim1.1.100.00 = sim1.1.100.00, sim1.1.100.25 = sim1.1.100.25, sim1.1.100.50 = sim1.1.100.50,
                 sim1.2.100.00 = sim1.2.100.00, sim1.2.100.25 = sim1.2.100.25, sim1.2.100.50 = sim1.2.100.50,
                 sim1.3.100.00 = sim1.3.100.00, sim1.3.100.25 = sim1.3.100.25, sim1.3.100.50 = sim1.3.100.50,
                 sim1.4.100.00 = sim1.4.100.00, sim1.4.100.25 = sim1.4.100.25, sim1.4.100.50 = sim1.4.100.50,
                 sim1.5.100.00 = sim1.5.100.00, sim1.5.100.25 = sim1.5.100.25, sim1.5.100.50 = sim1.5.100.50)

save(sim1.100, file = "sim1.100.forest.RData")


cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.Forest"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(randomForestSRC)))
invisible(clusterEvalQ(NULL, library(ranger)))

sim1.1.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 1.1))
sim1.1.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 1.1))
sim1.1.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 1.1))

sim1.2.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 1.2))
sim1.2.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 1.2))
sim1.2.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 1.2))

sim1.3.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 1.3))
sim1.3.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 1.3))
sim1.3.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 1.3))

sim1.4.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 1.4))
sim1.4.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 1.4))
sim1.4.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 1.4))

sim1.5.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 1.5))
sim1.5.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 1.5))
sim1.5.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 1.5))

stopCluster(cl)

rowMeans(sim1.1.200.00) ## 0.09511860 0.09513851 0.19322780 0.12996838
rowMeans(sim1.1.200.25) ## 0.1041594 0.1041606 0.2258248 0.1296001
rowMeans(sim1.1.200.50) ## 0.1204257 0.1203956 0.2315092 0.1345418

rowMeans(sim1.1.200.00) ## 0.06823050 0.06850718 0.21070689 0.12894672
rowMeans(sim1.1.200.25) ## 0.08239985 0.08239461 0.24581324 0.12815634
rowMeans(sim1.1.200.50) ## 0.09634210 0.09635808 0.25589998 0.13260609

sim1.200 <- list(sim1.1.200.00 = sim1.1.200.00, sim1.1.200.25 = sim1.1.200.25, sim1.1.200.50 = sim1.1.200.50,
                 sim1.2.200.00 = sim1.2.200.00, sim1.2.200.25 = sim1.2.200.25, sim1.2.200.50 = sim1.2.200.50,
                 sim1.3.200.00 = sim1.3.200.00, sim1.3.200.25 = sim1.3.200.25, sim1.3.200.50 = sim1.3.200.50,
                 sim1.4.200.00 = sim1.4.200.00, sim1.4.200.25 = sim1.4.200.25, sim1.4.200.50 = sim1.4.200.50,
                 sim1.5.200.00 = sim1.5.200.00, sim1.5.200.25 = sim1.5.200.25, sim1.5.200.50 = sim1.5.200.50)

save(sim1.200, file = "sim1.200.forest.RData")


cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.Forest"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(randomForestSRC)))
invisible(clusterEvalQ(NULL, library(ranger)))

sim2.1.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 2.1))
sim2.1.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 2.1))
sim2.1.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 2.1))

sim2.2.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 2.2))
sim2.2.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 2.2))
sim2.2.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 2.2))

sim2.3.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 2.3))
sim2.3.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 2.3))
sim2.3.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 2.3))

sim2.1.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 2.1))
sim2.1.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 2.1))
sim2.1.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 2.1))

sim2.2.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 2.2))
sim2.2.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 2.2))
sim2.2.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 2.2))

sim2.3.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 2.3))
sim2.3.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 2.3))
sim2.3.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 2.3))

stopCluster(cl)

sim2.100 <- list(sim2.1.100.00 = sim2.1.100.00, sim2.1.100.25 = sim2.1.100.25, sim2.1.100.50 = sim2.1.100.50,
                 sim2.2.100.00 = sim2.2.100.00, sim2.2.100.25 = sim2.2.100.25, sim2.2.100.50 = sim2.2.100.50,
                 sim2.3.100.00 = sim2.3.100.00, sim2.3.100.25 = sim2.3.100.25, sim2.3.100.50 = sim2.3.100.50)

sim2.200 <- list(sim2.1.200.00 = sim2.1.200.00, sim2.1.200.25 = sim2.1.200.25, sim2.1.200.50 = sim2.1.200.50,
                 sim2.2.200.00 = sim2.2.200.00, sim2.2.200.25 = sim2.2.200.25, sim2.2.200.50 = sim2.2.200.50,
                 sim2.3.200.00 = sim2.3.200.00, sim2.3.200.25 = sim2.3.200.25, sim2.3.200.50 = sim2.3.200.50)

save(sim2.100, file = "sim2.100.forest.RData")
save(sim2.200, file = "sim2.200.forest.RData")


cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.Forest"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(randomForestSRC)))
invisible(clusterEvalQ(NULL, library(ranger)))

sim3.1.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 3.1))
sim3.1.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 3.1))
sim3.1.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 3.1))

sim3.2.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 3.2))
sim3.2.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 3.2))
sim3.2.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 3.2))

sim3.3.100.00 <- parSapply(NULL, 1:500, function(z) do.Forest(100, 0, 3.3))
sim3.3.100.25 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .25, 3.3))
sim3.3.100.50 <- parSapply(NULL, 1:500, function(z) do.Forest(100, .5, 3.3))

sim3.1.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 3.1))
sim3.1.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 3.1))
sim3.1.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 3.1))

sim3.2.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 3.2))
sim3.2.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 3.2))
sim3.2.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 3.2))

sim3.3.200.00 <- parSapply(NULL, 1:500, function(z) do.Forest(200, 0, 3.3))
sim3.3.200.25 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .25, 3.3))
sim3.3.200.50 <- parSapply(NULL, 1:500, function(z) do.Forest(200, .5, 3.3))

stopCluster(cl)

sim3.100 <- list(sim3.1.100.00 = sim3.1.100.00, sim3.1.100.25 = sim3.1.100.25, sim3.1.100.50 = sim3.1.100.50,
                 sim3.2.100.00 = sim3.2.100.00, sim3.2.100.25 = sim3.2.100.25, sim3.2.100.50 = sim3.2.100.50,
                 sim3.3.100.00 = sim3.3.100.00, sim3.3.100.25 = sim3.3.100.25, sim3.3.100.50 = sim3.3.100.50)

sim3.200 <- list(sim3.1.200.00 = sim3.1.200.00, sim3.1.200.25 = sim3.1.200.25, sim3.1.200.50 = sim3.1.200.50,
                 sim3.2.200.00 = sim3.2.200.00, sim3.2.200.25 = sim3.2.200.25, sim3.2.200.50 = sim3.2.200.50,
                 sim3.3.200.00 = sim3.3.200.00, sim3.3.200.25 = sim3.3.200.25, sim3.3.200.50 = sim3.3.200.50)

save(sim3.100, file = "sim3.100.forest.RData")
save(sim3.200, file = "sim3.200.forest.RData")


debug(do.Forest)
do.Forest(200, 0, 1.9)

debug(predict)
system.time(predict(fit, dat.test[1:20000,]))

do <- function() {
     W <- oneWeight(ndInd, xlist, object$forest[[1]])
     for (i in 2:length(object$forest)) {
         tmp <- oneWeight(ndInd, xlist, object$forest[[i]])
         W <- lapply(1:nID, function(x) tmp[[x]] + W[[x]])
     }
}

debug(oneWeight)
oneWeight(ndInd, xlist, object$forest[[3]])


giveW2 <- function (ndi, idB2, ndInd2, ndTerm, szL2) {
    n <- length(ndi)
    w <- matrix(0, n, n)
    ind <- ndInd2 == ndi
    rs <- rowSums(ind)
    w[matrix(idB2, length(idB2), n)[t(ind)] + n * (rep(1:n, rs) - 1)] <- rep(1/t(szL2)[which(outer(ndTerm, ndi, "=="))], rs)
    w
}

library(microbenchmark)
microbenchmark(giveW(ndInd[,33], idB2, ndInd2, ndTerm, szL2), giveW2(ndInd[,33], idB2, ndInd2, ndTerm, szL2), times = 1e3)
