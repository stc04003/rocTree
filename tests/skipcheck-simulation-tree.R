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

#######################################################################
## Trees
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

do.tree <- function(n, cen, sce = 1.1) {
    ctrl <- sceCtrl(cen, sce)
    ## ctrl$minsp <- 20
    ## ctrl$minsp2 <- 0
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
    ## Fitting
    if (sce == 1.5) fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5
    if (sce %in% c(1.6, 1.7, 1.8, 1.9))
        fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
    if (!(sce %in% 15:19/10)) fm <- Surv(Y, death) ~ z1 + z2
    fit <- rocTree(fm, data = dat, id = id, control = ctrl)
    fit.dcon <- rocTree(fm, data = dat, id = id, control = c(splitBy = "dCON", ctrl))
    fit.ctree <- ctree(fm, data = dat0)
    fit.rpart <- rpart(fm, data = dat0)
    ft <- fit.rpart$cptable
    cp <- ft[which.min(ft[,4]), 1]
    fit.rpart <- prune(fit.rpart, cp)
    ## Testing
    pred <- predict(fit, dat.test)
    pred.dcon <- predict(fit.dcon, dat.test)
    pred.ctree <- predict(fit.ctree, dat0.test, type = "prob")
    pred.rpart <- predict(fit.rpart, dat0.test)
    fit.rpart$frame$nd <- 1:dim(fit.rpart$frame)[1]
    df.nd3 <- merge(data.frame(id = 1:n3, pred = pred.rpart), fit.rpart$frame,
                    by.x = "pred", by.y = "yval", sort = FALSE)
    df.nd3 <- df.nd3[order(df.nd3$id),]
    nd3 <- df.nd3$nd
    ## Get errors
    err.ctree <- err.rpart <- err <- err.dcon <- matrix(NA, length(tt), n3)
    for (i in 1:n3) {
        dat.tmp <- dat3[[i]]
        attr(dat.tmp, "prepBy") <- attr(dat, "prepBy")
        attr(dat.tmp, "scenario") <- attr(dat, "scenario")
        truth <- trueSurv(dat.tmp)(tt)
        err.ctree[,i] <- stepfun(pred.ctree[[i]]$time, c(1, pred.ctree[[i]]$surv))(tt)
        km <- survfit(Surv(Y, death) ~ 1, dat0[fit.rpart$where == nd3[i], ])
        err.rpart[,i] <- stepfun(km$time, c(1,km$surv))(tt)
        err[,i] <- stepfun(unique(dat.test$Y), c(1, exp(-cumsum(pred$dfPred[,i]))))(tt)
        err.dcon[,i] <- stepfun(unique(dat.test$Y), c(1, exp(-cumsum(pred.dcon$dfPred[,i]))))(tt)
        ## absolute error
        err.ctree[,i] <- abs(err.ctree[,i]- truth)
        err.rpart[,i] <- abs(err.rpart[,i]- truth)
        err[,i] <- abs(err[,i] - truth)
        err.dcon[,i] <- abs(err.dcon[,i] - truth)
        ## if (i %% 100 == 0) print(i)
    }
    ## sometimes I get errors from ctree
    c(mean(err), mean(err.dcon), mean(err.rpart), mean(err.ctree))
}

## running simulation

## ----------------------------------------------------------------------------------
## Scenario 1, with time-independent covariates
## ----------------------------------------------------------------------------------
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.tree"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(rpart)))
invisible(clusterEvalQ(NULL, library(party)))
invisible(clusterEvalQ(NULL, library(partykit)))
## n = 100
sim1.1.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 1.1))
sim1.1.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 1.1), error = function(e) rep(NA, 4)))
sim1.1.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 1.1), error = function(e) rep(NA, 4)))
##
sim1.2.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 1.2))
sim1.2.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 1.2), error = function(e) rep(NA, 4)))
sim1.2.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 1.2), error = function(e) rep(NA, 4)))
##
sim1.3.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 1.3))
sim1.3.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 1.3), error = function(e) rep(NA, 4)))
sim1.3.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 1.3), error = function(e) rep(NA, 4)))
##
sim1.4.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 1.4))
sim1.4.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 1.4), error = function(e) rep(NA, 4)))
sim1.4.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 1.4), error = function(e) rep(NA, 4)))
##
sim1.5.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 1.5))
sim1.5.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 1.5), error = function(e) rep(NA, 4)))
sim1.5.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 1.5), error = function(e) rep(NA, 4)))
## n = 200
sim1.1.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 1.1))
sim1.1.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 1.1), error = function(e) rep(NA, 4)))
sim1.1.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 1.1), error = function(e) rep(NA, 4)))
##
sim1.2.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 1.2))
sim1.2.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 1.2), error = function(e) rep(NA, 4)))
sim1.2.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 1.2), error = function(e) rep(NA, 4)))
##
sim1.3.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 1.3))
sim1.3.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 1.3), error = function(e) rep(NA, 4)))
sim1.3.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 1.3), error = function(e) rep(NA, 4)))
##
sim1.4.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 1.4))
sim1.4.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 1.4), error = function(e) rep(NA, 4)))
sim1.4.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 1.4), error = function(e) rep(NA, 4)))
##
sim1.5.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 1.5))
sim1.5.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 1.5), error = function(e) rep(NA, 4)))
sim1.5.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 1.5), error = function(e) rep(NA, 4)))
stopCluster(cl)

sim1.100 <- list(sim1.1.100.00 = sim1.1.100.00, sim1.1.100.25 = sim1.1.100.25, sim1.1.100.50 = sim1.1.100.50,
                 sim1.2.100.00 = sim1.2.100.00, sim1.2.100.25 = sim1.2.100.25, sim1.2.100.50 = sim1.2.100.50,
                 sim1.3.100.00 = sim1.3.100.00, sim1.3.100.25 = sim1.3.100.25, sim1.3.100.50 = sim1.3.100.50,
                 sim1.4.100.00 = sim1.4.100.00, sim1.4.100.25 = sim1.4.100.25, sim1.4.100.50 = sim1.4.100.50,
                 sim1.5.100.00 = sim1.5.100.00, sim1.5.100.25 = sim1.5.100.25, sim1.5.100.50 = sim1.5.100.50)

sim1.200 <- list(sim1.1.200.00 = sim1.1.200.00, sim1.1.200.25 = sim1.1.200.25, sim1.1.200.50 = sim1.1.200.50,
                 sim1.2.200.00 = sim1.2.200.00, sim1.2.200.25 = sim1.2.200.25, sim1.2.200.50 = sim1.2.200.50,
                 sim1.3.200.00 = sim1.3.200.00, sim1.3.200.25 = sim1.3.200.25, sim1.3.200.50 = sim1.3.200.50,
                 sim1.4.200.00 = sim1.4.200.00, sim1.4.200.25 = sim1.4.200.25, sim1.4.200.50 = sim1.4.200.50,
                 sim1.5.200.00 = sim1.5.200.00, sim1.5.200.25 = sim1.5.200.25, sim1.5.200.50 = sim1.5.200.50)

save(sim1.100, file = "sim1.100.RData")
save(sim1.200, file = "sim1.200.RData")


## ----------------------------------------------------------------------------------
## Scenario 2, with time-varying covariates
## ----------------------------------------------------------------------------------
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.tree"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(rpart)))
invisible(clusterEvalQ(NULL, library(party)))
invisible(clusterEvalQ(NULL, library(partykit)))
## n = 100
sim2.1.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 2.1))
sim2.1.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 2.1), error = function(e) rep(NA, 4)))
sim2.1.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 2.1), error = function(e) rep(NA, 4)))
##
sim2.2.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 2.2))
sim2.2.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 2.2), error = function(e) rep(NA, 4)))
sim2.2.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 2.2), error = function(e) rep(NA, 4)))
##
sim2.3.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 2.3))
sim2.3.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 2.3), error = function(e) rep(NA, 4)))
sim2.3.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 2.3), error = function(e) rep(NA, 4)))
## n = 200
sim2.1.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 2.1))
sim2.1.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 2.1), error = function(e) rep(NA, 4)))
sim2.1.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 2.1), error = function(e) rep(NA, 4)))
##
sim2.2.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 2.2))
sim2.2.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 2.2), error = function(e) rep(NA, 4)))
sim2.2.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 2.2), error = function(e) rep(NA, 4)))
##
sim2.3.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 2.3))
sim2.3.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 2.3), error = function(e) rep(NA, 4)))
sim2.3.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 2.3), error = function(e) rep(NA, 4)))
stopCluster(cl)

sim2.100 <- list(sim2.1.100.00 = sim2.1.100.00, sim2.1.100.25 = sim2.1.100.25, sim2.1.100.50 = sim2.1.100.50,
                 sim2.2.100.00 = sim2.2.100.00, sim2.2.100.25 = sim2.2.100.25, sim2.2.100.50 = sim2.2.100.50,
                 sim2.3.100.00 = sim2.3.100.00, sim2.3.100.25 = sim2.3.100.25, sim2.3.100.50 = sim2.3.100.50)

sim2.200 <- list(sim2.1.200.00 = sim2.1.200.00, sim2.1.200.25 = sim2.1.200.25, sim2.1.200.50 = sim2.1.200.50,
                 sim2.2.200.00 = sim2.2.200.00, sim2.2.200.25 = sim2.2.200.25, sim2.2.200.50 = sim2.2.200.50,
                 sim2.3.200.00 = sim2.3.200.00, sim2.3.200.25 = sim2.3.200.25, sim2.3.200.50 = sim2.3.200.50)

save(sim2.100, file = "sim2.100.RData")
save(sim2.200, file = "sim2.200.RData")


## ----------------------------------------------------------------------------------
## Scenario 3, with time-varying covariates
## ----------------------------------------------------------------------------------
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.tree"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(rpart)))
invisible(clusterEvalQ(NULL, library(party)))
invisible(clusterEvalQ(NULL, library(partykit)))
## n = 100
sim3.1.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 2.1))
sim3.1.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 2.1), error = function(e) rep(NA, 4)))
sim3.1.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 2.1), error = function(e) rep(NA, 4)))
##
sim3.2.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 2.2))
sim3.2.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 2.2), error = function(e) rep(NA, 4)))
sim3.2.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 2.2), error = function(e) rep(NA, 4)))
##
sim3.3.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 2.3))
sim3.3.100.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .25, 2.3), error = function(e) rep(NA, 4)))
sim3.3.100.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(100, .50, 2.3), error = function(e) rep(NA, 4)))
## n = 200
sim3.1.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 2.1))
sim3.1.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 2.1), error = function(e) rep(NA, 4)))
sim3.1.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 2.1), error = function(e) rep(NA, 4)))
##
sim3.2.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 2.2))
sim3.2.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 2.2), error = function(e) rep(NA, 4)))
sim3.2.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 2.2), error = function(e) rep(NA, 4)))
##
sim3.3.200.00 <- parSapply(NULL, 1:500, function(z) do.tree(200, 0, 2.3))
sim3.3.200.25 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .25, 2.3), error = function(e) rep(NA, 4)))
sim3.3.200.50 <- parSapply(NULL, 1:500, function(z)
    tryCatch(do.tree(200, .50, 2.3), error = function(e) rep(NA, 4)))
stopCluster(cl)

sim3.100 <- list(sim3.1.100.00 = sim3.1.100.00, sim3.1.100.25 = sim3.1.100.25, sim3.1.100.50 = sim3.1.100.50,
                 sim3.2.100.00 = sim3.2.100.00, sim3.2.100.25 = sim3.2.100.25, sim3.2.100.50 = sim3.2.100.50,
                 sim3.3.100.00 = sim3.3.100.00, sim3.3.100.25 = sim3.3.100.25, sim3.3.100.50 = sim3.3.100.50)

sim3.200 <- list(sim3.1.200.00 = sim3.1.200.00, sim3.1.200.25 = sim3.1.200.25, sim3.1.200.50 = sim3.1.200.50,
                 sim3.2.200.00 = sim3.2.200.00, sim3.2.200.25 = sim3.2.200.25, sim3.2.200.50 = sim3.2.200.50,
                 sim3.3.200.00 = sim3.3.200.00, sim3.3.200.25 = sim3.3.200.25, sim3.3.200.50 = sim3.3.200.50)

save(sim3.100, file = "sim3.100.RData")
save(sim3.200, file = "sim3.200.RData")


## ----------------------------------------------------------------------------------
## Additional simulation from scenario I with p = 10
## ----------------------------------------------------------------------------------

cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.tree"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(rpart)))
invisible(clusterEvalQ(NULL, library(party)))
invisible(clusterEvalQ(NULL, library(partykit)))

sim1.6.100.00 <- parSapply(NULL, 1:100, function(z) do.tree(100, 0, 1.6))
sim1.6.100.25 <- parSapply(NULL, 1:100, function(z) 
    tryCatch(do.tree(100, .25, 1.6), error = function(e) rep(NA, 4)))
sim1.6.100.50 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do.tree(100, .50, 1.6), error = function(e) rep(NA, 4)))

sim1.6.200.00 <- parSapply(NULL, 1:100, function(z) do.tree(200, 0, 1.6))
sim1.6.200.25 <- parSapply(NULL, 1:100, function(z) 
    tryCatch(do.tree(200, .25, 1.6), error = function(e) rep(NA, 4)))
sim1.6.200.50 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do.tree(200, .50, 1.6), error = function(e) rep(NA, 4)))

sim1.9.100.00 <- parSapply(NULL, 1:100, function(z) do.tree(100, 0, 1.9))
sim1.9.100.25 <- parSapply(NULL, 1:100, function(z) 
    tryCatch(do.tree(100, .25, 1.9), error = function(e) rep(NA, 4)))
sim1.9.100.50 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do.tree(100, .50, 1.9), error = function(e) rep(NA, 4)))

sim1.9.200.00 <- parSapply(NULL, 1:100, function(z) do.tree(200, 0, 1.9))
sim1.9.200.25 <- parSapply(NULL, 1:100, function(z) 
    tryCatch(do.tree(200, .25, 1.9), error = function(e) rep(NA, 4)))
sim1.9.200.50 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do.tree(200, .50, 1.9), error = function(e) rep(NA, 4)))

stopCluster(cl)

rbind(rowMeans(sim1.6.100.00),
      rowMeans(sim1.6.100.25, na.rm = T),
      rowMeans(sim1.6.100.50, na.rm = T))

rbind(rowMeans(sim1.6.200.00),
      rowMeans(sim1.6.200.25),
      rowMeans(sim1.6.200.50, na.rm = T))

rbind(rowMeans(sim1.9.100.00),
      rowMeans(sim1.9.100.25),
      rowMeans(sim1.9.100.50))

rbind(rowMeans(sim1.9.200.00),
      rowMeans(sim1.9.200.25),
      rowMeans(sim1.9.200.50))




set.seed(1234)
n <- 200
cen <- .25
sce <- 1.2
dat <- simu(n, cen, sce)
n3 <- 1000
dat3 <- lapply(1:n3, function(x) cbind(id = x, simuTest(dat)))
dat.test <- do.call(rbind, dat3)

## Fitting
if (sce == 1.5) fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5
if (sce %in% c(1.6, 1.7, 1.8, 1.9))
    fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
if (!(sce %in% 15:19/10)) fm <- Surv(Y, death) ~ z1 + z2

fit <- rocTree(fm, data = dat, id = id,
               control = list(minsp = 20, minsp2 = 0, CV = TRUE, parallel = TRUE, Trace = TRUE))
fit2 <- rocTree(fm, data = dat, id = id,
                control = list(minsp = 20, minsp2 = 5, CV = TRUE, parallel = TRUE, Trace = TRUE))
fit3 <- rocTree(fm, data = dat, id = id,
                control = list(minsp = 20, minsp2 = 10, CV = TRUE, parallel = TRUE, Trace = TRUE))
fit4 <- rocTree(fm, data = dat, id = id,
                control = list(minsp = 20, minsp2 = 20, CV = TRUE, parallel = TRUE, Trace = TRUE))

head(fit$treeMat, 10)
head(fit2$treeMat, 10)

fit$Frame
fit2$Frame
fit3$Frame
fit4$Frame

fit3 <- rocForest(fm, data = dat, id = id)
fit4 <- rocForest(fm, data = dat, id = id, control = list(minsp2 = 20))

fit3
print(fit3, 2)
print(fit3, 3)
print(fit3, 4)

fit3$forest[[1]]$treeMat
fit3$forest[[2]]$treeMat
fit3$forest[[3]]$treeMat
fit3$forest[[4]]$treeMat
fit3$forest[[5]]$treeMat
fit3$forest[[6]]$treeMat
fit3$forest[[7]]$treeMat
fit3$forest[[8]]$treeMat
fit3$forest[[9]]$treeMat
fit3$forest[[10]]$treeMat
fit3$forest[[11]]$treeMat

fit4$forest[[1]]$treeMat
fit4$forest[[2]]$treeMat
fit4$forest[[3]]$treeMat
fit4$forest[[4]]$treeMat
fit4$forest[[5]]$treeMat
fit4$forest[[6]]$treeMat
fit4$forest[[7]]$treeMat
fit4$forest[[8]]$treeMat
fit4$forest[[9]]$treeMat
fit4$forest[[10]]$treeMat
fit4$forest[[11]]$treeMat

pred <- predict(fit3)
predict(fit3, dat.test)


f1 <- function(A, B, C) {
    if ((A || B) & C) 1
    else 2
}

f2 <- function(B, C) {
    if (B & C) 1
    else 2
}

x1 <- expand.grid(c(TRUE, FALSE), c(TRUE, FALSE), c(TRUE, FALSE))
x2 <- expand.grid(c(TRUE, FALSE), c(TRUE, FALSE))

cbind(x1, apply(x1, 1, function(y) f1(y[1], y[2], y[3])))
cbind(x2, apply(x2, 1, function(y) f2(y[1], y[2])))
