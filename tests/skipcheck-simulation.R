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

do.tree <- function(n, cen, sce = 1.1) {
    ## Pre-determined control list
    ctrl <- list(CV = TRUE)
    if (sce %in% c(1.5, 1.1)) {
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
    fit <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = ctrl)
    fit.dcon <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = c(splitBy = "dCON", ctrl))
    fit.ctree <- ctree(Surv(Y, death) ~ z1 + z2, data = dat0)
    fit.rpart <- rpart(Surv(Y, death) ~ z1 + z2, data = dat0)
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

system.time(print(do.tree(100, 0)))
system.time(print(do.tree(200, 0)))

## running simulation
cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.tree"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(rpart)))
invisible(clusterEvalQ(NULL, library(party)))
invisible(clusterEvalQ(NULL, library(partykit)))
sim1.1.100.00 <- parSapply(NULL, 1:100, function(z) do.tree(100, 0, 1.1))
sim1.1.100.25 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do.tree(100, .25, 1.1), error = function(e) rep(NA, 4)))
sim1.1.100.50 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do.tree(100, .50, 1.1), error = function(e) rep(NA, 4)))
sim1.2.100.00 <- parSapply(NULL, 1:100, function(z) do.tree(100, 0, 1.2))
sim1.2.100.25 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do.tree(100, .25, 1.2), error = function(e) rep(NA, 4)))
sim1.2.100.50 <- parSapply(NULL, 1:100, function(z)
    tryCatch(do.tree(100, .50, 1.1), error = function(e) rep(NA, 4)))
stopCluster(cl)


do.Forest <- function(n, cen, sce = 1.1) {
    ## Pre-determined control list
    ctrl <- list(CV = TRUE)
    if (sce %in% c(1.5, 1.1)) {
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
    ## Fitting & predicting/testing
    fit <- rocForest(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = ctrl)
    fit.dcon <- rocForest(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = c(ctrl, splitBy = "dCON"))
    pred <- predict(fit, dat.test)
    pred.dcon <- predict(fit.dcon, dat.test)
    fit.grf <- regression_forest(subset(dat0, select = c(z1, z2)), dat0$Y, mtry = 2, honesty = TRUE)
    w <- get_sample_weights(fit.grf, newdata = dat0.test)
    fit.rf <- rfsrc(Surv(Y, death) ~ z1 + z2, data = dat0)
    pred.rf <- exp(-predict(fit.rf, newdata = subset(dat0.test, select = c("z1", "z2")))$chf)
    ## Get erros
    err.grf <- err.rfsrc <- err <- err.dcon <- matrix(NA, length(tt), n3)
    for (i in 1:n3) {
        dat.tmp <- dat3[[i]]
        attr(dat.tmp, "prepBy") <- attr(dat, "prepBy")
        attr(dat.tmp, "scenario") <- attr(dat, "scenario")
        truth <- trueSurv(dat.tmp)(tt)
        err[,i] <- with(pred$pred[[i]], stepfun(Time , c(1, Surv)))(tt)
        err.dcon[,i] <- with(pred.dcon$pred[[i]], stepfun(Time , c(1, Surv)))(tt)
        sw <- survfit(Surv(dat0$Y, dat0$death) ~ 1, weights = w[i,])
        err.grf[,i] <- stepfun(sw$time, c(1, sw$surv))(tt)
        err.rfsrc[,i] <- stepfun(unique(dat.test$Y), c(1, pred.rf[i,]))(tt)
        ## absolute error
        err.grf[,i] <- abs(err.grf[,i]- truth)
        err.rfsrc[,i] <- abs(err.rfsrc[,i]- truth)
        err[,i] <- abs(err[,i] - truth)
        err.dcon[,i] <- abs(err.dcon[,i] - truth)
    }
    c(mean(err), mean(err.dcon), mean(err.grf), mean(err.rfsrc))
}
