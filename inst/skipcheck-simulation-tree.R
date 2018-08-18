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
    dat <- simu(n, cen, sce)
    dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    n3 <- 1000
    dat3 <- lapply(1:n3, function(x) cbind(id = x, simuTest(dat)))
    dat.test <- do.call(rbind, dat3)
    dat0.test <- dat.test[cumsum(with(dat.test, unlist(lapply(split(id, id), length), use.names =FALSE))),]
    rownames(dat0.test) <- NULL
    dat0.test$Y <- dat0.test$id <- dat0.test$death <- NA
    ctrl <- list(CV = TRUE, tau = quantile(dat0$Y, .95), minsp = 20, minsp2 = 5)
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

## check package version:
set.seed(1)
do.tree(100, 0, 1.2)
## [1] 0.07637702 0.08631965 0.09369669 0.08908296

set.seed(1)
do.tree(100, 0, 1.3)
## [1] 0.08038948 0.07966596 0.08226746 0.08870899

set.seed(1)
do.tree(100, 0, 1.4)
## [1] 0.08448886 0.07240573 0.10411467 0.10370146

set.seed(1)
do.tree(100, 0, 1.6)
## [1] 0.09223166 0.07464071 0.08503459 0.09073412

set.seed(1)
do.tree(100, 0, 1.9)
## [1] 0.1357332 0.1211384 0.1267884 0.1239036

########################################################################################

cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.tree"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(rpart)))
invisible(clusterEvalQ(NULL, library(party)))
invisible(clusterEvalQ(NULL, library(partykit)))

sim1.1.100.00 <- parSapply(NULL, 1:500, function(z) do.tree(100, 0, 1.1))

stopCluster(cl)
