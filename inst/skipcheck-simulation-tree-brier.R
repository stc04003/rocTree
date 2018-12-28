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
## Trees using integrated Brier score as score function
#######################################################################

do.tree <- function(n, cen, sce = 1.1) {
    dat <- simu(n, cen, sce)
    dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    n3 <- 1000
    ## dat3 <- lapply(1:n3, function(x) cbind(id = x, simuTest(dat)))
    dat3 <- lapply(1:n3, function(x) {
        tmp <- subset(simu(n, cen, sce), id == sample(1:n, 1))
        tmp$id <- x
        return(tmp)
    })
    dat.test <- do.call(rbind, dat3)
    dat0.test <- dat.test[cumsum(with(dat.test, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    rownames(dat0.test) <- NULL
    dat0.test$Y <- dat0.test$id <- dat0.test$death <- NA
    ctrl <- list(CV = TRUE, tau = quantile(dat0$Y, .95), minsp = 20, minsp2 = 5)
    ## tt <- seq(0, ctrl$tau, length = 100)
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
    bs.ctree <- bs.rpart <- bs <- bs.dcon <- rep(NA, n3)
    for (i in 1:n3) {
        dat.tmp <- dat3[[i]]
        attr(dat.tmp, "prepBy") <- attr(dat, "prepBy")
        attr(dat.tmp, "scenario") <- attr(dat, "scenario")
        ## truth <- trueSurv(dat.tmp)(tt)
        bs.ctree[i] <- mean((1 - stepfun(pred.ctree[[i]]$time, c(1, pred.ctree[[i]]$surv))(dat.tmp$Y) - dat.tmp$death)^2)
        km <- survfit(Surv(Y, death) ~ 1, dat0[fit.rpart$where == nd3[i], ])
        bs.rpart[i] <- mean((1 - stepfun(km$time, c(1,km$surv))(dat.tmp$Y) - dat.tmp$death)^2)
        bs[i] <- mean((1 - stepfun(dat.tmp$Y, c(1, exp(-cumsum(pred$dfPred[1:nrow(dat.tmp),i]))))(dat.tmp$Y) - dat.tmp$death)^2)
        bs.dcon[i] <- mean((1 - stepfun(dat.tmp$Y, c(1, exp(-cumsum(pred.dcon$dfPred[1:nrow(dat.tmp),i]))))(dat.tmp$Y) - dat.tmp$death)^2)
        ## if (i %% 100 == 0) print(i)
    }
    c(mean(bs), mean(bs.dcon), mean(bs.rpart), mean(bs.ctree))
}


cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do.tree"))
invisible(clusterExport(NULL, "sceCtrl"))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(rpart)))
invisible(clusterEvalQ(NULL, library(party)))
invisible(clusterEvalQ(NULL, library(partykit)))

sim1.1.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 1.1), error = function(e) rep(NA, 4)))
sim1.1.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 1.1), error = function(e) rep(NA, 4)))
sim1.1.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 1.1), error = function(e) rep(NA, 4)))
sim1.1.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 1.1), error = function(e) rep(NA, 4)))
sim1.1.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 1.1), error = function(e) rep(NA, 4)))
sim1.1.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 1.1), error = function(e) rep(NA, 4)))
sim1.1 <- list(sim1.1.100.00 = sim1.1.100.00, sim1.1.100.25 = sim1.1.100.25, sim1.1.100.50 = sim1.1.100.50,
               sim1.1.200.00 = sim1.1.200.00, sim1.1.200.25 = sim1.1.200.25, sim1.1.200.50 = sim1.1.200.50)
save(sim1.1, file = "sim1.1.RData")

sim1.2.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 1.2), error = function(e) rep(NA, 4)))
sim1.2.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 1.2), error = function(e) rep(NA, 4)))
sim1.2.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 1.2), error = function(e) rep(NA, 4)))
sim1.2.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 1.2), error = function(e) rep(NA, 4)))
sim1.2.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 1.2), error = function(e) rep(NA, 4)))
sim1.2.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 1.2), error = function(e) rep(NA, 4)))
sim1.2 <- list(sim1.2.100.00 = sim1.2.100.00, sim1.2.100.25 = sim1.2.100.25, sim1.2.100.50 = sim1.2.100.50,
               sim1.2.200.00 = sim1.2.200.00, sim1.2.200.25 = sim1.2.200.25, sim1.2.200.50 = sim1.2.200.50)
save(sim1.2, file = "sim1.2.RData")

sim1.3.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 1.3), error = function(e) rep(NA, 4)))
sim1.3.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 1.3), error = function(e) rep(NA, 4)))
sim1.3.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 1.3), error = function(e) rep(NA, 4)))
sim1.3.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 1.3), error = function(e) rep(NA, 4)))
sim1.3.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 1.3), error = function(e) rep(NA, 4)))
sim1.3.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 1.3), error = function(e) rep(NA, 4)))
sim1.3 <- list(sim1.3.100.00 = sim1.3.100.00, sim1.3.100.25 = sim1.3.100.25, sim1.3.100.50 = sim1.3.100.50,
               sim1.3.200.00 = sim1.3.200.00, sim1.3.200.25 = sim1.3.200.25, sim1.3.200.50 = sim1.3.200.50)
save(sim1.3, file = "sim1.3.RData")

sim1.4.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 1.4), error = function(e) rep(NA, 4)))
sim1.4.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 1.4), error = function(e) rep(NA, 4)))
sim1.4.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 1.4), error = function(e) rep(NA, 4)))
sim1.4.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 1.4), error = function(e) rep(NA, 4)))
sim1.4.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 1.4), error = function(e) rep(NA, 4)))
sim1.4.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 1.4), error = function(e) rep(NA, 4)))
sim1.4 <- list(sim1.4.100.00 = sim1.4.100.00, sim1.4.100.25 = sim1.4.100.25, sim1.4.100.50 = sim1.4.100.50,
               sim1.4.200.00 = sim1.4.200.00, sim1.4.200.25 = sim1.4.200.25, sim1.4.200.50 = sim1.4.200.50)
save(sim1.4, file = "sim1.4.RData")

sim1.5.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 1.5), error = function(e) rep(NA, 4)))
sim1.5.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 1.5), error = function(e) rep(NA, 4)))
sim1.5.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 1.5), error = function(e) rep(NA, 4)))
sim1.5.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 1.5), error = function(e) rep(NA, 4)))
sim1.5.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 1.5), error = function(e) rep(NA, 4)))
sim1.5.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 1.5), error = function(e) rep(NA, 4)))
sim1.5 <- list(sim1.5.100.00 = sim1.5.100.00, sim1.5.100.25 = sim1.5.100.25, sim1.5.100.50 = sim1.5.100.50,
               sim1.5.200.00 = sim1.5.200.00, sim1.5.200.25 = sim1.5.200.25, sim1.5.200.50 = sim1.5.200.50)
save(sim1.5, file = "sim1.5.RData")

sim2.1.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 2.1), error = function(e) rep(NA, 4)))
sim2.1.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 2.1), error = function(e) rep(NA, 4)))
sim2.1.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 2.1), error = function(e) rep(NA, 4)))
sim2.1.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 2.1), error = function(e) rep(NA, 4)))
sim2.1.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 2.1), error = function(e) rep(NA, 4)))
sim2.1.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 2.1), error = function(e) rep(NA, 4)))
sim2.1 <- list(sim2.1.100.00 = sim2.1.100.00, sim2.1.100.25 = sim2.1.100.25, sim2.1.100.50 = sim2.1.100.50,
               sim2.1.200.00 = sim2.1.200.00, sim2.1.200.25 = sim2.1.200.25, sim2.1.200.50 = sim2.1.200.50)
save(sim2.1, file = "sim2.1.RData")

sim2.2.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 2.2), error = function(e) rep(NA, 4)))
sim2.2.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 2.2), error = function(e) rep(NA, 4)))
sim2.2.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 2.2), error = function(e) rep(NA, 4)))
sim2.2.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 2.2), error = function(e) rep(NA, 4)))
sim2.2.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 2.2), error = function(e) rep(NA, 4)))
sim2.2.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 2.2), error = function(e) rep(NA, 4)))
sim2.2 <- list(sim2.2.100.00 = sim2.2.100.00, sim2.2.100.25 = sim2.2.100.25, sim2.2.100.50 = sim2.2.100.50,
               sim2.2.200.00 = sim2.2.200.00, sim2.2.200.25 = sim2.2.200.25, sim2.2.200.50 = sim2.2.200.50)
save(sim2.2, file = "sim2.2.RData")

sim2.3.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 2.3), error = function(e) rep(NA, 4)))
sim2.3.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 2.3), error = function(e) rep(NA, 4)))
sim2.3.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 2.3), error = function(e) rep(NA, 4)))
sim2.3.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 2.3), error = function(e) rep(NA, 4)))
sim2.3.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 2.3), error = function(e) rep(NA, 4)))
sim2.3.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 2.3), error = function(e) rep(NA, 4)))
sim2.3 <- list(sim2.3.100.00 = sim2.3.100.00, sim2.3.100.25 = sim2.3.100.25, sim2.3.100.50 = sim2.3.100.50,
               sim2.3.200.00 = sim2.3.200.00, sim2.3.200.25 = sim2.3.200.25, sim2.3.200.50 = sim2.3.200.50)
save(sim2.3, file = "sim2.3.RData")

sim2.4.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 2.4), error = function(e) rep(NA, 4)))
sim2.4.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 2.4), error = function(e) rep(NA, 4)))
sim2.4.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 2.4), error = function(e) rep(NA, 4)))
sim2.4.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 2.4), error = function(e) rep(NA, 4)))
sim2.4.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 2.4), error = function(e) rep(NA, 4)))
sim2.4.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 2.4), error = function(e) rep(NA, 4)))
sim2.4 <- list(sim2.4.100.00 = sim2.4.100.00, sim2.4.100.25 = sim2.4.100.25, sim2.4.100.50 = sim2.4.100.50,
               sim2.4.200.00 = sim2.4.200.00, sim2.4.200.25 = sim2.4.200.25, sim2.4.200.50 = sim2.4.200.50)
save(sim2.4, file = "sim2.4.RData")

sim2.5.100.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, 0, 2.5), error = function(e) rep(NA, 4)))
sim2.5.100.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .25, 2.5), error = function(e) rep(NA, 4)))
sim2.5.100.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(100, .50, 2.5), error = function(e) rep(NA, 4)))
sim2.5.200.00 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, 0, 2.5), error = function(e) rep(NA, 4)))
sim2.5.200.25 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .25, 2.5), error = function(e) rep(NA, 4)))
sim2.5.200.50 <- parSapply(NULL, 1:500, function(z) tryCatch(do.tree(200, .50, 2.5), error = function(e) rep(NA, 4)))
sim2.5 <- list(sim2.5.100.00 = sim2.5.100.00, sim2.5.100.25 = sim2.5.100.25, sim2.5.100.50 = sim2.5.100.50,
               sim2.5.200.00 = sim2.5.200.00, sim2.5.200.25 = sim2.5.200.25, sim2.5.200.50 = sim2.5.200.50)
save(sim2.5, file = "sim2.5.RData")

stopCluster(cl)
