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
    dat <- simu(n, cen, sce)
    dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    n3 <- 500
    dat3 <- lapply(1:n3, function(x) {
        out <- data.frame(id = x, Y = sort(unique(dat$Y)), death = 0)
        tmp <- simu(1, cen, sce)
        if (sce < 2) out <- cbind(out, tmp[,-(1:3)])
        if (sce == 2.1) out <- cbind(out,
                                     z1 = with(tmp, e * (Y < u) + (1 - e) * (Y >= u)),
                                     z2 = tmp$z2, e = tmp$e, u = tmp$u)
        if (sce == 2.2)
            out <- cbind(out,
                         z1 = with(tmp, e * ((u1 <= Y) * (Y < u2) + (u3 <= Y)) + (1 - e) * 
                                        ((Y < u1) + (u2 <= Y) * (Y < u3))), 
                         z2 = tmp$z2, e = tmp$e, u1 = tmp$u1, u2 = tmp$u2, u3 = tmp$u3)
        if (sce == 2.3)
            out <- cbind(out,
                         z1 = with(tmp, k * Y + b), z2 = tmp$z2, k = tmp$k, b = tmp$b)
        if (sce == 2.4) 
            out <- cbind(out,
                         z1 = with(tmp, k * Y + b), z2 = tmp$z2, k = tmp$k, b = tmp$b)            
        if (sce == 2.5)
            out <- cbind(out,
                         z1 = with(tmp, k * Y * (2 * (Y > 5) - 1) + b),
                         z2 = tmp$z2, k = tmp$k, b = tmp$b)
        return(rbind(out, tmp))
    })
    dat.test <- do.call(rbind, lapply(dat3, function(x) x[-nrow(x),]))
    dat0.test <- dat.test[cumsum(with(dat.test, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    rownames(dat0.test) <- NULL
    dat0.test$Y <- dat0.test$id <- dat0.test$death <- NA
    ctrl <- list(CV = TRUE, tau = quantile(dat0$Y, .95), minsp = 3, minsp2 = 1)
    tt <- seq(0, ctrl$tau, length = 100)
    if (sce %in% c(1.1, 1.2))
        fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
    else fm <- Surv(Y, death) ~ z1 + z2
    fit <- rocForest(fm, data = dat, id = id, control = ctrl)
    fit.dcon <- rocForest(fm, data = dat, id = id, control = c(ctrl, splitBy = "dCON"))
    fit.ranger <- ranger(fm, data = dat0)
    fit.ranger.small <- ranger(fm, data = dat0, min.node.size = 10)
    pred <- predict(fit, dat.test)
    pred.dcon <- predict(fit.dcon, dat.test)
    pred.ranger <- predict(fit.ranger, dat0.test)
    pred.ranger.small <- predict(fit.ranger.small, dat0.test)
    fit.rf <- rfsrc(fm, data = dat0)
    fit.rf.small <- rfsrc(fm, data = dat0, nodesize = 10)
    if (sce %in% c(1.1, 1.2)) {
        pred.rf <- exp(-predict(fit.rf, newdata = subset(dat0.test, select = c(z1, z2, z3, z4, z5, z6, z7, z8, z9, z10)))$chf)
        pred.rf.small <- exp(-predict(fit.rf.small, newdata = subset(dat0.test, select = c(z1, z2, z3, z4, z5, z6, z7, z8, z9, z10)))$chf)
    }
    if (!(sce %in% c(1.1, 1.2))) {
        pred.rf <- exp(-predict(fit.rf, newdata = subset(dat0.test, select = c(z1, z2)))$chf)
        pred.rf.small <- exp(-predict(fit.rf.small, newdata = subset(dat0.test, select = c(z1, z2)))$chf)
    }
    err.ranger <- err.rfsrc <- err.ranger.small <- err.rfsrc.small <- err <- err.dcon <- matrix(NA, length(tt), n3)
    bs.ranger.small <- bs.rfsrc.small <- bs.ranger <- bs.rfsrc <- bs <- bs.dcon <- rep(NA, n3)
    int.bs.ranger.small <- int.bs.rfsrc.small <- int.bs.ranger <- int.bs.rfsrc <- int.bs <- int.bs.dcon <- rep(NA, n3)
    for (i in 1:n3) {
        dat.tmp <- dat.tmp2 <- dat.tmp3 <- dat3[[i]]
        dat.tmp <- dat.tmp[-nrow(dat.tmp),]
        attr(dat.tmp, "prepBy") <- attr(dat, "prepBy")
        attr(dat.tmp, "scenario") <- attr(dat, "scenario")
        truth <- trueSurv(dat.tmp)(tt)
        ## absolute error
        err[,i] <- abs(with(pred$pred[[i]], stepfun(Time , c(1, Surv)))(tt) - truth)
        err.dcon[,i] <- abs(with(pred.dcon$pred[[i]], stepfun(Time , c(1, Surv)))(tt) - truth)
        err.rfsrc[,i] <- abs(stepfun(sort(unique(subset(dat0, dat0$death > 0)$Y)), c(1, pred.rf[i,]))(tt) - truth)
        err.ranger[,i] <- abs(with(pred.ranger, stepfun(unique.death.times, c(1, survival[i,])))(tt) - truth)
        err.rfsrc.small[,i] <- abs(stepfun(sort(unique(subset(dat0, dat0$death > 0)$Y)), c(1, pred.rf.small[i,]))(tt) - truth)
        err.ranger.small[,i] <- abs(with(pred.ranger.small, stepfun(unique.death.times, c(1, survival[i,])))(tt) - truth)
        ## Brier score at event time
        dat.tmp2 <- dat.tmp2[nrow(dat.tmp2),]
        bs[i] <- (1 - with(pred$pred[[i]], stepfun(Time , c(1, Surv)))(dat.tmp2$Y) - dat.tmp2$death)^2
        bs.dcon[i] <- (1 - with(pred.dcon$pred[[i]], stepfun(Time , c(1, Surv)))(dat.tmp2$Y) - dat.tmp2$death)^2
        bs.rfsrc[i] <- (1 - stepfun(sort(unique(subset(dat0, dat0$death > 0)$Y)), c(1, pred.rf[i,]))(dat.tmp2$Y) - dat.tmp2$death)^2
        bs.ranger[i] <- (1 - with(pred.ranger, stepfun(unique.death.times, c(1, survival[i,])))(dat.tmp2$Y) - dat.tmp2$death)^2
        bs.rfsrc.small[i] <- (1 - stepfun(sort(unique(subset(dat0, dat0$death > 0)$Y)), c(1, pred.rf.small[i,]))(dat.tmp2$Y) - dat.tmp2$death)^2
        bs.ranger.small[i] <- (1 - with(pred.ranger.small, stepfun(unique.death.times, c(1, survival[i,])))(dat.tmp2$Y) - dat.tmp2$death)^2
        ## integrated Brier score up to event time
        dat.tmp3 <- dat.tmp3[dat.tmp3$Y <= dat.tmp2$Y,]
        int.bs[i] <- mean((1 - with(pred$pred[[i]], stepfun(Time , c(1, Surv)))(dat.tmp3$Y) - dat.tmp3$death)^2)
        int.bs.dcon[i] <- mean((1 - with(pred.dcon$pred[[i]], stepfun(Time , c(1, Surv)))(dat.tmp3$Y) - dat.tmp3$death)^2)
        int.bs.rfsrc[i] <- mean((1 - stepfun(sort(unique(subset(dat0, dat0$death > 0)$Y)), c(1, pred.rf[i,]))(dat.tmp3$Y) - dat.tmp3$death)^2)
        int.bs.ranger[i] <- mean((1 - with(pred.ranger, stepfun(unique.death.times, c(1, survival[i,])))(dat.tmp3$Y) - dat.tmp3$death)^2)
        int.bs.rfsrc.small[i] <- mean((1 - stepfun(sort(unique(subset(dat0, dat0$death > 0)$Y)), c(1, pred.rf.small[i,]))(dat.tmp3$Y) - dat.tmp3$death)^2)
        int.bs.ranger.small[i] <- mean((1 - with(pred.ranger.small, stepfun(unique.death.times, c(1, survival[i,])))(dat.tmp3$Y) - dat.tmp3$death)^2)
    }
    rm(fit)
    rm(fit.dcon)
    rm(fit.ranger)
    rm(fit.ranger.small)
    rm(fit.rf)
    rm(fit.rf.small)
    c(mean(err), mean(err.dcon), mean(err.ranger), mean(err.ranger.small), mean(err.rfsrc), mean(err.rfsrc.small),
      mean(bs), mean(bs.dcon), mean(bs.ranger), mean(bs.ranger.small), mean(bs.rfsrc), mean(bs.rfsrc.small),
      mean(int.bs), mean(int.bs.dcon), mean(int.bs.ranger), mean(int.bs.ranger.small), mean(int.bs.rfsrc), mean(int.bs.rfsrc.small))
}
