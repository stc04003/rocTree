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

#######################################################################
## Trees
#######################################################################

do.tree <- function(n, cen, sce = 1.1) {
    dat <- simu(n, cen, sce)
    dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    n3 <- 1000
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
    ctrl <- list(CV = TRUE, tau = quantile(dat0$Y, .95), minsp = 20, minsp2 = 5)
    tt <- seq(0, ctrl$tau, length = 100)
    ## Fitting
    if (sce %in% c(1.1, 1.2))
        fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
    else fm <- Surv(Y, death) ~ z1 + z2
    fit <- rocTree(fm, data = dat, id = id, control = ctrl)
    fit.dcon <- rocTree(fm, data = dat, id = id, control = c(splitBy = "dCON", ctrl))
    fit.ctree <- ctree(fm, data = dat0)
    fit.rpart <- rpart(fm, data = dat0)
    dat$Y0 <- with(dat, unlist(lapply(split(Y, id), function(x) c(0, x[-length(x)]))))
    if (sce %in% c(1.1, 1.2))
        fit.cox <- coxph(Surv(Y0, Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10, data = dat)
    else fit.cox <- coxph(Surv(Y0, Y, death) ~ z1 + z2, data = dat)
    ft <- fit.rpart$cptable
    cp <- ft[which.min(ft[,4]), 1]
    fit.rpart <- prune(fit.rpart, cp)
    ## Testing
    pred <- predict(fit, dat.test)
    pred.dcon <- predict(fit.dcon, dat.test)
    pred.ctree <- predict(fit.ctree, dat0.test, type = "prob")
    pred.rpart <- predict(fit.rpart, dat0.test)
    pred.cox <- basehaz(fit.cox, centered = FALSE)
    fit.rpart$frame$nd <- 1:dim(fit.rpart$frame)[1]
    df.nd3 <- merge(data.frame(id = 1:n3, pred = pred.rpart), fit.rpart$frame,
                    by.x = "pred", by.y = "yval", sort = FALSE)
    df.nd3 <- df.nd3[order(df.nd3$id),]
    nd3 <- df.nd3$nd
    ## Get errors
    err.cox <- err.ctree <- err.rpart <- err.dcon <- err <- matrix(NA, length(tt), n3)
    bs.cox <- bs.ctree <- bs.rpart <- bs.dcon <- bs <- rep(NA, n3)
    int.bs.cox <- int.bs.ctree <- int.bs.rpart <- int.bs.dcon <- int.bs <- rep(NA, n3)
    for (i in 1:n3) {
        dat.tmp <- dat.tmp2 <- dat.tmp3 <- dat3[[i]]
        dat.tmp <- dat.tmp[-nrow(dat.tmp),]
        attr(dat.tmp, "prepBy") <- attr(dat, "prepBy")
        attr(dat.tmp, "scenario") <- attr(dat, "scenario")
        ## absolute error
        truth <- trueSurv(dat.tmp)(tt)
        err.ctree[,i] <- abs(stepfun(pred.ctree[[i]]$time, c(1, pred.ctree[[i]]$surv))(tt) - truth)
        km <- survfit(Surv(Y, death) ~ 1, dat0[fit.rpart$where == nd3[i], ])
        err.rpart[,i] <- abs(stepfun(km$time, c(1,km$surv))(tt) - truth)
        err[,i] <- abs(stepfun(unique(dat.test$Y), c(1, exp(-cumsum(pred$dfPred[,i]))))(tt) - truth)
        err.dcon[,i] <- abs(stepfun(unique(dat.test$Y), c(1, exp(-cumsum(pred.dcon$dfPred[,i]))))(tt) - truth)
        ## Brier score at event time
        dat.tmp2 <- dat.tmp2[nrow(dat.tmp2),]
        bs[i] <- (1 - stepfun(unique(dat.test$Y), c(1, exp(-cumsum(pred$dfPred[,i]))))(dat.tmp2$Y) - dat.tmp2$death)^2
        bs.dcon[i] <- (1 - stepfun(unique(dat.test$Y), c(1, exp(-cumsum(pred.dcon$dfPred[,i]))))(dat.tmp2$Y) - dat.tmp2$death)^2
        bs.ctree[i] <- (1 - stepfun(pred.ctree[[i]]$time, c(1, pred.ctree[[i]]$surv))(dat.tmp2$Y) - dat.tmp2$death)^2
        bs.rpart[i] <- (1 - stepfun(km$time, c(1,km$surv))(dat.tmp2$Y) - dat.tmp2$death)^2
        ## integrated Brier score up to event time
        dat.tmp3 <- dat.tmp3[dat.tmp3$Y <= dat.tmp2$Y,]
        int.bs[i] <- mean((1 - stepfun(unique(dat.test$Y), c(1, exp(-cumsum(pred$dfPred[,i]))))(dat.tmp3$Y) - dat.tmp3$death)^2)
        int.bs.dcon[i] <- mean((1 - stepfun(unique(dat.test$Y), c(1, exp(-cumsum(pred.dcon$dfPred[,i]))))(dat.tmp3$Y) - dat.tmp3$death)^2)
        int.bs.ctree[i] <- mean((1 - stepfun(pred.ctree[[i]]$time, c(1, pred.ctree[[i]]$surv))(dat.tmp3$Y) - dat.tmp3$death)^2)
        int.bs.rpart[i] <- mean((1 - stepfun(km$time, c(1,km$surv))(dat.tmp3$Y) - dat.tmp3$death)^2)
        ## Time varying coxph
        if (sce %in% c(1.1, 1.2)) 
            pred.cox.surv <- stepfun(pred.cox$time, c(1, exp(-pred.cox$hazard * exp(with(dat.tmp, cbind(z1, z2, z3, z4, z5, z6, z7, z8, z9, z10)) %*% coef(fit.cox)))))
        else pred.cox.surv <- stepfun(pred.cox$time, c(1, exp(-pred.cox$hazard * exp(with(dat.tmp, cbind(z1, z2)) %*% coef(fit.cox)))))
        err.cox[,i] <- abs(pred.cox.surv(tt) - truth)
        bs.cox[i] <- (1 - pred.cox.surv(dat.tmp2$Y) - dat.tmp2$death)^2
        int.bs.cox[i] <- mean((1 - pred.cox.surv(dat.tmp3$Y) - dat.tmp3$death)^2)
    }
    rm(fit)
    rm(fit.dcon)
    rm(fit.rpart)
    rm(fit.ctree)
    c(mean(err), mean(err.dcon), mean(err.rpart), mean(err.ctree), mean(err.cox),
      mean(bs), mean(bs.dcon), mean(bs.rpart), mean(bs.ctree), mean(bs.cox),
      mean(int.bs), mean(int.bs.dcon), mean(int.bs.rpart), mean(int.bs.ctree), mean(int.bs.cox))
}
