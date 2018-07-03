#######################################################################
## Load package
#######################################################################

library(rocTree)
library(survival)
library(rpart)
library(party)
library(partykit)

#######################################################################
## Scenario 1.1
#######################################################################

## Preparing training data and testing data
## n3 is the number of subjects in testing data

set.seed(123)
dat <- simu(200, 0, 1.1)
dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
n3 <- 1000
dat3 <- lapply(1:n3, function(x) cbind(id = x, simuTest(dat)))
dat.test <- do.call(rbind, dat3)
dat0.test <- dat.test[cumsum(with(dat.test, unlist(lapply(split(id, id), length), use.names =FALSE))),]
rownames(dat0.test) <- NULL
tt <- seq(0, 1.5, .01)

## rocTree and forest    
system.time(fit <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id))
system.time(fit.cv <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id,
                              control = list(CV = TRUE)))
system.time(fit.dcon <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id,
                                control = list(splitBy = "dCON")))
system.time(fit.dcon.cv <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id,
                                   control = list(splitBy = "dCON", CV = TRUE)))

## ctree from party
fit.ctree <- ctree(Surv(Y, death) ~ z1 + z2, data = dat0)
fit.rpart <- rpart(Surv(Y, death) ~ z1 + z2, data = dat0)
ft <- fit.rpart$cptable
cp <- ft[which.min(ft[,4]), 1]
fit.rpart <- prune(fit.rpart, cp)

## Predicting
system.time(pred.rt <- predict(fit, dat.test))
system.time(pred.ctree <- predict(fit.ctree, dat0.test, type = "prob"))
system.time(pred.fit.rpart <- predict(fit.rpart, dat0.test))

fit.rpart$frame$nd <- 1:dim(fit.rpart$frame)[1]
df.nd3 <- merge(data.frame(id = 1:n3, pred = pred.fit.rpart),
                fit.rpart$frame, by.x = "pred", by.y = "yval", sort = FALSE)
df.nd3 <- df.nd3[order(df.nd3$id),]
nd3 <- df.nd3$nd


## Get error
absErr.ctree <- absErr.rpart <- absErr.rocTree <- matrix(NA, length(tt), n3)
for (i in 1:n3) {
    dat.tmp <- dat3[[i]]
    attr(dat.tmp, "prepBy") <- attr(dat, "prepBy")
    attr(dat.tmp, "scenario") <- attr(dat, "scenario")
    truth <- trueSurv(dat.tmp)(tt)
    absErr.ctree[,i] <- stepfun(pred.ctree[[i]]$time, c(1, pred.ctree[[i]]$surv))(tt)
    km <- survfit(Surv(Y, death) ~ 1, dat0[fit.rpart$where == nd3[i], ])
    absErr.rpart[,i] <- stepfun(km$time, c(1,km$surv))(tt)
    absErr.rocTree[,i] <- stepfun(unique(dat.test$Y), c(1, exp(-cumsum(pred.rt$dfPred[,i]))))(tt)
    absErr.ctree[,i] <- abs(absErr.ctree[,i]- truth)
    absErr.rpart[,i] <- abs(absErr.rpart[,i]- truth)
    absErr.rocTree[,i] <- abs(absErr.rocTree[,i] - truth)
    if (i %% 100 == 0) print(i)
}

rowMeans(absErr.ctree)
rowMeans(absErr.rpart)
rowMeans(absErr.rocTree)



## plots
tt <- seq(0, 1.5, .01)
with(predict(fit, dat.test)$pred, plot(Time, Surv, 's', col = 2))
with(predict(fit.cv, dat.test)$pred, lines(Time, Surv, 's', col = 3))
with(predict(fit.dcon, dat.test)$pred, lines(Time, Surv, 's', col = 4))
with(predict(fit.dcon.cv, dat.test)$pred, lines(Time, Surv, 's', col = 5))
lines(tt, trueSurv(dat.test)(tt), col = 1, lwd = 2)
