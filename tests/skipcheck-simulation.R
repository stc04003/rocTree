#######################################################################
## Load package
#######################################################################

library(rocTree)
library(survival)

## scenario 1.1
set.seed(123)
dat <- simu(200, 0, 1.1)
dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    
system.time(foo <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id))
system.time(foo.cv <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = list(CV = TRUE)))

system.time(foo.dcon <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = list(splitBy = "dCON")))
system.time(foo.dcon.cv <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id,
                                   control = list(splitBy = "dCON", CV = TRUE)))



dat.test <- simuTest(dat)
system.time(pred.rt <- predict(foo, dat.test))

tt <- seq(0, 1.5, .01)
with(predict(foo, dat.test)$pred, plot(Time, Surv, 's', col = 2))
with(predict(foo.cv, dat.test)$pred, lines(Time, Surv, 's', col = 3))
with(predict(foo.dcon, dat.test)$pred, lines(Time, Surv, 's', col = 4))
with(predict(foo.dcon.cv, dat.test)$pred, lines(Time, Surv, 's', col = 5))
lines(tt, trueSurv(dat.test)(tt), col = 1, lwd = 2)
