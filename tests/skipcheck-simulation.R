#######################################################################
## Load package
#######################################################################

library(rocTree)
library(survival)

## scenario 1.1
dat <- simu(200, 0, 1.1)
dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]

    
system.time(foo <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id))
system.time(foo.cv <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = list(CV = TRUE)))

system.time(foo.dcon <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = list(splitBy = "dCON")))
system.time(foo.dcon.cv <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id,
                                   control = list(splitBy = "dCON", CV = TRUE)))

foo
foo.dcon

foo.cv
foo.dcon.cv

system.time(pred.rt <- predict(foo, dat.test))
system.time(pred.rt2 <- predict(foo, dat.test0))


tt <- seq(0, 1.5, .01)
with(pred.rt$pred, plot(Time, Surv, 's'))
lines(tt, trueSurv(dat.test)(tt), col = 2)

with(dat.test, plot(Y, exp(-Y^2 * exp(2 * z1 + 2 * z2)), col = 2, cex = .2))
with(pred.rt$pred, lines(Time, Surv, 's', lwd = 2))



with(pred.rt2$pred, lines(Time, Surv, 's', col = 3))
with(pred.rt2$pred, plot(Time, Surv, 's', col = 2))
e
########################################################################

options(error = recover)
options(error = stop)
rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = list(splitBy = "dCON"))


tt <- seq(0, 1, .01)

