#######################################################################
## Load package
#######################################################################

library(rocTree)
library(survival)


dat <- simu(200, 0, 1.1)
head(dat)
summary(dat)

system.time(foo <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id))
system.time(foo.cv <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = list(CV = TRUE)))

system.time(foo.dcon <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = list(splitBy = "dCON")))
system.time(foo.dcon.cv <- rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id,
                                   control = list(splitBy = "dCON", CV = TRUE)))

foo
foo.dcon

foo.cv
foo.dcon.cv


e
########################################################################

options(error = recover)
options(error = stop)
rocTree(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = list(splitBy = "dCON"))
