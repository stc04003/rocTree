library(rocTree)
library(survival)

data(simudat)
system.time(fit <- rocTree(Surv(Time, Status) ~ X1 + X2 + X3, id = ID,
                           data = simudat, control = list(CV = TRUE, nflds = 10)))
fit


set.seed(123)
dat <- simu(40, 0, 1.1)
system.time(fit <- rocTree(Surv(Y, death) ~ z1 + z2, id = id, 
                           data = dat, control = list(CV = TRUE, nflds = 10)))
fit

system.time(fit <- rocTree(Surv(Y, death) ~ z1 + z2, id = id, data = dat, control = list(CV = TRUE)))
fit


dat <- simu(50, 0, 1.1)


system.time(fit <- rocForest(Surv(Y, death) ~ z1 + z2, id = id, data = dat, control = list(minsp = 3, minsp2 = 1)))
