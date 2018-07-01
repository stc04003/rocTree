#######################################################################
## Load package
#######################################################################

library(rocTree)
library(survival)

dat <- simu(200, 0, 1.1)
dat

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


tt <- seq(0, 1, .01)

plot(tt, trueHaz(dat)(tt), 's')
for (i in 11:15/10) {
    for (j in 0:2 * .25) {
        lines(tt, trueHaz(simu(200, j, i))(tt) , lwd = 1.2, col = "gray65")
    }
}

