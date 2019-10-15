set.seed(1)
dat <- simu(100, 0, 1.3)

fit <- rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = dat, ensemble = FALSE)

## Plot tree
plot(fit)

## Plot survival estimates at terminal nodes
plot(fit, type = "survival")

## Plot hazard estimates at terminal nodes
plot(fit, type = "haz")
