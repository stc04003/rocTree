data(simDat)
fit <- rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = simDat, ensemble = FALSE)

## Plot tree
plot(fit)

## Plot survival estimates at terminal nodes
plot(fit, type = "survival")

## Plot hazard estimates at terminal nodes
plot(fit, type = "haz")
