library(rocTree)
library(survival)

set.seed(1)
dat <- simu(100, 0, 1.3)

## Fitting a pruned survival tree
rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = dat, ensemble = FALSE)

## Fitting a unpruned survival tree
rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = dat, ensemble = FALSE,
        control = list(numFold = 0))

## Fitting the ensemble algorithm (default)
rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = dat, ensemble = TRUE)

