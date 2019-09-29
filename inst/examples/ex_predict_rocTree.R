library(rocTree)
library(survival)

set.seed(1)
dat <- simu(100, 0, 1.3)

fit <- rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = dat, ensemble = FALSE)

## testing data
newdat <- dplyr::tibble(Time = sort(unique(dat$Time)), 
                        z1 = 1 * (Time < median(Time)), 
                        z2 = runif(1))
newdat

## Predict survival 
pred <- predict(fit, newdat)
plot(pred)

## Predict hazard
pred <- predict(fit, newdat, type = "hazard")
plot(pred)
