data(simDat)
fit <- rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = simDat, ensemble = FALSE)

## testing data
newdat <- data.frame(Time = sort(unique(simDat$Time)), 
                     z2 = runif(1))
newdat$z1 <- 1 * (newdat$Time < median(newdat$Time))
head(newdat)

## Predict survival 
pred <- predict(fit, newdat)
plot(pred)

## Predict hazard
pred <- predict(fit, newdat, type = "hazard")
plot(pred)
