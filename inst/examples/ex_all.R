library(rocTree)

data(simDat)

simDat <- subset(simDat, id %in% c(91, 16, 5, 11))
simDat$id <- c(rep(3, 3), rep(4, 4), rep(2:1, 2:1))
simDat <- simDat[order(simDat$id),]

## simDat <- subset(simDat, id %in% c(91, 16, 5))
## simDat$id <- rep(3:1, 3:1)
## simDat <- simDat[order(simDat$id),]

## ########################################################
## ex_rocTree.R
## ########################################################

rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = simDat, ensemble = FALSE)

## Root                          
##  ¦--2) z1 <= 0.32338          
##  ¦   ¦--4) z1 <= 0.15423*     
##  ¦   °--5) z1 > 0.15423*      
##  °--3) z1 > 0.32338           
##      ¦--6) z2 <= 0.64179      
##      ¦   ¦--12) z2 <= 0.22388*
##      ¦   °--13) z2 > 0.22388* 
##      °--7) z2 > 0.64179*

rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = simDat, ensemble = TRUE)
newdat$z1 <- 1 * (newdat$Time < median(newdat$Time))

head(newdat)

## Predict survival 
pred <- predict(fit, newdat)
plot(pred)

## Predict hazard
pred <- predict(fit, newdat, type = "hazard")
plot(pred)

## ########################################################
## ex_plot_rocTree
## ########################################################
fit <- rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = simDat,
               ensemble = FALSE)
## Plot tree
plot(fit)
## Plot survival estimates at terminal nodes
plot(fit, type = "survival")
## Plot hazard estimates at terminal nodes
plot(fit, type = "haz")
