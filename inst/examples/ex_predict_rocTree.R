library(rocTree)
library(survival)

set.seed(1)
dat <- simu(100, 0, 1.3)

fit <- rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = dat, ensemble = FALSE)
plot(fit)

plot(fit, type = "survival")
plot(fit, type = "haz")

debug(plot)

newdat <- dplyr::tibble(Time = sort(unique(dat$Time)), 
                        z1 = 1 * (Time < median(Time)), 
                        z2 = runif(1))
newdat

pred <- predict(fit, newdat)
plot(pred)

pred <- predict(fit, newdat, type = "hazard")
plot(pred)



###############

library(ggplot2)

## survival
tmp <- data.frame(fit$data$.Y0, fit$data$.X0)
colnames(tmp) <- c(fit$rName, fit$vNames)
atTerm <- lapply(split(tmp, fit$nodeLabel + 1), function(x) predict(fit, newdata = x)$pred)
atTerm <- do.call(rbind, atTerm)
atTerm$nd <- as.factor(rep(sort(unique(fit$nodeLabel)) + 1, table(fit$nodeLabel)))
rownames(atTerm) <- NULL

ggplot(atTerm, aes(x = Time, y = Survival, col = nd)) + geom_step(lwd = I(1.1)) +
    xlab("Time") + ylab("Survival probabilities") + labs(col = "Node")


## Hazard
atTerm <- lapply(split(tmp, fit$nodeLabel + 1), function(x)
    predict(fit, newdata = x, type = "haz", control = list(h = .4))$pred)
atTerm <- do.call(rbind, atTerm)
atTerm$nd <- as.factor(rep(sort(unique(fit$nodeLabel)) + 1, table(fit$nodeLabel)))
rownames(atTerm) <- NULL

ggplot(atTerm, aes(x = Time, y = hazard, col = nd)) +
    geom_smooth(method = "loess", se = FALSE) + 
    xlab("Time") + ylab("Hazard estimates") + labs(col = "Nodes")

ggplot(atTerm, aes(x = Time, y = hazard, col = nd)) +
    geom_step(lwd = I(1.1)) +
    xlab("Time") + ylab("Hazard estimates") + labs(col = "Nodes")

