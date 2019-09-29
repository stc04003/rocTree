library(rocTree)
library(survival)

set.seed(1)
dat <- simu(100, 0, 1.3)

fit <- rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = dat, ensemble = FALSE)

## Plot tree
plot(fit)

## Plot survival estimates at terminal nodes
plot(fit, type = "survival")

## Plot hazard estimates at terminal nodes
plot(fit, type = "haz")


dat <- simu(n = 100, cen = 0.25, sce = 2.1, summary = TRUE)
system.time(fit2 <- rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = dat, 
                            ensemble = FALSE, control = list(numFold = 10)))
fit2

plot(fit2, type = "survival")

plot(fit2, type = "haz")



tmp <- data.frame(x$data$.Y0, x$data$.X0)[x$data$.D0 > 0,]
colnames(tmp) <- c(x$rName, x$vNames)
t0 <- tmp[,1]

atTerm <- lapply(split(tmp, x$nodeLabel), function(xx) {
    xx <- xx[findInt(t0, xx[,1]),]
    xx[,1] <- t0
    rownames(xx) <- NULL
    return(xx)
})

head(atTerm[[1]], 20)
head(atTerm[[2]], 20)
head(atTerm[[3]], 20)

str(lapply(atTerm, function(xx)
    data.frame(Time = t0, Survival = predict(x, newdata = xx)$survFun(t0))))

str(lapply(atTerm, function(xx) predict(x, newdata = xx)$pred))


atTerm <- do.call(rbind, atTerm)
atTerm$nd <- as.factor(rep(sort(unique(x$nodeLabel)) + 1, each = length(t0)))

