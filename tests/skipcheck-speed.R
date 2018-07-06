library(rocTree)
library(survival)
n <- 100
cen <- 0
sce <- 1.1

ctrl <- list(CV = TRUE)
if (sce %in% c(1.5, 1.1)) {
    if (cen == 0) ctrl <- c(ctrl, tau = .8)
    if (cen == .25) ctrl <- c(ctrl, tau = .6)
    if (cen == .50) ctrl <- c(ctrl, tau = .5)
}
if (sce == 1.2) {
    if (cen == 0) ctrl <- c(ctrl, tau = 2)
    if (cen == .25) ctrl <- c(ctrl, tau = 1.7)
    if (cen == .50) ctrl <- c(ctrl, tau = 1)
}
if (sce == 1.3) {
    if (cen == 0) ctrl <- c(ctrl, tau = 3)
    if (cen == .25) ctrl <- c(ctrl, tau = 2.4)
    if (cen == .50) ctrl <- c(ctrl, tau = 1.7)
}
if (sce == 1.4) {
    if (cen == 0) ctrl <- c(ctrl, tau = 2.5)
    if (cen == .25) ctrl <- c(ctrl, tau = 2)
    if (cen == .50) ctrl <- c(ctrl, tau = 1.2)
}

set.seed(1234)
dat <- simu(n, cen, sce)
dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
n3 <- 1000
dat3 <- lapply(1:n3, function(x) cbind(id = x, simuTest(dat)))
dat.test <- do.call(rbind, dat3)
dat0.test <- dat.test[cumsum(with(dat.test, unlist(lapply(split(id, id), length), use.names =FALSE))),]
rownames(dat0.test) <- NULL
dat0.test$Y <- dat0.test$id <- dat0.test$death <- NA
tt <- seq(0, ctrl$tau, length = 100)

fit <- rocForest(Surv(Y, death) ~ z1 + z2, data = dat, id = id, control = ctrl)
system.time(pred2 <- predict(fit, dat.test))
## 137.272 
