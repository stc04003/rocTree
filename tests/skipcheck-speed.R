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



library(MASS)


sim1.8 <- function(n, cen = 0) {
    V <- .75^as.matrix(dist(1:10, upper = TRUE))
    z <- mvrnorm(n = n, mu = rep(0, 10), Sigma = V)
    Time <- sqrt(rexp(n) * exp(-rowSums(z) * .5))
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 7.83)
    if (cen == .50) cens <- runif(n, 0, 1.95)
    Y <- pmin(Time, cens)
    d <- 1 * (Time <= cens) 
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), d[x]))))
    dat <- cbind(dat, z[dat$id,])
    names(dat)[4:13] <- paste("z", 1:10, sep = "")
    dat   
}

y <- unique(sim1.8(10000, .5)$Y)
quantile(y, prob = 93:100 / 100)
