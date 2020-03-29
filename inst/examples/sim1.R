## Assumes time-independent covariates

## General
library(survival)
library(rocTree)
library(MASS)

## Tree packages
library(rpart)
library(party)
library(partykit)

## Forest packages
library(randomForestSRC)
library(ranger)


indCov1 <- function (n, cen, tdForm = TRUE, order = TRUE) {
    W <- matrix(mvrnorm(n, rep(0, 25), 0.9^outer(1:25, 1:25, 
        function(x, y) abs(x - y))), n)
    Y0 <- rexp(n, exp(-0.1 * rowSums(W[, 11:25])))
    if (cen == 0) 
        cc <- Inf
    if (cen == 0.25) 
        cc <- rexp(n, 0.25)
    if (cen == 0.5) 
        cc <- rexp(n, 1)
    death <- 1 * (Y0 <= cc)
    Y <- pmin(Y0, cc)
    dat <- data.frame(ID = 1:n, Time = Y, death = death, X = W)
    if (order) 
        dat <- ord(dat)
    if (tdForm) 
        dat <- ti2td(dat)
    attr(dat, "from") <- "indCov1"
    return(dat)
}

indCov2 <- function (n, cen, tdForm = TRUE, order = TRUE) {
    W <- matrix(runif(25 * n), ncol = 25)
    Y0 <- rexp(n, 1/(sin(W[, 1] * pi) + 2 * abs(W[, 2] - 0.5) + 
        W[, 3]^3))
    cc <- runif(n, 0, cen)
    if (cen == 0) 
        cc <- Inf
    if (cen == 0.25) 
        cc <- runif(n, 0, 5.3)
    if (cen == 0.5) 
        cc <- runif(n, 0, 2)
    death <- 1 * (Y0 <= cc)
    Y <- pmin(Y0, cc)
    dat <- data.frame(ID = 1:n, Time = Y, death = death, X = W)
    if (order) 
        dat <- ord(dat)
    if (tdForm) 
        dat <- ti2td(dat)
    attr(dat, "from") <- "indCov2"
    return(dat)
}

indCov3 <- function (n, cen, tdForm = TRUE, order = TRUE) {
    W <- matrix(mvrnorm(n, rep(0, 25), 0.75^outer(1:25, 1:25, 
        function(x, y) abs(x - y))), n)
    Y0 <- rgamma(n, shape = 0.5 + 0.3 * abs(rowSums(W[, 11:15])), 
        scale = 2)
    if (cen == 0) 
        cc <- Inf
    if (cen == 0.25) 
        cc <- runif(n, 0, 11.7)
    if (cen == 0.5) 
        cc <- runif(n, 0, 4.9)
    death <- 1 * (Y0 <= cc)
    Y <- pmin(Y0, cc)
    dat <- data.frame(ID = 1:n, Time = Y, death = death, X = W)
    if (order) 
        dat <- ord(dat)
    if (tdForm) 
        dat <- ti2td(dat)
    attr(dat, "from") <- "indCov3"
    return(dat)
}

ord <- function (dat) {
    dat <- dat[order(dat$ID, dat$Time), ]
    ord0 <- aggregate(Time ~ ID, dat, max)$Time
    dat <- dat[order(rep(ord0, aggregate(Time ~ ID, dat, length)$Time)), ]
    dat$ID <- rep(1:length(unique(dat$ID)), table(dat$ID)[unique(dat$ID)])
    rownames(dat) <- NULL
    return(dat)
}

ti2td <- function(dat) {
    dat2 <- data.frame(ID = dat$ID,
                       Time = c(t(dat$Time * ifelse(outer(dat$Time, dat$Time, "<="), TRUE, NA))),
                       death = c(diag(dat$death)), 
                       dat[, grep("X.", names(dat))])
    dat2 <- dat2[complete.cases(dat2), ]
    dat <- dat2[order(dat2$ID), ]
    rownames(dat) <- NULL
    return(dat)
}

trueSurv <- function (dat, n = 5000) {
    .case <- attr(dat, "from")
    if (!is.null(dat$death)) 
        dat <- subset(dat, death > 0)
    .X <- dat[, grep("X.", names(dat))]
    if (.case == "indCov1") 
        .tt <- replicate(n, rexp(nrow(.X), exp(-0.1 * rowSums(.X[, 11:25]))))
    if (.case == "indCov2") 
        .tt <- replicate(n, rexp(nrow(.X), 1/(sin(.X[, 1] * pi) + 2 * abs(.X[, 2] - 0.5) + .X[, 3]^3)))
    if (.case == "indCov3") 
        .tt <- replicate(n, rgamma(nrow(.X), shape = 0.5 + 0.3 * abs(rowSums(.X[, 11:15])), scale = 2))
    if (.case == "indCov4") 
        .tt <- replicate(n, rlnorm(n, meanlog = 0.1 * (rowSums(.X[,  1:5]) + rowSums(.X[, 21:25]))))
    .tt <- rowMeans(.tt)
    list(Time = .tt[order(.tt)], Surv = 1 - ecdf(.tt)(.tt[order(.tt)]))
}

## ###################################################################################################
## Function to fit all models
## Fit rocForest and rocTree only
## ###################################################################################################

n <- 4
cen <- .25
sce <- 1

do <- function(n, cen, sce) {
    txt1 <- paste("dat <- indCov", sce, "(", n, ",", cen, ")", sep = "")
    eval(parse(text = txt1))
    if (sce == 1) .t0 <- seq(0, 8.211596, length = 1000)
    if (sce == 2) .t0 <- seq(0, 4.45324, length = 1000)
    if (sce == 3) .t0 <- seq(0, 8.677644, length = 1000)
    if (sce == 4) .t0 <- seq(0, 6.736605, length = 1000)
    fit1 <- rocTree(Surv(Time, death) ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7 + X.8 + X.9 + X.10 +
                        X.11 + X.12 + X.13 + X.14 + X.15 + X.16 + X.17 + X.18 + X.19 + X.20 +
                        X.21 + X.22 + X.23 + X.24 + X.25, id = ID, data = dat, ensemble = FALSE)
    fit2 <- rocTree(Surv(Time, death) ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7 + X.8 + X.9 + X.10 +
                        X.11 + X.12 + X.13 + X.14 + X.15 + X.16 + X.17 + X.18 + X.19 + X.20 +
                        X.21 + X.22 + X.23 + X.24 + X.25, id = ID, data = dat)
    err1 <- err2 <- rep(NA, 200)
    for (i in 1:200) {
        if (sce == 1) newX <- t(matrix(mvrnorm(1, rep(0,25), 0.9^outer(1:25, 1:25, function(x, y)
            abs(x - y))), 25, n))
        if (sce == 2) newX <- t(matrix(runif(25), 25, n))
        if (sce == 3) newX <- t(matrix(mvrnorm(1, rep(0,25), 0.75^outer(1:25, 1:25, function(x, y)
            abs(x - y))), 25, n))
        if (sce == 4) newX <- t(matrix(mvrnorm(1, rep(0,25), 0.75^outer(1:25, 1:25, function(x, y)
            abs(x - y))), 25, n))
        newX <- as.data.frame(newX)
        names(newX) <- names(dat)[grep("X.", names(dat))]
        newX$Time <- sort(unique(dat$Time))
        if (sce == 1) trueSv <- 1 - pexp(.t0, exp(-.1 * sum(newX[1, 11:25])))
        if (sce == 2) trueSv <- 1 - pexp(.t0, 1 / (sin(newX[1, 1] * pi) +
                                                   2 * abs(newX[1, 2] - .5) + newX[1, 3]^3))
        if (sce == 3) trueSv <- 1 - pgamma(.t0, shape = .5 + .3 * abs(sum(newX[1, 11:15])), scale = 2)
        if (sce == 4) trueSv <- 1 - plnorm(.t0, meanlog = .1 * sum(newX[1, c(1:5, 21:25)]))        
        ref <- stepfun(.t0, c(1, trueSv))(.t0)
        pred1 <- predict(fit1, newX)
        pred2 <- predict(fit2, newX)
        err1[i] <- mean(abs(with(pred1$pred, stepfun(Time, c(1, Survival))(.t0) - ref)))
        err2[i] <- mean(abs(with(pred2$pred, stepfun(Time, c(1, Survival))(.t0) - ref)))
    }
    return(c(mean(err1), mean(err2)))
}

set.seed(1)
system.time(print(do(100, 0, 1)))
## [1] 0.09775445 0.06651857
##    user  system elapsed 
##   8.084   0.012   8.094 

set.seed(1)
system.time(print(do(500, 0, 1)))

library(parallel)

cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do"))
invisible(clusterExport(NULL, "indCov1"))
invisible(clusterExport(NULL, "indCov2"))
invisible(clusterExport(NULL, "indCov3"))
invisible(clusterExport(NULL, "ord"))
invisible(clusterExport(NULL, "trueSurv"))
invisible(clusterExport(NULL, "ti2td"))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(MASS)))
invisible(clusterEvalQ(NULL, library(rpart)))
invisible(clusterEvalQ(NULL, library(ranger)))
invisible(clusterEvalQ(NULL, library(randomForestSRC)))

sim1.100.0 <- parSapply(NULL, 1:50, function(z) do(500, 0, 1))
sim2.100.0 <- parSapply(NULL, 1:100, function(z) do(500, 0, 2))
sim3.100.0 <- parSapply(NULL, 1:100, function(z) do(500, 0, 3))

stopCluster(cl)

## With transformation
rowMeans(sim1.100.0) # 0.09120784 0.06411867
rowMeans(sim2.100.0) # 0.10937177 0.08108885
rowMeans(sim3.100.0) # 0.1380290 0.1176399
