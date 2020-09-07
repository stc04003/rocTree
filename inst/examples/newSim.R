## #############################################################################
## Simulation settings in Steingrimsson et al (2018).
## Section 5.1
## Time-independent covariates:
## #############################################################################

#' Setting 1 in Steingrimsson et al (2018).
#'
#' 25 covariates from multivariate normal with mean 0 and covariance matrix \eqn{0.9^|i - j|}
#' This is proportional hazard model.
#' censoring from exponential distribution
#'
#'
#' @param n is sample size
#' @param cen is censoring percentage (0\%, 25\%, 50\% are allowed)
#' @param tdForm to print the data in time dependent form?
#' @param order to order data by survival time (ascending order)
#'
#' @importFrom MASS mvrnorm
#' @export
#'
indCov1 <- function(n, cen, tdForm = TRUE, order = TRUE) {
    W <- matrix(mvrnorm(n, rep(0, 25), .9^outer(1:25, 1:25, function(x, y) abs(x - y))), n)
    Y0 <- rexp(n, exp(-.1 * rowSums(W[,11:25])))
    if (cen == 0) cc <- Inf
    if (cen == .25) cc <- rexp(n, .25)
    if (cen == .50) cc <- rexp(n, 1)
    death <- 1 * (Y0 <= cc)
    Y <- pmin(Y0, cc)
    dat <- data.frame(ID = 1:n, Time = Y, death = death, X = W)
    if (order) dat <- ord(dat)
    if (tdForm) dat <- ti2td(dat)
    attr(dat, "from") <- "indCov1"
    return(dat)
}

#' Computes the empirical survival curve given the covariate matrix
#'
#' Gives survival curve
#'
#' @param dat is the data.frame generated from indCov1~4
#' @param n is the number of replicates used in computing ecdf
#'
#' @export
trueSurv <- function(dat, n = 5e3) {
    .case <- attr(dat, "from")
    if (!is.null(dat$death)) dat <- subset(dat, death > 0)
    .X <- dat[,grep("X.", names(dat))]
    ## .Y <- dat$Time
    if (.case == "indCov1")
        .tt <- replicate(n, rexp(nrow(.X), exp(-.1 * rowSums(.X[,11:25]))))
    if (.case == "indCov2")
        .tt <- replicate(n, rexp(nrow(.X), 1 / (sin(.X[,1] * pi) + 2 * abs(.X[,2] - .5) + .X[,3]^3)))
    if (.case == "indCov3")
        .tt <- replicate(n, rgamma(nrow(.X),
                                   shape = .5 + .3 * abs(rowSums(.X[,11:15])),
                                   scale = 2))
    if (.case == "indCov4")
        .tt <- replicate(n, rlnorm(n, meanlog = .1 * (rowSums(.X[,1:5]) + rowSums(.X[,21:25]))))
    .tt <- rowMeans(.tt)
    list(Time = .tt[order(.tt)], Surv = 1 - ecdf(.tt)(.tt[order(.tt)]))
    ## if (.case == "indCov1") .sv <- 1 - pexp(.Y, exp(-.1 * rowSums(.X[,11:25])))
    ## if (.case == "indCov2") .sv <- 1 - pexp(.Y, 1 / (sin(.X[,1] * pi) + 2 * abs(.X[,2] - .5) + .X[, 3]^3))
    ## if (.case == "indCov3") .sv <- 1 - pgamma(.Y, shape = .5 + .3 * abs(rowSums(.X[,11:15])), scale = 2)
    ## if (.case == "indCov4") .sv <- 1 - plnorm(.Y, meanlog = .1 * (rowSums(.X[,1:5]) + rowSums(.X[,21:25])))
}

#' Give .z0
#'


#' Setting 2 in Steingrimsson et al (2018).
#'
#' 25 covariates from Uniform(0, 1)
#' Proportional hazards assumption is mildly violated
#' censoring from Uniform(0, 6)
#'
#' @param n is sample size
#' @param cen is censoring percentage (0\%, 25\%, 50\% are allowed)
#' @param tdForm to print the data in time dependent form?
#' @param order to order data by survival time (ascending order)
#'
#' @export
indCov2 <- function(n, cen, tdForm = TRUE, order = TRUE) {
    W <- matrix(runif(25 * n), ncol = 25)
    Y0 <- rexp(n, 1 / (sin(W[,1] * pi) + 2 * abs(W[,2] - .5) + W[,3]^3))
    cc <- runif(n, 0, cen)
    if (cen == 0) cc <- Inf
    if (cen == .25) cc <- runif(n, 0, 5.3)
    if (cen == .50) cc <- runif(n, 0, 2)
    death <- 1 * (Y0 <= cc)
    Y <- pmin(Y0, cc)
    dat <- data.frame(ID = 1:n, Time = Y, death = death, X = W)
    if (order) dat <- ord(dat)
    if (tdForm) dat <- ti2td(dat)
    attr(dat, "from") <- "indCov2"
    return(dat)
}


#' Setting 3 in Steingrimsson et al (2018).
#'
#' 25 covariates from multivariate normal with mean 0 and covariance matrix \eqn{0.75^|i - j|}
#' Survival times are gamma distributed with shape = 0.5 + 0.3 * abs(W_11, ..., W_15) and scale 2
#' Censoring time is uniform
#'
#' @param n is sample size
#' @param cen is censoring percentage (0\%, 25\%, 50\% are allowed)
#' @param tdForm to print the data in time dependent form?
#' @param order to order data by survival time (ascending order)
#'
#' @export
indCov3 <- function(n, cen, tdForm = TRUE, order = TRUE) {
    W <- matrix(mvrnorm(n, rep(0, 25), .75^outer(1:25, 1:25, function(x, y) abs(x - y))), n)
    Y0 <- rgamma(n, shape = 0.5 + 0.3 * abs(rowSums(W[,11:15])), scale = 2)
    if (cen == 0) cc <- Inf
    if (cen == .25) cc <- runif(n, 0, 11.7)
    if (cen == .50) cc <- runif(n, 0, 4.9)
    death <- 1 * (Y0 <= cc)
    Y <- pmin(Y0, cc)
    dat <- data.frame(ID = 1:n, Time = Y, death = death, X = W)
    if (order) dat <- ord(dat)
    if (tdForm) dat <- ti2td(dat)
    attr(dat, "from") <- "indCov3"
    return(dat)
}


#' Setting 4 in Steingrimsson et al (2018).
#'
#' 25 covariates from multivariate normal with mean 0 and covariance matrix \eqn{0.75^|i - j|}
#' Survival times are log-normal distributed with mean = 0.1|W_1 + ... + W_5 + W_21 + ... W_25.
#' shape = 0.5 + 0.3 * abs(W_11, ..., W_15) and scale 2
#' Censoring time is uniform
#'
#' @param n is sample size
#' @param cen is censoring percentage (0\%, 25\%, 50\% are allowed)
#' @param tdForm to print the data in time dependent form?
#' @param order to order data by survival time (ascending order)
#'
#' @export
indCov4 <- function(n, cen, tdForm = TRUE, order = TRUE) {
    W <- matrix(mvrnorm(n, rep(0, 25), .75^outer(1:25, 1:25, function(x, y) abs(x - y))), n)
    Y0 <- rlnorm(n, meanlog = 0.1 * (rowSums(W[,1:5]) + rowSums(W[,21:25])))
    if (cen == 0) cc <- Inf
    if (cen == .25) cc <- rlnorm(n, meanlog = 0.1 * (rowSums(W[,1:5]) + rowSums(W[,21:25])) + 1)
    if (cen == .50) cc <- rlnorm(n, meanlog = 0.1 * (rowSums(W[,1:5]) + rowSums(W[,21:25])) + .01)
    death <- 1 * (Y0 <= cc)
    Y <- pmin(Y0, cc)
    dat <- data.frame(ID = 1:n, Time = Y, death = death, X = W)
    if (order) dat <- ord(dat)
    if (tdForm) dat <- ti2td(dat)
    attr(dat, "from") <- "indCov4"
    return(dat)
}


#' @noRd
#' @keywords internal
#' Internal function to convert time independent data to time dependent data
#'
#' @param dat data.frame prepared from one of the `indCov` function.
ti2td <- function(dat) {
    dat2 <- data.frame(ID = dat$ID,
                       Time = c(t(dat$Time * ifelse(outer(dat$Time, dat$Time, "<="), TRUE, NA))),
                       death = c(diag(dat$death)),
                       ## death = c(t(dat$death * ifelse(outer(dat$ID, dat$ID, "=="), 1, 0))),
                       dat[,grep("X.", names(dat))])
    dat2 <- dat2[complete.cases(dat2),]
    dat <- dat2[order(dat2$ID),]
    rownames(dat) <- NULL
    return(dat)
}

#' @noRd
#' @keywords internal
#' Internal function to order the survival times from small to large
#'
#' @param dat data.frame prepared from one of the `indCov` function.
ord <- function(dat) {
    dat <- dat[order(dat$ID, dat$Time), ]
    ord0 <- aggregate(Time ~ ID, dat, max)$Time
    dat <- dat[order(rep(ord0, aggregate(Time ~ ID, dat, length)$Time)),]
    dat$ID <- rep(1:length(unique(dat$ID)), table(dat$ID)[unique(dat$ID)])
    rownames(dat) <- NULL
    return(dat)
}

## #############################################################################
## Time-dependent covariates
## #############################################################################

#' Time-dependent covariate case; modified form the original `sim2.3` of the first `rocTree` pkg.
#' Continuous time-dependent with Cox model
#' This gives 1 time-independent covariate (`X.1`) and 10 time-dependent covariate (`Z.1` to `Z.10`)
#'
#'
#' @param n is the sample size
#' @param cen is censoring percentage (0\%, 25\%, 50\% are allowed)
#' @param order to order data by survival time (ascending order)
#'
#' @export
dCov1 <- function (n, cen, order = TRUE) {
    k <- runif(n, 1, 2)
    b <- runif(n, 1, 2)
    z2 <- runif(n)
    Time <- log(10 * rexp(n) * exp(-z2 - b) * k + 1)/k
    if (cen == 0)
        cens <- rep(Inf, n)
    if (cen == 0.25)
        cens <- runif(n, 0, 2.48)
    if (cen == 0.5)
        cens <- runif(n, 0, 1.23)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(ID = x,
                   Time = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   X.1 = z2[x], k = k[x], b = b[x])))
    dat <- dat[order(dat$ID, dat$Time), ]
    tmp <- t(apply(dat, 1, function(zz)
        mvrnorm(1, rep(zz[2] * zz[5] + zz[6], 10),
                .5 * .9^outer(1:10, 1:10, function(x, y) abs(x - y)))))
    colnames(tmp) <- paste("Z", 1:10, sep = ".")
    dat <- cbind(dat, tmp)
    if (order) dat <- ordDat(dat)
    attr(dat, "from") <- "dCov1"
    return(dat)
}

ordDat <- function(dat) {
    dat <- dat[order(dat$ID, dat$Time), ]
    ord0 <- aggregate(Time ~ ID, dat, max)$Time
    dat <- dat[order(rep(ord0, aggregate(Time ~ ID, dat, length)$Time)), ]
    dat$ID <- rep(1:length(unique(dat$ID)), table(dat$ID)[unique(dat$ID)])
    rownames(dat) <- NULL
    return(dat)
}


#' Time-dependent covariate case; modified form the original `sim2.4` of the first `rocTree` pkg.
#' Continuous time-dependent with non-Cox model
#' This gives 1 time-independent covariate (`X.1`) and 10 time-dependent covariate (`Z.1` to `Z.10`)
#'
#' @param n is the sample size
#' @param cen is censoring percentage (0\%, 25\%, 50\% are allowed)
#' @param order to order data by survival time (ascending order)
#'
#' @export
dCov2 <- function(n, cen, order = TRUE) {
    k <- runif(n, 1, 2)
    b <- runif(n, 1, 2)
    z2 <- runif(n)
    u <- runif(n)
    invF <- function(x, k, b, z2, u) {
        ## all assume to be 1 dimensional
        10 * log(u) + x - cos(k * x + b + z2) / k + cos(b + z2) / k
    }
    Time <- sapply(1:n, function(y)
        uniroot(f = invF, interval = c(0, 500), k = k[y], b = b[y], z2 = z2[y], u = u[y])$root)
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 39.5)
    if (cen == .50) cens <- runif(n, 0, 16.5)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(ID = x, Time = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   X.1 = z2[x], k = k[x], b = b[x])))
    dat <- dat[order(dat$ID, dat$Time),]
    tmp <- t(apply(dat, 1, function(zz)
        mvrnorm(1, rep(zz[2] * zz[5] + zz[6], 10),
                .5 * .9^outer(1:10, 1:10, function(x, y) abs(x - y)))))
    colnames(tmp) <- paste("Z", 1:10, sep = ".")
    dat <- cbind(dat, tmp)
    if (order) dat <- ord(dat)
    attr(dat, "from") <- "dCov2"
    return(dat)
}

#' Time-dependent covariate case; modified form the original `sim2.5` of the first `rocTree` pkg.
#' Continuous time-dependent with non-Cox model and non-monotonic
#' This gives 1 time-independent covariate (`X.1`) and 10 time-dependent covariate (`Z.1` to `Z.10`)
#'
#' @param n is the sample size
#' @param cen is censoring percentage (0\%, 25\%, 50\% are allowed)
#' @param order to order data by survival time (ascending order)
#'
#' @export
dCov3 <- function(n, cen, order = TRUE) {
    k <- runif(n, 1, 2)
    b <- runif(n, 1, 2)
    z2 <- runif(n)
    u <- runif(n)
    invF <- function(x, k, b, z2, u) {
        ## all assume to be 1 dimensional
        10 * log(u) + x - cos(k * x * (2 * (x > 5) - 1) + b + z2) / k + cos(b + z2) / k
    }
    Time <- sapply(1:n, function(y)
        uniroot(f = invF, interval = c(0, 500), k = k[y], b = b[y], z2 = z2[y], u = u[y])$root)
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 42)
    if (cen == .50) cens <- runif(n, 0, 17)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(ID = x, Time = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   X.1 = z2[x], k = k[x], b = b[x])))
    dat <- dat[order(dat$ID, dat$Time),]
    tmp <- t(apply(dat, 1, function(zz)
        mvrnorm(1, rep(zz[2] * zz[5] * (2 * (zz[2] > 5) - 1)+ zz[6], 10),
                .5 * .9^outer(1:10, 1:10, function(x, y) abs(x - y)))))
    colnames(tmp) <- paste("Z", 1:10, sep = ".")
    dat <- cbind(dat, tmp)
    if (order) dat <- ord(dat)
    attr(dat, "from") <- "dCov3"
    return(dat)
}


#' Time-dependent covariate case; 
#' Continuous time-dependent with non-Cox model and non-monotonic
#' This gives 10 time-independent covariate (`X.1` to `X.10`) and
#' 10 time-dependent covariate (`Z.1` to `Z.10`).
#'
#' The hazard function is \eqn{\sum_{i = 1}^{10}(Z_i(t) - X_i)^2} and
#' \eqn{Z_i(t) = k_i t}, where \eqn{k_i} is independent unif(0, 1)
#'
#' @param n is the sample size
#' @param cen is censoring percentage (0\%, 25\%, 50\% are allowed)
#' @param order to order data by survival time (ascending order)
#'
#' @export
dCov4 <- function(n, cen, order = TRUE) {
    k <- matrix(runif(10 * n), n)
    X <- matrix(runif(10 * n), n)
    u <- rexp(n)
    invF <- function(x, k2, x2, kx, u) {
        ## x is survival times
        ## k2 is sum(k^2)
        ## x2 is sum(x^2)
        ## kx is sum(k * x)
        ## u is rexp r.v.
        x <- x / 10
        (x^3 / 3) * k2 - (x^2) * kx + x * x2 - u
    }
    Time <- sapply(1:n, function(y)
        uniroot(f = invF, interval = c(0, 500), k2 = sum(k[y,]^2),
                x2 = sum(X[y,]^2), kx = sum(k[y,] * X[y,]), u = u[y])$root)
    if (cen == 0) 
        cens <- rep(Inf, n)
    if (cen == 0.25) 
        cens <- runif(n, 0, 18.63)
    if (cen == 0.5) 
        cens <- runif(n, 0, 6.6)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(ID = x, 
                   Time = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])))))
    dat <- dat[order(dat$ID, dat$Time), ]
    dat <- cbind(dat, X = X[dat$ID, ], k = k[dat$ID,])
    dat <- cbind(dat, Z = sapply(1:10, function(i)
        eval(parse(text = paste("with(dat, Z.", i, " <- k.", i, " * Time)", sep = "")))))
    if (order) dat <- ord(dat)
    attr(dat, "from") <- "dCov4"
    return(dat)
}

#' Time-dependent covariate case; 
#' Continuous time-dependent with non-Cox model and non-monotonic
#' This gives 10 time-independent covariate (`X.1` to `X.10`) and
#' 10 time-dependent covariate (`Z.1` to `Z.10`).
#'
#' The hazard function is \eqn{\sum_{i = 1}^{10}(Z_i(t) - X_i)^2} and
#' \eqn{Z_i(t) = k_i t}, where \eqn{k_i} is multivariate normal with
#' exchangeable covariance variance with pho = 0.9
#' 
#' @param n is the sample size
#' @param cen is censoring percentage (0\%, 25\%, 50\% are allowed)
#' @param order to order data by survival time (ascending order)
#'
#' @export
dCov5 <- function(n, cen, order = TRUE) {
    k <- mvrnorm(n, rep(1, 10), .9^outer(1:10, 1:10, function(a, b) abs(a - b)))
    X <- matrix(rnorm(10 * n), n)
    u <- rexp(n)
    invF <- function(x, k2, x2, kx, u) {
        x <- x / 20
        (x^3 / 3) * k2 - (x^2) * kx + x * x2 - u
    }
    Time <- sapply(1:n, function(y)
        uniroot(f = invF, interval = c(0, 500), k2 = sum(k[y,]^2),
                x2 = sum(X[y,]^2), kx = sum(k[y,] * X[y,]), u = u[y])$root)
    if (cen == 0) 
        cens <- rep(Inf, n)
    if (cen == 0.25) 
        cens <- runif(n, 0, 8.25)
    if (cen == 0.5) 
        cens <- runif(n, 0, 3.34)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(ID = x, 
                   Time = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])))))
    dat <- dat[order(dat$ID, dat$Time), ]
    dat <- cbind(dat, X = X[dat$ID, ], k = k[dat$ID,])
    dat <- cbind(dat, Z = sapply(1:10, function(i)
        eval(parse(text = paste("with(dat, Z.", i, " <- k.", i, " * Time)", sep = "")))))
    if (order) dat <- ord(dat)
    attr(dat, "from") <- "dCov5"
    return(dat)
}
