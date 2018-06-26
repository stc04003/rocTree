#######################################################################
## Simulation codes
#######################################################################

## --------------------------------------------------------------------
## Simulation 1 with time independent covariates
## --------------------------------------------------------------------

## \lambda(t, Z) = \lambda_0(t) exp(2Z_1 + 2Z_2), \lambda_0(t) = 2t
sim1.1 <- function(n, cen = 0) {
    if (!(cen %in% c(0, .25, .50)))
        stop("Only 3 levels of censoring rates (0%, 25%, 50%) are allowed.")
    z1 <- runif(n)
    z2 <- runif(n)
    Time <- sqrt(rexp(n) * exp(-2 * z1 - 2 * z2))
    if (cen == 0) cens <- Inf
    if (cen == .25) cens <- runif(n, 0, 1.41)
    if (cen == .50) cens <- runif(n, 0, 0.66)
    data.frame(Y = pmin(Time, cens), delta = 1 * (Time <= cens), z1 = z1, z2 = z2)
}

## \lambda(t, Z) = \lambda_0(t) exp\{2sin(2\pi Z_1) + 2|Z_2 - 0.5|\}, \lambda_0(t) = 2t
sim1.2 <- function(n, cen = 0) {
    if (!(cen %in% c(0, .25, .50)))
        stop("Only 3 levels of censoring rates (0%, 25%, 50%) are allowed.")
    z1 <- runif(n)
    z2 <- runif(n)
    Time <- sqrt(rexp(n) * exp(-2 * sin(2 * pi * z1) - 2 * abs(z2 - .5)))
    if (cen == 0) cens <- Inf
    if (cen == .25) cens <- runif(n, 0, 1.23)
    if (cen == .50) cens <- runif(n, 0, 0.59)
    data.frame(Y = pmin(Time, cens), delta = 1 * (Time <= cens), z1 = z1, z2 = z2)    
}

## same as sim1.1 with 3 noise variables
sim1.5 <- function(n, cen = 0) {
    if (!(cen %in% c(0, .25, .50)))
        stop("Only 3 levels of censoring rates (0%, 25%, 50%) are allowed.")
    z1 <- runif(n)
    z2 <- runif(n)
    Time <- sqrt(rexp(n) * exp(-2 * z1 - 2 * z2))
    if (cen == 0) cens <- Inf
    if (cen == .25) cens <- runif(n, 0, 1.41)
    if (cen == .50) cens <- runif(n, 0, 0.66)
    data.frame(Y = pmin(Time, cens), delta = 1 * (Time <= cens),
               z1 = z1, z2 = z2, z3 = runif(n), z4 = runif(n), z5 = runif(n))
}

## AFT model: log(T) = -2 + 2z1 + 2z2 + e, e ~ N(0, 0.5^2)
sim1.3 <- function(n, cen = 0) {
    if (!(cen %in% c(0, .25, .50)))
        stop("Only 3 levels of censoring rates (0%, 25%, 50%) are allowed.")
    z1 <- runif(n)
    z2 <- runif(n)
    Time <- exp(-2 + 2 * z1 + 2 * z2 + rnorm(n, sd = .5))
    cens <- runif(n, 0, cen)
    if (cen == 0) cens <- Inf
    if (cen == .25) cens <- runif(n, 0, 6.00)
    if (cen == .50) cens <- runif(n, 0, 2.40)
    data.frame(Y = pmin(Time, cens), delta = 1 * (Time <= cens), z1 = z1, z2 = z2)    
}

## Generalized gamma family: T = exp(\sigma \omega), \omega = log(Q^2g)/Q, g~Gamma(Q^-2 1), \sigma = 2*Z1, Q = 2*Z2
sim1.4 <- function(n, cen = 0) {
    if (!(cen %in% c(0, .25, .50)))
        stop("Only 3 levels of censoring rates (0%, 25%, 50%) are allowed.")
    z1 <- runif(n)
    z2 <- runif(n)    
    sig <- 2 * z1
    Q <- 2 * z2
    g <- rgamma(n, Q^-2, 1)
    w <- log(Q^2 * g) / Q
    Time <- exp(sig * w)
    cens <- runif(n, 0, cen)
    if (cen == 0) cens <- Inf
    if (cen == .25) cens <- runif(n, 0, 4.12)
    if (cen == .50) cens <- runif(n, 0, 1.63)
    data.frame(Y = pmin(Time, cens), delta = 1 * (Time <= cens), z1 = z1, z2 = z2)    
}

## --------------------------------------------------------------------
## Simulation 2 with time dependent covariates
## \lambda(t, Z(t)) = \lambda_0(t) exp(Z_1(t) + Z2)
## --------------------------------------------------------------------

## Dichotomous time dependent covariate with at most one change in value:
## Z1(t) = E * I(t\ge U_0) + (1 - E) * I(t < U_0), U_0 ~ rexp(5), E ~ Bernoulli(.5)
sim2.1 <- function(n, cen = 0) {

}

## Dichotomous time dependent covariate with multiple changes:
## Z_1(t) = E * (I(U_1\le t < U_2) + I(U_3\le t)) + (1 - E) * (I(t < U_1) + I(U_2\le t < U_3));
## U1, U2, U3 ~ exp(10), E ~ Bernoulli(.5)
sim2.2 <- function(n, cen = 0) {

}

## Continuous time dependent covaraite:
## Z_1(t) = kt + b, k ~ runif(1, 2), b ~ runif(1, 2)
sim2.3 <- function(n, cen = 0) {

}

##########################################################################
## Starting to put things in package form
## Scenario III on page 18
datGen <- function(n) {
    ## 3 covariates, X1 and x2 are used in the paper, X3 is a noise
    lambda <- .1
    P <- 2
    U <- runif(n)
    X <- matrix(runif(n * P, 0, 1), ncol = P, nrow = n)
    a <- X[, 1]
    b <- X[, 2]
    betat <- 2
    beta <- 1
    ks <- runif(n, 1, 2)
    it <- runif(n, 1, 2)
    T <- log(1 + (betat * ks * (-log(U))) / lambda / exp(beta * a + betat * it)) / betat / ks
    C <- runif(n, 0, .8)
    Y <- pmin(T, C)
    E <- T <= C ## censoring indicator
    dat <- NULL ## ID, Time, Status, X1, X2, X3
    for (i in 1:n) {
        dat <- rbind(dat, cbind(rep(i, rank(Y)[i]), Y[which(Y <= Y[i])], E[i], a[i], Y[which(Y <= Y[i])] * ks[i] + it[i], runif(1)))
    }
    dat <- data.frame(dat)
    colnames(dat) <- c("ID", "Time", "Status", "X1", "X2", "X3")
    dat <- dat[order(dat$ID, dat$Time),]
    rownames(dat) <- NULL
    dat
}

set.seed(1234)
dat <- datGen(100)
