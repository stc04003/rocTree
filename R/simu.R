globalVariables(c("n", "cen", "Y", "id")) ## global variables for simu

#' Function to generate simulated data used in the manuscript.
#'
#' This function is used to generate simulated data in the manuscript.
#' The underlying model is assumed to be
#' \deqn{\lambda(t, Z) = \lambda_0(t) e^{(\beta_1 * Z_1 + \beta_2 * Z_2)},}
#' where \eqn{\lambda_0(t)} is the baseline hazard function,
#' (\eqn{b_1, b_2}) is the regression coefficient, and
#' \eqn{Z_1} and \eqn{Z_2} are the covariates.
#' When time varying covariate is involved (scenario 2), only \eqn{Z_1} is assumed to dependent on time.
#'
#' The following scenarios are included.
#' \describe{
#' \item{Scenario 1}{assumes time independent covariates with \eqn{\Lambda(0) = 2t}:}
#' \describe{
#' \item{1.1}{Proportional hazards model:
#' \eqn{\lambda(t, Z) = \lambda_0(t) e^{(2Z_1 + 2Z_2)}}.}
#' \item{1.2}{Proportional hazards model with nonlinear covariate effects:
#' \eqn{\lambda(t, Z) = \lambda_0(t) e^{[2\sin(2\pi Z_1) + 2|Z_2 - 0.5|]}}.}
#' \item{1.3}{Accelerated failure time model:
#' \eqn{\log(T) = -2 + 2Z_1 + 2Z_2 + \epsilon}, where \eqn{\epsilon} follows \eqn{N(0, 0.5^2)}.}
#' \item{1.4}{Generalized gamma family:
#' \eqn{T = e^{\sigma\omega}} with \eqn{\omega = \log(Q^2 g) / Q}, \eqn{g} follows Gamma(\eqn{Q^{-2}, 1}),
#' \eqn{\sigma = 2Z_1, Q = 2Z_2}.}
#' \item{1.5}{Proportional hazards model with noise variable:
#' \eqn{\lambda(t, Z) = \lambda_0(t) e^{(2Z_1 + 2Z_2 + 0Z_3 + 0Z_4 + 0Z_5)}}.}
#' \item{1.6}{Proportional hazards model with 10 covaraites:
#' \eqn{\lambda(t, Z) = \lambda_0(t) e^{\sum_{j = 1}^pZ_j}}, for \eqn{\eta_j = 0.1, j = 1, ..., 10.}}
#' \item{1.7}{Proportional hazards model with noise variable:
#' \eqn{\lambda(t, Z) = \lambda_0(t) e^{(2Z_1 + 2Z_2 + 0Z_3 + 0Z_4 + ... + 0Z_{10})}}.}
#' }
#' \item{Scenario 2}{assumes time dependent covariate (\eqn{Z_1}) with \eqn{\Lambda(0) = 2t}.
#' The survival times are generated from the hazard \eqn{\lambda(t, Z(t)) = \lambda_0(t)e^{Z_1(t) + Z_2}},
#' where \eqn{Z_1(t)} is the time-dependent covariate.}
#' \describe{
#' \item{2.1}{assumes}
#' \item{2.2}{assumes}
#' }
#' }
#' 
#' @param n is the number of subject
#' @param cen is the censoring percentage; right now it can be either 0\%, 25\%, or 50\%.
#' @param scenario can be numeric or character string.
#' This indicates the simulation scenario noted in the manuscript.
#' See \bold{Details} for all options.
#' @param summary a logical value indicating whether a brief data summary will be printed.
#' 
#' @importFrom stats delete.response rexp rgamma rnorm runif rbinom uniroot
#' @importFrom tibble as.tibble
#' @importFrom dplyr "%>%" arrange select
#' 
#' @return \code{simu} returns a \code{data.frame} in the class of "roc.simu".
#' This is needed for \code{trueHaz} and \code{trueSurv}.
#' The returned data.frame consists of columns:
#' \describe{
#' \item{Y}{is the observed follow-up time.}
#' \item{death}{is the death indicator; death = 0 if censored.}
#' \item{z1}{is the time-independent covariate.}
#' \item{z2}{is the secondary covariate. It is time-independent in Scenario 1 and time-dependent in Scenario 2 and 3.}
#' }
#'
#' @name simu
#' @rdname simu
#' @export
#' 
simu <- function(n, cen, scenario, summary = FALSE) {
    if (!(cen %in% c(0, .25, .50)))
        stop("Only 3 levels of censoring rates (0%, 25%, 50%) are allowed.")
    ## if (!(scenario %in% paste(rep(1:3, each = 6), rep(1:7, 3), sep = ".")))
    ##     stop("See ?simu for scenario definition.")
    if (n > 1e4) stop("Sample size too large.")
    dat <- as.tibble(eval(parse(text = paste("sim", scenario, "(n = ", n, ", cen = ", cen, ")", sep = ""))))
    if (summary) {
        cat("\n")
        cat("Summary results:\n")
        cat("Number of subjects:", n)
        cat("\nNumber of subjects experienced death:", sum(dat$death) / n)
        if (substr(scenario, 1, 1) == "1") {
            cat("\nTime independent covaraites: z1 and z2.")
        } else {
            cat("\nTime independent covaraites: z1.")
            cat("\nTime dependent covaraites: z2.")
        }
        cat("\n\n")
    }
    attr(dat, "scenario") <- scenario
    attr(dat, "prepBy") <- "rocSimu"
    return(dat %>% arrange(id, Y))
}

#' Function to generate the true hazard used in the simulation.
#'
#' This function is used to generate the true cumulative hazard function for ONE ID used in the simulation.
#'
#' @param dat is a data.frame prepared by \code{simu}.
#'
#' @importFrom stats approxfun complete.cases
#' @importFrom utils head
#' @rdname simu
#' @export
trueHaz <- function(dat) {
    if (attr(dat, "prepBy") != "rocSimu") stop("Inputed data must be prepared by \"simu\".")
    scenario <- attr(dat, "scenario")
    dat <- dat %>% arrange(Y)
    cumHaz <- eval(parse(text = paste("trueHaz", scenario, "(dat)", sep = "")))
    approxfun(x = dat$Y, y = cumHaz, method = "constant", yleft = 0, yright = max(cumHaz))
}

#' Function to generate the true survival used in the simulation.
#'
#' This function is used to generate the true survival function for ONE ID used in the simulation.
#'
#' @rdname simu
#' @export
#' 
trueSurv <- function(dat) {
    if (attr(dat, "prepBy") != "rocSimu") stop("Inputed data must be prepared by \"simu\".")
    scenario <- attr(dat, "scenario")
    dat <- dat %>% arrange(Y)
    Surv <- eval(parse(text = paste("trueSurv", scenario, "(dat)", sep = "")))
    approxfun(x = dat$Y, y = Surv, method = "constant", yleft = 1, yright = min(Surv))
}

#' Function to generate testing sets given training data.
#'
#' This function is used to generate testing sets for each scenario.
#'
#' @rdname simu
#' @export
simuTest <- function(dat) {
    if (attr(dat, "prepBy") != "rocSimu") stop("Inputed data must be prepared by \"simu\".")
    scenario <- attr(dat, "scenario")
    dat <- as.tibble(eval(parse(text = paste("simuTest", scenario, "(dat)", sep = ""))))
    attr(dat, "prepBy") <- "rocSimu"
    attr(dat, "scenario") <- scenario
    return(dat)
}

#' ##########################################################################################
#' ##########################################################################################
#' Background functions for simulation
#' @keywords internal
#' @noRd
sim1.1 <- function(n, cen = 0) {
    z1 <- runif(n)
    z2 <- runif(n)
    Time <- sqrt(rexp(n) * exp(-2 * z1 - 2 * z2))
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 1.41)
    if (cen == .50) cens <- runif(n, 0, 0.66)
    Y <- pmin(Time, cens)
    d <- 1 * (Time <= cens)
    do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), d[x]),  z1 = z1[x], z2 = z2[x])))
}

sim1.2 <- function(n, cen = 0) {
    z1 <- runif(n)
    z2 <- runif(n)
    Time <- sqrt(rexp(n) * exp(-2 * sin(2 * pi * z1) - 2 * abs(z2 - .5)))
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 3.48) ## 1.23)
    if (cen == .50) cens <- runif(n, 0, 1.50) ## 0.59)
    Y <- pmin(Time, cens)
    do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])), z1 = z1[x], z2 = z2[x])))
}

sim1.5 <- function(n, cen = 0) {
    z1 <- runif(n)
    z2 <- runif(n)
    z3 <- runif(n)
    z4 <- runif(n)
    z5 <- runif(n)
    Time <- sqrt(rexp(n) * exp(-2 * z1 - 2 * z2))
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 1.41)
    if (cen == .50) cens <- runif(n, 0, 0.66)
    Y <- pmin(Time, cens)
    do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   z1 = z1[x], z2 = z2[x], z3 = z3[x], z4 = z4[x], z5 = z5[x])))
}

sim1.3 <- function(n, cen = 0) {
    z1 <- runif(n)
    z2 <- runif(n)
    Time <- exp(-2 + 2 * z1 + 2 * z2 + rnorm(n, sd = .5))
    cens <- runif(n, 0, cen)
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 6.00)
    if (cen == .50) cens <- runif(n, 0, 2.40)
    Y <- pmin(Time, cens)
    do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])), z1 = z1[x], z2 = z2[x])))
}

sim1.4 <- function(n, cen = 0) {
    z1 <- runif(n)
    z2 <- runif(n)    
    sig <- 2 * z1
    Q <- 2 * z2
    g <- rgamma(n, Q^-2, 1)
    w <- log(Q^2 * g) / Q
    Time <- exp(sig * w)
    cens <- runif(n, 0, cen)
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 4.12)
    if (cen == .50) cens <- runif(n, 0, 1.63)
    Y <- pmin(Time, cens)
    do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])), z1 = z1[x], z2 = z2[x])))
}

sim1.6 <- function(n, cen = 0) {
    z <- matrix(runif(n * 10), n, 10)
    ## z <- matrix(rnorm(n * 10), n, 10)
    Time <- sqrt(rexp(n) * exp(-rowSums(z[,1:2]) * 2))
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 2.76)
    if (cen == .50) cens <- runif(n, 0, 1.36)
    Y <- pmin(Time, cens)
    d <- 1 * (Time <= cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), d[x]))))
    dat <- cbind(dat, z[dat$id,])
    names(dat)[4:13] <- paste("z", 1:10, sep = "")
    dat
}

sim1.7 <- function(n, cen = 0) {
    z <- matrix(runif(n * 10), n, 10)
    Time <- sqrt(rexp(n) * exp(-rowSums(z[,1:2]) * 2))
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 1.41)
    if (cen == .50) cens <- runif(n, 0, 0.66)
    Y <- pmin(Time, cens)
    d <- 1 * (Time <= cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), d[x]))))
    dat <- cbind(dat, z[dat$id,])
    names(dat)[4:13] <- paste("z", 1:10, sep = "")
    dat    
}

#' @importFrom MASS mvrnorm
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

sim1.9 <- function(n, cen = 0) {
    V <- .75^as.matrix(dist(1:10, upper = TRUE))
    z <- mvrnorm(n = n, mu = rep(0, 10), Sigma = V)
    Time <- sqrt(rexp(n) * exp(- z %*% rep(c(-.5, .5), 5)))
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 3.77)
    if (cen == .50) cens <- runif(n, 0, 1.78)
    Y <- pmin(Time, cens)
    d <- 1 * (Time <= cens) 
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), d[x]))))
    dat <- cbind(dat, z[dat$id,])
    names(dat)[4:13] <- paste("z", 1:10, sep = "")
    dat   
}

sim2.1 <- function(n, cen = 0) {
    e <- rbinom(n, 1, .5)
    u <- rexp(n, 5)
    z2 <- runif(n)
    Time <- rep(NA, n)
    for (i in 1:n) {
        sol <- rexp(1)
        if (e[i] == 1)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u[i]) * exp(z2[i]) * x^2 -
                (x >= u[i]) * exp(z2[i]) * (x^2 * exp(1) + u[i]^2 * (1 - exp(1))),
                interval = c(0, 50))$root
        if (e[i] == 0)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u[i]) * exp(z2[i] + 1) * x^2 -
                (x > u[i]) * exp(z2[i]) * (exp(1) * x^2 + x^2 - u[i]^2),
                interval = c(0, 50))$root
    }
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 1.73)
    if (cen == .50) cens <- runif(n, 0, 0.83)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   z2 = z2[x], e = e[x], u = u[x])))
    dat$z1 <- with(dat, e * (Y < u) + (1 - e) * (Y >= u))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "e", "u")])
}

sim2.2 <- function(n, cen = 0) {
    e <- rbinom(n, 1, .5)
    u <- matrix(rexp(3 * n, 5), n)
    u <- t(apply(u, 1, sort))
    u1 <- u[,1]
    u2 <- u[,2]
    u3 <- u[,3]
    z2 <- runif(n)
    Time <- rep(NA, n)
    for (i in 1:n) {
        sol <- rexp(1)
        if (e[i] == 1)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u1[i]) * exp(z2[i]) * x^2 - 
                (x >= u1[i]) * (x < u2[i]) * exp(z2[i]) * (u1[i]^2 + exp(1) * (x^2 - u1[i]^2)) -
                (x >= u2[i]) * (x < u3[i]) * exp(z2[i]) * (u1[i]^2 + exp(1) * (u2[i]^2 - u1[i]^2) + x^2 - u2[i]^2) -
                (x >= u3[i]) * exp(z2[i]) * (u1[i]^2 + u3[i]^2 - u2[i]^2 + exp(1) * (u2[i]^2 - u1[i]^2 + x^2 - u3[i]^2)),
                interval = c(0, 50))$root
        if (e[i] == 0)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u1[i]) * exp(z2[i]) * exp(1) * x^2 -
                (x >= u1[i]) * (x < u2[i]) * exp(z2[i]) * (exp(1) * u1[i]^2 + x^2 - u1[i]^2) - 
                (x >= u2[i]) * (x < u3[i]) * exp(z2[i]) * (exp(1) * (u1[i]^2 + x^2 - u2[i]^2) + u2[i]^2 - u1[i]^2) -
                (x >= u3[i]) * exp(z2[i]) * (exp(1) * (u1[i]^2 + u3[i]^2 - u2[i]^2) + u2[i]^2 - u1[i]^2 - u3[i]^2 + x^2),
                interval = c(0, 50))$root
    }
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 2.11) 
    if (cen == .50) cens <- runif(n, 0, 1.02) 
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   z2 = z2[x], e = e[x], u1 = u1[x], u2 = u2[x], u3 = u3[x])))
    dat$z1 <- with(dat, e * (u1 <= Y) * (Y < u2) + e * (u3 <= Y) + (1 - e) * (Y < u1) + (1 - e) * (u2 <= Y) * (Y < u3))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "e", "u1", "u2", "u3")])
}

sim2.3 <- function(n, cen = 0) {
    k <- runif(n, 1, 2)
    b <- runif(n, 1, 2)
    z2 <- runif(n)
    Time <- rep(NA, n)
    for (i in 1:n) {
        sol <- rexp(1)
        Time[i] <- uniroot(f = function(x)
            sol - 2 * exp(z2[i] + b[i]) * (x * exp(k[i] * x) / k[i] - (exp(k[i] * x) - 1) / k[i]^2),
            interval = c(0, 50))$root
    }
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 1.11)
    if (cen == .50) cens <- runif(n, 0, 0.56)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   z2 = z2[x], k = k[x], b = b[x])))
    dat$z1 <- with(dat, k * Y + b)
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "k", "b")])
}

sim3.1 <- function(n, cen = 0) {
    e <- rbinom(n, 1, .5)
    u <- rexp(n, 5)
    z2 <- runif(n)
    Time <- rep(NA, n)    
    for (i in 1:n) {
        sol <- rexp(1)
        if (e[i] == 1)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u[i]) * .1 * exp(z2[i]) * x -
                (x >= u[i]) * .1 * exp(z2[i]) * (u[i] + exp(1) * (x - u[i])),
                interval = c(0, 100))$root
        if (e[i] == 0)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u[i]) * .1 * exp(z2[i] + 1) * x -
                (x > u[i]) * .1 * exp(z2[i]) * (exp(1) * u[i] + x + u[i]),
                interval = c(0, 100))$root
    }
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 14.99)
    if (cen == .50) cens <- runif(n, 0, 5.73)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   z2 = z2[x], e = e[x], u = u[x])))
    dat$z1 <- with(dat, e * (Y < u) + (1 - e) * (Y >= u))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "e", "u")])
}

sim3.2 <- function(n, cen = 0) {
    e <- rbinom(n, 1, .5)
    u <- matrix(rexp(3 * n, 5), n)
    u <- t(apply(u, 1, sort))
    u1 <- u[,1]
    u2 <- u[,2]
    u3 <- u[,3]
    z2 <- runif(n)
    Time <- rep(NA, n)
    for (i in 1:n) {
        sol <- rexp(1)
        if (e[i] == 1)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u1[i]) * .1 * exp(z2[i]) * x - 
                (x >= u1[i]) * (x < u2[i]) * .1 * exp(z2[i]) * (u1[i] + exp(1) * (x - u1[i])) -
                (x >= u2[i]) * (x < u3[i]) * .1 * exp(z2[i]) * (u1[i] + exp(1) * (u2[i] - u1[i]) + x - u2[i]) -
                (x >= u3[i]) * .1 * exp(z2[i]) * (u1[i] + u3[i] - u2[i] + exp(1) * (u2[i] - u1[i] + x - u3[i])),
                interval = c(0, 100))$root
        if (e[i] == 0)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u1[i]) * .1 * exp(z2[i]) * exp(1) * x -
                (x >= u1[i]) * (x < u2[i]) * .1 * exp(z2[i]) * (exp(1) * u1[i] + x - u1[i]) - 
                (x >= u2[i]) * (x < u3[i]) * .1 * exp(z2[i]) * (exp(1) * (u1[i] + x - u2[i]) + u2[i] - u1[i]) -
                (x >= u3[i]) * .1 * exp(z2[i]) * (exp(1) * (u1[i] + u3[i] - u2[i]) + u2[i] - u1[i] - u3[i] + x),
                interval = c(0, 100))$root
    }
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 15.54)
    if (cen == .50) cens <- runif(n, 0, 5.78)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   z2 = z2[x], e = e[x], u1 = u1[x], u2 = u2[x], u3 = u3[x])))
    dat$z1 <- with(dat, e * (u1 <= Y) * (Y < u2) + e * (u3 <= Y) + (1 - e) * (Y < u1) + (1 - e) * (u2 <= Y) * (Y < u3))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "e", "u1", "u2", "u3")])
}

sim3.3 <- function(n, cen = 0) {
    k <- runif(n, 1, 2)
    b <- runif(n, 1, 2)
    z2 <- runif(n)
    Time <- log(10 * rexp(n) * exp(-z2 - b) * k + 1) / k
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 2.48)
    if (cen == .50) cens <- runif(n, 0, 1.23)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   z1 = sort(Y[Y <= Y[x]]) * k[x] + b[x],
                   z2 = z2[x], k = k[x], b = b[x])))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "k", "b")])
}

#' #############################################################################################################
#' #############################################################################################################
#' Background functions for true survival and cumulative hazard curves given case and datasets
#' @keywords internal
#' @noRd
trueHaz1.1 <- function(dat) with(dat, Y^2 * exp(2 * z1 + 2 * z2))
trueSurv1.1 <- function(dat) exp(-with(dat, Y^2 * exp(2 * z1 + 2 * z2)))

trueHaz1.2 <- function(dat) with(dat, Y^2 * exp(2 * sin(2 * pi * dat$z1) + 2 * abs(dat$z2 - .5)))
trueSurv1.2 <- function(dat) with(dat, exp(-Y^2 * exp(2 * sin(2 * pi * dat$z1) + 2 * abs(dat$z2 - .5))))

trueHaz1.3 <- function(dat) with(dat, -log(pnorm(Y / exp(-2 + 2 * dat$z1 + 2 * dat$z2), sd = .5)))
trueSurv1.3 <- function(dat) with(dat, pnorm(Y / exp(-2 + 2 * dat$z1 + 2 * dat$z2), sd = .5))

trueHaz1.4 <- function(dat) with(dat, -log(pgamma(exp(z2 * Y / z1) / 4 / z2^2, 1 / 4 / z2^2, 1)))
trueSurv1.4 <- function(dat) with(dat, pgamma(exp(z2 * Y / z1) / 4 / z2^2, 1 / 4 / z2^2, 1))

trueHaz1.5 <- function(dat) trueHaz1.1(dat)
trueSurv1.5 <- function(dat) trueSurv1.1(dat)

trueHaz1.7 <- function(dat) trueHaz1.1(dat)
trueSurv1.7 <- function(dat) trueSurv1.1(dat)

trueHaz1.6 <- function(dat) {
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    dat$Y^2 * exp(2 * rowSums(z[,1:2]))
}
trueSurv1.6 <- function(dat) {
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    exp(-dat$Y^2 * exp(2 * rowSums(z[,1:2])))
}

trueHaz1.8 <- function(dat) {
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    dat$Y^2 * exp(.5 * rowSums(z))
}
trueSurv1.8 <- function(dat) {
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    exp(-dat$Y^2 * exp(.5 * rowSums(z)))
}

trueHaz1.9 <- function(dat) {
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    dat$Y^2 * exp(as.matrix(z) %*% rep(c(-.5, .5), 5))
}
trueSurv1.9 <- function(dat) {
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    exp(-dat$Y^2 * exp(as.matrix(z) %*% rep(c(-.5, .5), 5)))
}

trueSurv2.1 <- function(dat) {
    e <- dat$e[1]
    u <- dat$u[1]
    z2 <- dat$z2[1]
    Y <- dat$Y
    oneSurv <- function(Y) {
        if (e == 1) {
            if (Y < u) return(exp(-exp(z2) * Y^2))
            else return(exp(-(exp(z2) * u^2 + exp(z2 + 1) * (Y^2 - u^2))))
        }
        if (e == 0) {
            if (Y < u) return(exp(-exp(z2 + 1) * Y^2))
            else return(exp(-(exp(z2 + 1) * u^2 + exp(z2) * (Y^2 - u^2))))
        }    
    }
    return(sapply(Y, oneSurv))
}

trueHaz2.1 <- function(dat) {
    e <- dat$e[1]
    u <- dat$u[1]
    Y <- dat$Y
    z2 <- dat$z2[1]
    oneHaz <- function(Y) {
        if (e == 1) {
            if (Y < u) return(exp(z2) * Y^2)
            else return(exp(z2) * u^2 + exp(z2 + 1) * (Y^2 - u^2))
        }
        if (e == 0) {
            if (Y < u) return(exp(z2 + 1) * Y^2)
            else return(exp(z2 + 1) * u^2 + exp(z2) * (Y^2 - u^2))
        }
    }
    return(sapply(Y, oneHaz))
}
        
trueHaz2.2 <- function(dat) {
    e <- dat$e[1]
    u1 <- dat$u1[1]
    u2 <- dat$u2[1]
    u3 <- dat$u3[1]
    Y <- dat$Y
    z2 <- dat$z2[1]
    oneHaz <- function(Y) {
        if (e == 1) {
            if (Y < u1) return(exp(z2 + 1) * Y^2)
            if (Y >= u1 & Y < u2) return(exp(z2) * (exp(1) * u1^2 + Y^2 - u1^2))
            if (Y >= u2 & Y < u3) return(exp(z2) * (exp(1) * u1^2 + (u2^2 - u1^2) + exp(1) * (Y^2 - u2^2)))
            if (Y > u3) return(exp(z2) * (exp(1) * u1^2 + u2^2 - u1^2 + exp(1) * (u3^2 - u2^2) + (Y^2 - u3^2)))
        }
        if (e == 0) {
            if (Y < u1) return(exp(z2) * Y^2 )
            if (Y >= u1 & Y < u2) return(exp(z2) * (u1^2 + exp(1) * (Y^2 - u1^2)))
            if (Y >= u2 & Y < u3) return(exp(z2) * (u1^2 + exp(1) * (u2^2 - u1^2) + (Y^2 - u2^2)))
            if (Y > u3) return(exp(z2) * (u1^2 + exp(1) * (u2^2 - u1^2) + (u3^2 - u2^2) + exp(1) * (Y^2 - u3^2)))
        }
    }
    return(sapply(Y, oneHaz))
}

trueSurv2.2 <- function(dat) {
    e <- dat$e[1]
    u1 <- dat$u1[1]
    u2 <- dat$u2[1]
    u3 <- dat$u3[1]
    Y <- dat$Y
    z2 <- dat$z2[1]
    oneSurv <- function(Y) {
        if (e == 1) {
            if (Y < u1) return(exp(-exp(z2 + 1) * Y^2))
            if (Y >= u1 & Y < u2) return(exp(-(exp(z2) * (exp(1) * u1^2 + Y^2 - u1^2))))
            if (Y >= u2 & Y < u3) return(exp(-(exp(z2) * (exp(1) * u1^2 + (u2^2 - u1^2) + exp(1) * (Y^2 - u2^2)))))
            if (Y > u3) return(exp(-(exp(z2) * (exp(1) * u1^2 + u2^2 - u1^2 + exp(1) * (u3^2 - u2^2) + (Y^2 - u3^2)))))
        }
        if (e == 0) {
            if (Y < u1) return(exp(-(exp(z2) * Y^2 )))
            if (Y >= u1 & Y < u2) return(exp(-(exp(z2) * (u1^2 + exp(1) * (Y^2 - u1^2)))))
            if (Y >= u2 & Y < u3) return(exp(-(exp(z2) * (u1^2 + exp(1) * (u2^2 - u1^2) + (Y^2 - u2^2)))))
            if (Y > u3) return(exp(-(exp(z2) * (u1^2 + exp(1) * (u2^2 - u1^2) + (u3^2 - u2^2) + exp(1) * (Y^2 - u3^2)))))
        }
    }
    return(sapply(Y, oneSurv))
}
    
trueHaz3.1 <- function(dat) {
    e <- dat$e[1]
    u <- dat$u[1]
    Y <- dat$Y
    z2 <- dat$z2[1]
    oneHaz <- function(Y) {
        if (e == 1) {
            if (Y < u) return(.1 * exp(z2) * Y)
            else return(.1 * exp(z2) * u + .1 * exp(z2 + 1) * (Y - u))
        }
        if (e == 0) {
            if (Y < u) return(exp(z2 + 1) * Y^2)
            else return(.1 * exp(z2 + 1) * u + .1 * exp(z2) * (Y - u))
        }
    }
    return(sapply(Y, oneHaz))
}
    
trueSurv3.1 <- function(dat) {
    e <- dat$e[1]
    u <- dat$u[1]
    Y <- dat$Y
    z2 <- dat$z2[1]
    oneSurv <- function(Y) {
        if (e == 1) {
            if (Y < u) return(exp(-.1 * exp(z2) * Y))
            else return(exp(-(.1 * exp(z2) * u + .1 * exp(z2 + 1) * (Y - u))))
        }
        if (e == 0) {
            if (Y < u) return(exp(-exp(z2 + 1) * Y^2))
            else return(exp(-(.1 * exp(z2 + 1) * u + .1 * exp(z2) * (Y - u))))
        }
    }
    return(sapply(Y, oneSurv))
}

trueHaz3.2 <- function(dat) {
    e <- dat$e[1]
    u1 <- dat$u1[1]
    u2 <- dat$u2[1]
    u3 <- dat$u3[1]
    Y <- dat$Y
    z2 <- dat$z2[1]
    oneHaz <- function(Y) {
        if (e == 1) {
            if (Y < u1) return(.1 * exp(z2 + 1) * Y)
            if (Y >= u1 & Y < u2) return(.1 * exp(z2) * (exp(1) * u1 + Y - u1))
            if (Y >= u2 & Y < u3) return(.1 * exp(z2) * (exp(1) * u1 + (u2 - u1) + exp(1) * (Y - u2)))
            if (Y > u3) return(.1 * exp(z2) * (exp(1) * u1 + (u2 - u1) + exp(1) * (u3 - u2) + (Y - u3)))
        }
        if (e == 0) {
            if (Y < u1) return(.1 * exp(z2) * Y)
            if (Y >= u1 & Y < u2) return(.1 * exp(z2) * (u1 + exp(1) * (Y - u1)))
            if (Y >= u2 & Y < u3) return(.1 * exp(z2) * (u1 + exp(1) * (u2 - u1) + (Y - u2)))
            if (Y > u3) return(.1 * exp(z2) * (u1 + exp(1) * (u2 - u1) + (u3 - u2) + exp(1) * (Y - u3)))
        }
    }
    return(sapply(Y, oneHaz))
}

trueSurv3.2 <- function(dat) {
    e <- dat$e[1]
    u1 <- dat$u1[1]
    u2 <- dat$u2[1]
    u3 <- dat$u3[1]
    Y <- dat$Y
    z2 <- dat$z2[1]
    oneSurv <- function(Y) {
        if (e == 1) {
            if (Y < u1) return(exp(-.1 * exp(z2 + 1) * Y))
            if (Y >= u1 & Y < u2) return(exp(-.1 * exp(z2) * (exp(1) * u1 + Y - u1)))
            if (Y >= u2 & Y < u3) return(exp(-.1 * exp(z2) * (exp(1) * u1 + (u2 - u1) + exp(1) * (Y - u2))))
            if (Y > u3) return(exp(-.1 * exp(z2) * (exp(1) * u1 + (u2 - u1) + exp(1) * (u3 - u2) + (Y - u3))))
        }
        if (e == 0) {
            if (Y < u1) return(exp(-.1 * exp(z2) * Y))
            if (Y >= u1 & Y < u2) return(exp(-.1 * exp(z2) * (u1 + exp(1) * (Y - u1))))
            if (Y >= u2 & Y < u3) return(exp(-.1 * exp(z2) * (u1 + exp(1) * (u2 - u1) + (Y - u2))))
            if (Y > u3) return(exp(-.1 * exp(z2) * (u1 + exp(1) * (u2 - u1) + (u3 - u2) + exp(1) * (Y - u3))))
        }
    }
    return(sapply(Y, oneSurv))
}

trueHaz2.3 <- function(dat) {
    Y <- dat$Y
    k <- dat$k[1]
    b <- dat$b[1]
    z2 <- dat$z2[1]
    return(2 * exp(b + z2) * (Y * exp(k * Y) / k - (exp(k * Y) - 1) / k^2))
}

trueSurv2.3 <- function(dat) {
    Y <- dat$Y
    k <- dat$k[1]
    b <- dat$b[1]
    z2 <- dat$z2[1]
    return(exp(-(2 * exp(b + z2) * (Y * exp(k * Y) / k - (exp(k * Y) - 1) / k^2))))
}

trueHaz3.3 <- function(dat) {
    Y <- dat$Y
    k <- dat$k[1]
    b <- dat$b[1]
    z2 <- dat$z2[1]
    return(.1 * exp(z2 + b) * (exp(k * Y) - 1) / k)
}

trueSurv3.3 <- function(dat) {
    Y <- dat$Y
    k <- dat$k[1]
    b <- dat$b[1]
    z2 <- dat$z2[1]
    return(exp(-.1 * exp(z2 + b) * (exp(k * Y) - 1) / k))
}

#' Background functions for generating testing sets given training data
#'
#' These functions give ONE draw of test set
#' 
#' @keywords internal
#' @noRd
simuTest1.1 <- function(dat) {
    Y <- sort(unique(dat$Y))
    data.frame(Y = Y, z1 = runif(1), z2 = runif(1))
}

simuTest1.2 <- function(dat) simuTest1.1(dat)
simuTest1.3 <- function(dat) simuTest1.1(dat)
simuTest1.4 <- function(dat) simuTest1.1(dat)

simuTest1.6 <- function(dat) {
    Y <- sort(unique(dat$Y))
    data.frame(Y = Y, z1 = runif(1), z2 = runif(1), z3 = runif(1), z4 = runif(1), z5 = runif(1),
               z6 = runif(1), z7 = runif(1), z8 = runif(1), z9 = runif(1), z10 = runif(1))
}


simuTest1.7 <- function(dat) simuTest1.6(dat)
simuTest1.8 <- function(dat) {
    Y <- sort(unique(dat$Y))
    V <- .75^as.matrix(dist(1:10, upper = TRUE))
    z <- mvrnorm(n = 1, mu = rep(0, 10), Sigma = V)
    as.data.frame(cbind(Y = Y, z1 = z[1], z2 = z[2], z3 = z[3], z4 = z[4], z5 = z[5],
                        z6 = z[6], z7 = z[7], z8 = z[8], z9 = z[9], z10 = z[10]))
}
simuTest1.9 <- function(dat) simuTest1.8(dat)

simuTest1.5 <- function(dat) {
    Y <- sort(unique(dat$Y))
    data.frame(Y = Y, z1 = runif(1), z2 = runif(1), z3 = runif(1), z4 = runif(1), z5 = runif(1))
}

simuTest2.1 <- function(dat) {
    Y <- sort(unique(dat$Y))
    e <- rbinom(1, 1, .5)
    u <- rexp(1, 5)
    data.frame(Y = Y, z1 = e * (Y < u) + (1 - e) * (Y >= u), z2 = runif(1), e = e, u = u)
}

simuTest2.2 <- function(dat) {
    Y <- sort(unique(dat$Y))
    e <- rbinom(1, 1, .5)
    u <- sort(rexp(3 * 1, 5))
    u1 <- u[1]
    u2 <- u[2]
    u3 <- u[3]
    z1 <- e * ((u1 <= Y) * (Y < u2) + (u3 <= Y)) + (1 - e) * ((Y < u1) + (u2 <= Y) * (Y < u3))
    data.frame(Y = Y, z1 = z1, z2 = runif(1), e = e, u1 = u1, u2 = u2, u3 = u3)
}

simuTest2.3 <- function(dat) {
    Y <- sort(unique(dat$Y))
    k <- runif(1, 1, 2)
    b <- runif(1, 1, 2)
    data.frame(Y = Y, z1 = k * Y + b, z2 = runif(1), k = k, b = b)
}

simuTest3.1 <- function(dat) simuTest2.1(dat)
simuTest3.2 <- function(dat) simuTest2.2(dat)
simuTest3.3 <- function(dat) simuTest2.3(dat)
