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
#' @importFrom dplyr "%>%" arrange
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
    if (!(scenario %in% paste(rep(1:3, each = 6), rep(1:6, 3), sep = ".")))
        stop("See ?simu for scenario definition.")
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
#' This function is used to generate the true cumulative hazard function used in the simulation.
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
    if (!any(names(dat) == "id")) {
        dat <- dat %>% arrange(Y)
        dat$id <- unlist(lapply(split(dat$Y, dat$Y), function(x) 1:length(x)))
    }
    ## if (substr(scenario, 1, 1) == "1") dat <- do.call(rbind, lapply(split(dat, dat$id), function(x) x[which.max(x$Y),]))
    eval(parse(text = paste("trueHaz", scenario, "(dat)", sep = "")))    
}

#' Function to generate the true survival used in the simulation.
#'
#' This function is used to generate the true survival function used in the simulation.
#'
#' @rdname simu
#' @export
#' 
trueSurv <- function(dat) {
    if (attr(dat, "prepBy") != "rocSimu") stop("Inputed data must be prepared by \"simu\".")
    scenario <- attr(dat, "scenario")
    if (!any(names(dat) == "id")) {
        dat <- dat %>% arrange(Y)
        dat$id <- unlist(lapply(split(dat$Y, dat$Y), function(x) 1:length(x)))
    }
    ## if (substr(scenario, 1, 1) == "1") dat <- do.call(rbind, lapply(split(dat, dat$id), function(x) x[which.max(x$Y),]))
    eval(parse(text = paste("trueSurv", scenario, "(dat)", sep = "")))    
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

## sim3.1 <- function(n, cen = 0) {
##     e <- rbinom(n, 1, .5)
##     u0 <- rexp(n, 5)
##     z2 <- runif(n)
##     Time <- rep(NA, n)
##     u <- rexp(n)
##     T1 <- (u / .1 / exp(z2) < u0) * (u / .1 / exp(z2)) + (u / .1 / exp(z2) >= u0) * ((u / .1 / exp(z2) - u0) / exp(1) + u)
##     T2 <- (u / .1 / exp(z2) / exp(1) < u0) * (u / .1 / exp(z2) / exp(1)) + (u / .1 / exp(z2) / exp(1) >= u0) * ((u / .1 / exp(z2) - exp(1) * u) + u)
##     Time <- T1 * e + T2 * (1 - e)
##     if (cen == 0) cens <- rep(Inf, n)
##     if (cen == .25) cens <- runif(n, 0, 14.99)
##     if (cen == .50) cens <- runif(n, 0, 5.76)
##     Y <- pmin(Time, cens)
##     dat <- do.call(rbind, lapply(1:n, function(x)
##         data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
##                    z1 = e[x] * (sort(Y[Y <= Y[x]]) < u[x]) + (1 - e[x]) * (sort(Y[Y <= Y[x]]) >= u[x]),
##                    z2 = z2[x], e = e[x], u = u[x])))
##     return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "e", "u")])
## }

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

#' Background functions for true survival and cumulative hazard curves given case and datasets
#' @keywords internal
#' @noRd
trueHaz1.1 <- function(dat) {
    cumHaz <- with(dat, Y^2 * exp(2 * z1 + 2 * z2))
    approxfun(x = dat$Y, y = cumHaz, method = "constant", yleft = 0, yright = max(cumHaz))
}

trueSurv1.1 <- function(dat) {
    Surv <- exp(-with(dat, Y^2 * exp(2 * z1 + 2 * z2)))
    approxfun(x = dat$Y, y = Surv, method = "constant", yleft = 1, yright = min(Surv))
}

trueHaz1.2 <- function(dat) {
    cumHaz <- with(dat, Y^2 * exp(2 * sin(2 * pi * dat$z1) + 2 * abs(dat$z2 - .5)))
    approxfun(x = dat$Y, y = cumHaz, method = "constant", yleft = 0, yright = max(cumHaz)) 
}

trueSurv1.2 <- function(dat) {
    Surv <- with(dat, exp(-Y^2 * exp(2 * sin(2 * pi * dat$z1) + 2 * abs(dat$z2 - .5))))
    approxfun(x = dat$Y, y = Surv, method = "constant", yleft = 1, yright = min(Surv))
}

trueHaz1.3 <- function(dat) {
    cumHaz <- with(dat, -log(pnorm(Y / exp(-2 + 2 * dat$z1 + 2 * dat$z2), sd = .5)))
    approxfun(x = dat$Y, y = cumHaz, method = "constant", yleft = 0, yright = max(cumHaz))
}

trueSurv1.3 <- function(dat) {
    Surv <- with(dat, pnorm(Y / exp(-2 + 2 * dat$z1 + 2 * dat$z2), sd = .5))
    approxfun(x = dat$Y, y = Surv, method = "constant", yleft = 1, yright = min(Surv))
}

trueHaz1.4 <- function(dat) {
    cumHaz <- with(dat, -log(pgamma(exp(z2 * Y / z1) / 4 / z2^2, 1 / 4 / z2^2, 1)))
    approxfun(x = dat$Y, y = cumHaz, method = "constant", yleft = 0, yright = max(cumHaz))
}

trueSurv1.4 <- function(dat) {
    Surv <- with(dat, pgamma(exp(z2 * Y / z1) / 4 / z2^2, 1 / 4 / z2^2, 1))
    approxfun(x = dat$Y, y = Surv, method = "constant", yleft = 1, yright = min(Surv))
}

trueHaz1.5 <- function(dat) trueHaz1.1(dat)
trueSurv1.5 <- function(dat) trueSurv1.1(dat)

#' Background functions for generating testing sets given training data
#' @keywords internal
#' @noRd
simuTest1.1 <- function(dat) {
    Y <- sort(unique(dat$Y))
    data.frame(Y = Y, z1 = runif(1), z2 = runif(1))
}

simuTest1.2 <- function(dat) simuTest1.1(dat)
simuTest1.3 <- function(dat) simuTest1.1(dat)
simuTest1.4 <- function(dat) simuTest1.1(dat)
simuTest1.5 <- function(dat) simuTest1.1(dat)

simuTest2.1 <- function(dat) {
    Y <- sort(unique(dat$Y))
    n <- length(Y)
    e <- rbinom(n, 1, .5)
    u <- rexp(n, 5)
    data.frame(Y = Y, z1 = e * (Y < u) + (1 - e) * (Y >= u), z2 = runif(1))
}

simuTest2.2 <- function(dat) {
    Y <- sort(unique(dat$Y))
    e <- rbinom(n, 1, .5)
    u <- matrix(rexp(3 * n, 5), n)
    u1 <- u[,1]
    u2 <- u[,2]
    u3 <- u[,3]
    z1 <- e * ((u1 <= Y) * (Y < u2) + (u3 <= Y)) + (1 - e) * ((Y < u1) + (u2 <= Y) * (Y < u3))
    data.frame(Y = Y, z1 = z1, z2 = runif(1))
}

simuTest2.3 <- function(dat) {
    Y <- sort(unique(dat$Y))
    n <- length(Y)
    k <- runif(n, 1, 2)
    b <- runif(n, 1, 2)
    data.frame(Y = Y, z1 = k * Y + b, z2 = runif(1))
}

simuTest3.1 <- function(dat) simuTest2.1(dat)
simuTest3.2 <- function(dat) simuTest2.2(dat)
simuTest3.3 <- function(dat) simuTest2.3(dat)
