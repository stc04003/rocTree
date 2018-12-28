globalVariables(c("n", "cen", "Y", "id", "z2", "e", "u", "y", "z2", "u1", "u2", "u3")) ## global variables for simu

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
#' The following scenarios form the paper are included.
#' \describe{
#' \item{Scenario 1}{assumes time independent covariates with \eqn{\Lambda(0) = 2t}:}
#' \describe{
#' \item{1.1}{Proportional hazards model:
#' \deqn{\lambda(t, Z) = \lambda_0(t) e^{(-0.5 Z_1 + 0.5 Z_2 - 0.5 Z_3 ... + 0.5 Z_{10})}.}}
#' \item{1.2}{Proportional hazards model with noise variable:
#' \deqn{\lambda(t, Z) = \lambda_0(t) e^{(2Z_1 + 2Z_2 + 0Z_3 + ... + 0Z_{10})}.}}
#' \item{1.3}{Proportional hazards model with nonlinear covariate effects:
#' \deqn{\lambda(t, Z) = \lambda_0(t) e^{[2\sin(2\pi Z_1) + 2|Z_2 - 0.5|]}.}}
#' \item{1.4}{Accelerated failure time model:
#' \deqn{\log(T) = -2 + 2Z_1 + 2Z_2 + \epsilon,} where \eqn{\epsilon} follows \eqn{N(0, 0.5^2).}}
#' \item{1.5}{Generalized gamma family:
#' \deqn{T = e^{\sigma\omega},} with \eqn{\omega = \log(Q^2 g) / Q}, \eqn{g} follows Gamma(\eqn{Q^{-2}, 1}),
#' \eqn{\sigma = 2Z_1, Q = 2Z_2.}}
#' }
#' \item{Scenario 2}{assumes a time dependent covariate, \eqn{Z_1(t)}:}
#' \describe{
#' \item{2.1}{Dichotomous time dependent covariate with at most one change in value:
#'  \deqn{\lambda(t, Z(t)) = \lambda_0(t)e^{2Z_1(t) + 2Z_2},}
#' where \eqn{Z_1(t)} is the time-dependent covariate: \eqn{Z_1(t) = \theta I(t \ge U_0) + (1 - \theta) I(t < U_0)},
#' and \eqn{\theta} is a Bernoulli variable with equal probability.}
#' \item{2.2}{Dichotomous time dependent covariate with multiple changes:
#'  \deqn{\lambda(t, Z(t)) = e^{2Z_1(t) + 2Z_2},}
#' where \eqn{Z_1(t) = \theta[I(U_1\le t < U_2) + I(U_3 \le t)] + (1 - \theta)[I(t < U_1) + I(U_2\le t < U_3)]},
#' \eqn{\theta} is a Bernoulli variable with equal probability, and \eqn{U_1\le U_2\le U_3}
#' are the first three terms of a stationary Poisson process with rate 10.}
#' \item{2.3}{Continuous time dependent covariate:
#' \deqn{\lambda(t, Z(t)) = 0.1 e^{Z_1(t) + Z_2},} where \eqn{Z_1(t) = kt + b},
#' \eqn{k} and \eqn{b} follow independent uniform distributions over \eqn{(1, 2)}.
#' }
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
#' @importFrom MASS mvrnorm
#' @importFrom stats dist quantile rlnorm
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
#' This function is used to generate the true cumulative hazard function for
#' ONE ID (the 1st ID) used in the simulation.
#' This function is usually called for the testing set.
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

sim1.3 <- function(n, cen = 0) {
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

sim1.4 <- function(n, cen = 0) {
    z1 <- runif(n)
    z2 <- runif(n)
    Time <- rlnorm(n, meanlog = -2 + 2 * z1 + 2 * z2, sdlog = .5)
    ## exp(-2 + 2 * z1 + 2 * z2 + rnorm(n, sd = .5))
    cens <- runif(n, 0, cen)
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 6.00)
    if (cen == .50) cens <- runif(n, 0, 2.40)
    Y <- pmin(Time, cens)
    do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])), z1 = z1[x], z2 = z2[x])))
}

#' @importFrom flexsurv pgengamma rgengamma 
sim1.5 <- function(n, cen = 0) {
    z1 <- runif(n)
    z2 <- runif(n)
    Time <- rgengamma(n, mu = 0, sigma = 2 * z1, Q = 2 * z2)
    cens <- runif(n, 0, cen)
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 4.12)
    if (cen == .50) cens <- runif(n, 0, 1.63)
    Y <- pmin(Time, cens)
    do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])), z1 = z1[x], z2 = z2[x])))
}

sim1.2 <- function(n, cen = 0) {
    z <- matrix(runif(n * 10), n, 10)
    if (n == 1) Time <- sqrt(rexp(n) * exp(-sum(z[,1:2]) * 2))
    else Time <- sqrt(rexp(n) * exp(-rowSums(z[,1:2]) * 2))
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

sim1.1 <- function(n, cen = 0) {
    V <- .75^as.matrix(dist(1:10, upper = TRUE))
    z <- matrix(mvrnorm(n = n, mu = rep(0, 10), Sigma = V), ncol = 10)
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

sim2.1 <- function(n, cen) {
    z2 <- runif(n)
    p1 <- rbinom(n, 1, .5)
    t0 <- rexp(n, 5)
    nu <- 2
    U <- runif(n)
    b1 <- 2
    b2 <- 2
    a <- b2 * z2
    Tp0 <- ifelse(-log(U) < exp(a) * t0^nu, 
    (-log(U) * exp(-a))^(1/nu),
    ((-log(U) - exp(a) * t0^nu + exp(a + b1) * t0^nu) / exp(a + b1))^(1/nu))
    Tp1 <- ifelse(-log(U) < exp(a + b1) * t0^nu,
    (-log(U) / exp(a + b1))^(1/nu),
    ((-log(U) - exp(a + b1) * t0^nu + exp(a) * t0^nu) / exp(a))^(1/nu))
    Time <- Tp0 * p1 + Tp1 * (1 - p1)
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 1.8)
    if (cen == .50) cens <- runif(n, 0, 0.8)
    Y <- pmin(Time, cens)
    d <- 1 * (Time <= cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), d[x]),
                   z2 = z2[x], e = p1[x], u = t0[x])))
    dat$z1 <- with(dat, e * (Y < u) + (1 - e) * (Y >= u))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "e", "u")])        
}

sim2.2 <- function(n, cen) {
    U <- runif(n)
    z2 <- runif(n)
    b1 <- 2
    b2 <- 2
    a <- b2 * z2
    p1 <- rbinom(n, 1, .5)
    t1 <- rexp(n, 10)
    t2 <- t1 + rexp(n, 10)
    t3 <- t2 + rexp(n, 10)
    Tp0 <- (-log(U) / exp(a)) * (-log(U) < exp(a) * t1) + 
        ((-log(U) - exp(a) * t1 + exp(a + b1) * t1) / exp(a + b1)) *
        (-log(U) >= exp(a) * t1 & -log(U) < exp(a) * (t1 + exp(b1) * (t2 - t1))) +
        ((-log(U) - exp(a) * t1 - exp(a + b1) * (t2 - t1) + exp(a) * t2) / exp(a)) *
        (-log(U) >= exp(a) * (t1 + exp(b1) * (t2-t1)) & -log(U) < exp(a) *
         (t1 + exp(b1) * (t2 - t1) + (t3 - t2))) +
        ((-log(U) - exp(a) * t1 - exp(a + b1) * (t2 - t1) -
          exp(a) * (t3 - t2) + exp(a + b1) * t3) / exp(a + b1)) *
        (-log(U) >= exp(a) * (t1 + exp(b1) * (t2 - t1) + (t3 - t2)))
    Tp1 <- (-log(U) / exp(a + b1)) * (-log(U) < exp(a + b1) * t1) + 
        ((-log(U) - exp(a + b1) * t1 + exp(a) * t1) / exp(a)) *
        (-log(U) >= exp(a + b1) * t1 & -log(U) < exp(a) * (exp(b1) * t1 + (t2 - t1)))+
        ((-log(U) - exp(a + b1) * t1 - exp(a) * (t2 - t1) + exp(a + b1) * t2) / exp(a + b1)) *
        (-log(U) >= exp(a) * (exp(b1) * t1 + (t2 - t1)) & -log(U) < exp(a) *
         (exp(b1) * t1 + (t2 - t1) + exp(b1) * (t3 - t2))) +
        ((-log(U) - exp(a + b1) * t1 - exp(a) * (t2 - t1) -
          exp(a + b1) * (t3 - t2) + exp(a) * t3) / exp(a)) *
        (-log(U) >= exp(a) * (exp(b1) * t1 + (t2 - t1) + exp(b1) * (t3 - t2)))
    Time <- Tp0 * p1 + Tp1 * (1 - p1)
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 0.7)
    if (cen == .50) cens <- runif(n, 0, 0.275)
    Y <- pmin(Time, cens)
    d <- 1 * (Time <= cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]), death = c(rep(0, sum(Y < Y[x])), d[x]),
                   z2 = z2[x], e = p1[x], u1 = t1[x], u2 = t2[x], u3 = t3[x])))
    dat$z1 <- with(dat, e * (u1 <= Y) * (Y < u2) + e * (u3 <= Y) + (1 - e) * (Y < u1) + (1 - e) * (u2 <= Y) * (Y < u3))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "e", "u1", "u2", "u3")])      
}

sim2.3 <- function(n, cen = 0) {
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
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    dat$Y^2 * exp(as.matrix(z) %*% rep(c(-.5, .5), 5))
}
trueSurv1.1 <- function(dat) {
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    exp(-dat$Y^2 * exp(as.matrix(z) %*% rep(c(-.5, .5), 5)))
}
    
trueHaz1.2 <- function(dat) {
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    dat$Y^2 * exp(2 * rowSums(z[,1:2]))
}
trueSurv1.2 <- function(dat) {
    z <- dat %>% select(paste("z", 1:10, sep = ""))
    exp(-dat$Y^2 * exp(2 * rowSums(z[,1:2])))
}

trueHaz1.3 <- function(dat) with(dat, Y^2 * exp(2 * sin(2 * pi * z1) + 2 * abs(z2 - .5)))
trueSurv1.3 <- function(dat) with(dat, exp(-Y^2 * exp(2 * sin(2 * pi * z1) + 2 * abs(z2 - .5))))

trueHaz1.4 <- function(dat) with(dat, -log(1 - plnorm(Y, meanlog = -2 + 2 * z1 + 2 * z2, sdlog = .5)))
trueSurv1.4 <- function(dat) with(dat, 1 - plnorm(Y, meanlog = -2 + 2 * z1 + 2 * z2, sdlog = .5))

trueHaz1.5 <- function(dat) with(dat, -log(1 - pgengamma(Y, mu = 0, sigma = 2 * z1, Q = 2 * z2)))
trueSurv1.5 <- function(dat) with(dat, 1 - pgengamma(Y, mu = 0, sigma = 2 * z1, Q = 2 * z2))

trueHaz2.1 <- function(dat) {
    attach(dat)
    b1 <- 2
    b2 <- 2
    a <- b2 * z2
    nu <- 2
    Surv <- exp(a) * Y^nu * (e == 1 & Y < u) +
        (exp(a) * (u^nu + exp(b1) * Y^nu - exp(b1) * u^nu)) * (e == 1 & Y >= u) +
        (exp(a + b1) * Y^nu) * (e == 0 & Y < u) +
        (exp(a) * (exp(b1) * u^nu + Y^nu - u^nu)) * (e == 0 & Y >= u)
    detach(dat)
    return(Surv)
}
trueSurv2.1 <- function(dat) {
    haz <- trueHaz2.1(dat)
    return(exp(-haz))
}

trueHaz2.2 <- function(dat) {
    attach(dat)
    b1 <- 2
    b2 <- 2
    a <- b2 * z2
    haz <- (exp(a) * Y) * (e == 1 & Y < u1) +
        (exp(a) * (u1 + exp(b1) * Y - exp(b1) * u1)) * (e == 1 & Y >= u1 & Y < u2) + 
        (exp(a) * (u1 + exp(b1) * u2 - exp(b1) * u1 + Y - u2)) * (e == 1 & Y >= u2 & Y < u3) + 
        (exp(a) * (u1 + exp(b1) * u2 - exp(b1) * u1 + u3 - u2 + exp(b1) * (Y - u3))) *
        (e == 1 & Y >= u3) + 
        (exp(a + b1) * Y) * (e == 0 & Y < u1) + (exp(a) * (exp(b1) * u1 + Y - u1)) *
        (e == 0 & Y >= u1 & Y < u2) +
        (exp(a) * (exp(b1) * u1 + u2 - u1 + exp(b1) * (Y - u2))) * (e == 0 & Y >= u2 & Y < u3) + 
        (exp(a) * (exp(b1) * u1 + u2 - u1 + exp(b1) * (u3 - u2) + (Y - u3))) * (e == 0 & Y >= u3)    
    detach(dat)
    haz
}

trueSurv2.2 <- function(dat) {
    haz <- trueHaz2.2(dat)
    return(exp(-haz))
}

trueSurv2.3 <- function(dat) {
    Y <- dat$Y
    k <- dat$k[1]
    b <- dat$b[1]
    z2 <- dat$z2[1]
    return(exp(-.1 * exp(z2 + b) * (exp(k * Y) - 1) / k))
}
trueHaz2.3 <- function(dat) {
    Y <- dat$Y
    k <- dat$k[1]
    b <- dat$b[1]
    z2 <- dat$z2[1]
    return(.1 * exp(z2 + b) * (exp(k * Y) - 1) / k)
}

#' Background functions for generating testing codes
#' @keywords internal
#' @noRd

simuTest1.1 <- function(dat) {
    Y <- sort(unique(dat$Y))
    V <- .75^as.matrix(dist(1:10, upper = TRUE))
    z <- mvrnorm(n = 1, mu = rep(0, 10), Sigma = V)
    as.data.frame(cbind(Y = Y, z1 = z[1], z2 = z[2], z3 = z[3], z4 = z[4], z5 = z[5],
                        z6 = z[6], z7 = z[7], z8 = z[8], z9 = z[9], z10 = z[10]))
}
simuTest1.2 <- function(dat) {
    Y <- sort(unique(dat$Y))
    data.frame(Y = Y, z1 = runif(1), z2 = runif(1), z3 = runif(1), z4 = runif(1), z5 = runif(1),
               z6 = runif(1), z7 = runif(1), z8 = runif(1), z9 = runif(1), z10 = runif(1))
}
simuTest1.3 <- function(dat) {
    Y <- sort(unique(dat$Y))
    data.frame(Y = Y, z1 = runif(1), z2 = runif(1))
}
simuTest1.4 <- function(dat) simuTest1.3(dat)
simuTest1.5 <- function(dat) simuTest1.3(dat)

simuTest2.1 <- function(dat) {
    Y <- sort(unique(dat$Y))
    e <- rbinom(1, 1, .5)
    u <- rexp(1, 5)
    data.frame(Y = Y, z1 = e * (Y < u) + (1 - e) * (Y >= u), z2 = runif(1), e = e, u = u)
}
simuTest2.2 <- function(dat) {
    Y <- sort(unique(dat$Y))
    e <- rbinom(1, 1, .5)
    u1 <- rexp(1, 10)
    u2 <- u1 + rexp(1, 10)
    u3 <- u2 + rexp(1, 10)
    z1 <- e * ((u1 <= Y) * (Y < u2) + (u3 <= Y)) + (1 - e) * ((Y < u1) + (u2 <= Y) * (Y < u3))
    data.frame(Y = Y, z1 = z1, z2 = runif(1), e = e, u1 = u1, u2 = u2, u3 = u3)
}
simuTest2.3 <- function(dat) {
    Y <- sort(unique(dat$Y))
    k <- runif(1, 1, 2)
    b <- runif(1, 1, 2)
    data.frame(Y = Y, z1 = k * Y + b, z2 = runif(1), k = k, b = b)
}


#' Continuious time varying covariate with non-Cox model
#'
#' \lambda(t, Z(t)) = 0.1 \cdot \left[\sin(Z_1(t) + Z_2) + 1\right],
#' where \eqn{Z_1(t) = kt + b},
#' \eqn{k} and \eqn{b} follow independent uniform distributions over \eqn{(1, 2)}.
sim2.4 <- function(n, cen = 0) {
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
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   z1 = sort(Y[Y <= Y[x]]) * k[x] + b[x],
                   z2 = z2[x], k = k[x], b = b[x])))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "k", "b")])
}
trueSurv2.4 <- function(dat) {
    haz <- trueHaz2.4(dat)
    return(exp(-haz))
}
trueHaz2.4 <- function(dat) {
    Y <- dat$Y
    k <- dat$k[1]
    b <- dat$b[1]
    z2 <- dat$z2[1]
    return(.1 * (Y - cos(k * Y + b + z2) / k + cos(b + z2) / k))
}
simuTest2.4 <- function(dat) {
    Y <- sort(unique(dat$Y))
    k <- runif(1, 1, 2)
    b <- runif(1, 1, 2)
    data.frame(Y = Y, z1 = k * Y + b, z2 = runif(1), k = k, b = b)
}


#' Continuious time varying covariate with non-Cox model and non-monotonic Z_1(t)
#'
#' \lambda(t, Z(t)) = 0.1 \cdot \left[\sin(Z_1(t) + Z_2) + 1\right],
#' where \eqn{Z_1(t) = kt * (2 * I(t > 5) - 1) + b},
#' \eqn{k} and \eqn{b} follow independent uniform distributions over \eqn{(1, 2)}.
sim2.5 <- function(n, cen = 0) {
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
        data.frame(id = x, Y = sort(Y[Y <= Y[x]]),
                   death = c(rep(0, sum(Y < Y[x])), 1 * (Time[x] <= cens[x])),
                   z1 = sort(Y[Y <= Y[x]]) * k[x] * (2 * (sort(Y[Y <= Y[x]]) > 5) - 1) + b[x],
                   z2 = z2[x], k = k[x], b = b[x])))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "k", "b")])
}

trueSurv2.5 <- function(dat) {
    haz <- trueHaz2.5(dat)
    return(exp(-haz))
}

trueHaz2.5 <- function(dat) {
    Y <- dat$Y
    k <- dat$k[1]
    b <- dat$b[1]
    z2 <- dat$z2[1]
    return(.1 * (Y - cos(k * Y * (2 * (Y > 5) - 1) + b + z2) / k + cos(b + z2) / k))
}

simuTest2.5 <- function(dat) {
    Y <- sort(unique(dat$Y))
    k <- runif(1, 1, 2)
    b <- runif(1, 1, 2)
    data.frame(Y = Y, z1 = k * Y * (2 * (Y > 5) - 1) + b, z2 = runif(1), k = k, b = b)
}
