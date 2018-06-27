globalVariables(c("n", "cen")) ## global variables for simu

#' Function to generate simulated data used in the manuscript.
#'
#' This function is used to generate simulated data in the manuscript.
#' The underlying model is
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
#' \item{1.1}{assumes \eqn{\lambda(t, Z) = \lambda_0(t) e^{(2Z_1 + 2Z_2)}}.}
#' \item{1.2}{assumes \eqn{\lambda(t, Z) = \lambda_0(t) e^{[2sin(2\pi Z_1) + 2|Z_2 - 0.5|]}}.}
#' }
#' \item{Scenario 2}{assumes time dependent covariate (\eqn{Z_1}) with \eqn{\Lambda(0) = 2t}:}
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
#' @importFrom stats delete.response rexp rgamma rnorm runif
#' @export
#' 
simu <- function(n, cen, scenario) {
    if (!(cen %in% c(0, .25, .50)))
        stop("Only 3 levels of censoring rates (0%, 25%, 50%) are allowed.")
    if (!(scenario %in% paste(rep(1:3, each = 6), rep(1:6, 3), sep = ".")))
        stop("See ?simu for scenario definition")
    eval(parse(text = paste("simu", scenario, "(n = ", n, ", cen = ", cen, ")", sep = "")))
}

#' Function to generate the true hazard used in the simulation.
#'
#' This function is used to generate the true cumulative hazard function used in the simulation.
#'
#' @param tt is a numerical vector, specifying where to evaluate at.
#' @param scenario can be numeric or character string.
#' This indicates the simulation scenario noted in the manuscript. See \code{simu} for details.
#'
#' @export
trueHaz <- function(tt, scenario) {
    if (!(scenario %in% paste(rep(1:3, each = 6), rep(1:6, 3), sep = ".")))
        stop("See ?simu for scenario definition")
    eval(parse(text = paste("trueHaz", scenario, "(n = ", n, ", cen = ", cen, ")", sep = "")))    
}

#' Function to generate the true survival used in the simulation.
#'
#' This function is used to generate the true survival function used in the simulation.
#'
#' @param tt is a numerical vector, specifying where to evaluate at.
#' @param scenario can be numeric or character string.
#' This indicates the simulation scenario noted in the manuscript. See \code{simu} for details.
#'
#' @export
#' 
trueSurv <- function(tt, scenario) {
    if (!(scenario %in% paste(rep(1:3, each = 6), rep(1:6, 3), sep = ".")))
        stop("See ?simu for scenario definition")
    eval(parse(text = paste("trueSurv", scenario, "(n = ", n, ", cen = ", cen, ")", sep = "")))        
}

## Backgroud functions
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




##############################################################################################################################
## setClass("dataSetting",
##          representation(n = "numeric", cen = "numeric"),
##          prototype(n = 100))
## setGeneric("simu", function(dataSetting) standardGeneric("simu"))

