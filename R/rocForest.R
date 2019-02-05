#' ROC-guided Regression Forest
#'
#' Fits a "\code{rocForest}" model.
#' 
#' @param formula a formula object, with the response on the left of a '~' operator,
#' and the terms on the right. The response must be a survival object returned by the
#' 'Surv' function.
#' @param data an optional data frame in which to interpret the variables occurring
#' in the 'formula'.
#' @param id an optional vector used to identify time dependent covariate.
#' If missing, then each individual row of 'data' is presumed to represent a distinct
#' subject and each covariate is treated as a baseline covariate. The length of 'id'
#' should be the same as the number of observations.
#' @param subset an optional vector specifying a subset of observations to be used in
#' the fitting process.
#' @param splitBy a character string specifying the splitting algorithm. See *Details*.
#' @param control a list of control parameters. See 'details' for important special
#' features of control parameters.
#' @export
#'
#' @examples
#' library(survival)
#' set.seed(123)
#' dat <- simu(100, 0, 1.1)
#' fit <- rocForest(Surv(Y, death) ~ z1 + z2, id = id, data = dat,
#'         control = list(minsp = 3, minsp2 = 1, B = 50))
#' fit
#'
#' ## Print individual trees
#' print(fit, 1)
#' print(fit, 2)
#'
#' @return An object of S3 class "\code{rocForest}" representing the fit,
#' with the following components:
rocForest <- function(formula, data, id, subset, splitBy = c("CON", "dCON"), control = list()) {
    splitBy <- match.arg(splitBy)
    ## ctrl <- rocTree.control(control)
    Call <- match.call()
    indx <- match(c("formula", "data", "id", "subset"), names(Call), nomatch = 0L)
    if (indx[1] == 0L) stop("A 'formula' argument is required")
    tmp <- Call[c(1L, indx)]
    tmp[[1L]] <- quote(stats::model.frame)
    m <- eval.parent(tmp)
    Y <- model.response(m)[,1]
    Status <- model.response(m)[,2]
    id <- model.extract(m, id)
    if (is.null(id)) id <- 1:nrow(m)
    Y0 <- unlist(lapply(split(Y, id), max))
    DF <- m[unlist(lapply(unique(id)[order(Y0)], function(x) which(id == x))),]
    rownames(DF) <- NULL
    Y <- model.response(DF)[,1]
    Status <- model.response(DF)[,2]
    id <- model.extract(DF, id)
    if (is.null(id)) id <- 1:nrow(m)
    Y0 <- unlist(lapply(split(Y, id), max))
    if (!is.null(control$fsz) && is.function(control$fsz)) control$fsz <- control$fsz(length(Y0))
    if (is.null(control$fsz)) control$fsz <- round(length(Y0) / 2)
    if (is.numeric(control$fsz) && control$fsz > length(Y)) stop("Invalid split size in forest.")
    parm.control <- control[names(control) %in% names(attr(getClass(splitBy), "slots"))]
    parm <- do.call("new", c(list(Class = splitBy), parm.control))
    parm@splitBy <- splitBy
    if (parm@h <= 0) parm@h <- parm@tau/20
    if (parm@hN <= 0) parm@hN <- parm@tau/20
    if (!any(Y %in% Y0)) {
        DF <- DF[Y %in% Y0,]
        Y <- model.response(DF)[,1]
        id <- model.extract(DF, id)
    }
    if (is.null(id)) {id <- 1:nrow(DF)
    } else {id <- rep(1:length(unique(id)), table(factor(id, levels = unique(id))))}
    X <- model.matrix(attr(m, "terms"), DF)
    if (sum(colnames(X) == "(Intercept)") > 0) 
        X <- as.matrix(X[, -which(colnames(X) == "(Intercept)")])
    p <- ncol(X)
    if (parm@disc[1] < 0) parm@disc <- rep(0, p)
    if (length(parm@disc) == 1) parm@disc <- rep(parm@disc, p)    
    vNames <- colnames(X)
    xlist <- sapply(1:p, function(z) rocTree.Xlist(X[,z], parm@disc[z], Y, id), simplify = FALSE)
    xlist0 <- sapply(1:p, function(z) rocTree.Xlist(X[,z], 1, Y, id), simplify = FALSE) ## for prediction
    Y0 <- unlist(lapply(split(Y, id), max), use.names = FALSE)
    E0 <- unlist(lapply(split(Status, id), max), use.names = FALSE)
    out <- NULL
    out$Y0 <- Y0
    out$E0 <- E0
    out$xlist <- xlist
    out$xlist0 <- xlist0
    out$vNames <- vNames
    out$parm <- parm
    out$terms <- attr(m, "terms")
    out$call <- match.call()
    attr(out$terms, "id") <- Call[[match("id", names(Call))]]
    if (!parm@parallel) out$forest <- lapply(1:parm@B, function(x) forest1(Y0, E0, xlist, parm))
    if (parm@parallel) {
        cl <- makeCluster(parm@parCluster)
        clusterExport(cl = cl, 
                      varlist = c("Y0", "E0", "xlist", "parm"),
                      envir = environment())
        out$forest <- parLapply(cl, 1:parm@B, function(x) forest1(Y0, E0, xlist, parm))
        stopCluster(cl)
    }
    class(out) <- "rocForest"
    return(out)
}

is.rocForest <- function(x) inherits(x, "rocForest")

#' This function provides one tree in \code{rocForest}
#'
#' @keywords internal
#' @noRd
forest1 <- function(Y, E, xlist, parm) {
    n <- length(Y)
    ## S <- 2 * (1 + (n - 1) %/% 4)
    if (is.numeric(parm@fsz)) S <- parm@fsz
    else S <- parm@fsz(n)
    idB <- sample(1:n, S)
    idB1 <- sort(idB[1:(S/2)])
    idB2 <- sort(idB[-(1:(S/2))])
    Y1 <- Y[idB1]
    E1 <- E[idB1]
    Y2 <- Y[idB2]
    X1.list <- lapply(xlist, function(x) x[idB1, idB1])
    X2.list <- lapply(xlist, function(x) x[,idB2])
    if (parm@splitBy == "CON") out <- growSeq2(Y1, E1, X1.list, Y2, X2.list, Y, parm)
    if (parm@splitBy == "dCON") out <- growRP2(Y1, E1, X1.list, Y2, X2.list, Y, parm)
    out$idB2 <- idB2
    return(out)
}

