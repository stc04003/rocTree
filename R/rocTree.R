#' Class definition
#' @importFrom methods getClass
#' @noRd 
setClass("splitClass",
         representation(tau = "numeric", M = "numeric", hN = "numeric", h = "numeric",
                        minsp = "numeric", minsp2 = "numeric", disc = "numeric", nflds = "numeric",
                        CV = "logical", Trace = "logical", parallel = "logical", parCluster = "numeric",
                        B = "numeric", ghN = "numeric", fsz = "numeric", splitBy = "character"),
         prototype(tau = 0.4, M = 1000, hN = 0, h = 0, minsp = 20, minsp2 = 5, disc = -1, nflds = 10,
                   CV = FALSE, Trace = FALSE, parallel = FALSE,
                   parCluster = parallel::detectCores() / 2,
                   B = 500, ghN = 0.2, fsz = 0, splitBy = "CON"),
         contains = "VIRTUAL")
setClass("CON", contains = "splitClass")
setClass("dCON", contains = "splitClass")
#' Method dispatch
#' @noRd
setGeneric("grow", function(Y, E, X.list, parm) standardGeneric("grow"))
setMethod("grow", signature(parm = "CON"), growSeq)
setMethod("grow", signature(parm = "dCON"), growRP)

## setGeneric("CV", function(Y1, E1, X1.list, Y2, E2, X2.list, X12.list, beta.seq, parm) standardGeneric("CV"))
## setMethod("CV", signature(parm = "CON"), CVSeq)
## setMethod("CV", signature(parm = "dCON"), CVRP)

#' ROC-guided Regression Trees
#'
#' Fits a "\code{rocTree}" model.
#'
#' The argument "control" defaults to a list with the following values:
#' \describe{
#'   \item{\code{tau}}{maximum follow-up time; default value is 0.4.}
#'   \item{\code{M}}{maximum node number allowed to be in the tree; the default value is 1000.}
#'   \item{\code{hN}}{smoothing parameter; the default value is "tau / 20".}
#'   \item{\code{minsp}, \code{minsp2}}{The interval (\code{minsp2}, \code{minsp}) denotes the
#' range for the number of the minimum risk observations required to split; the default value is (5, 20).}
#'   \item{\code{disc}}{a logical vector specifying whether the input covariate are discrete (\code{disc} = 0).
#' The length of "disc" should be the same as the number of covariates.}
#'   \item{\code{CV}}{a logical vector specifying whether a cross-validation is performed; the default value is FALSE.}
#'   \item{\code{nflds}}{number of folds; the default value is 10. This argument is only needed if \code{CV} = TRUE.}
#'   \item{\code{Trace}}{a logical vector specifying whether to display the splitting path; the default value is FALSE.}
#'   \item{\code{parallel}}{a logical vector specifying whether parallel computing is applied when \code{CV} = TRUE;
#' the default value is FALSE.}
#'   \item{\code{parCluster}}{an integer value specifying the number of CPU cores to be used when \code{CV} = TRUE and
#' \code{parallel} = TRUE. The default value is half of the CPU cores detected.}
#' }
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
#' @return An object of S3 class "\code{rocTree}" representing the fit, with the
#' following components:
#' \describe{
#' \item{Frame}{is a data frame describe the resulting tree.}
#' \item{r2Final}{estimated hazards at all terminal nodes.}
#' \item{tt}{time values to estimate hazard. Default length is 500. }
#' }
#' @references Sun Y. and Wang, M.C. (2018+). ROC-guided classification and survival trees. \emph{Technical report}.
#' @keywords rocTree
#' @seealso See \code{\link{print.rocTree}} and \code{\link{plot.rocTree}} for printing and plotting an \code{rocTree}, respectively.
#' @examples
#' library(survival)
#' set.seed(1)
#' dat <- simu(100, 0, 1.3)
#' fit <- rocTree(Surv(Y, death) ~ z1 + z2, id = id, data = dat,
#'        control = list(CV = TRUE, nflds = 5))
#' fit
rocTree <- function(formula, data, id, subset, splitBy = c("CON", "dCON"), control = list()) {
    splitBy <- match.arg(splitBy)
    parm.control <- control[names(control) %in% names(attr(getClass(splitBy), "slots"))]
    parm <- do.call("new", c(list(Class = splitBy), parm.control))
    parm@splitBy <- splitBy
    Call <- match.call()
    indx <- match(c("formula", "data", "id", "subset"), names(Call), nomatch = 0L)
    if (indx[1] == 0L) stop("A 'formula' argument is required")
    tmp <- Call[c(1L, indx)]
    tmp[[1L]] <- quote(stats::model.frame)
    ## prepare data
    m <- eval.parent(tmp)
    Y <- model.response(m)[,1]
    Status <- model.response(m)[,2]
    id <- model.extract(m, id)
    if (is.null(id)) id <- 1:nrow(m)
    Y0 <- unlist(lapply(split(Y, id), max))
    ## sort m by Y0
    ## do.call(rbind, split(m, id)[unique(id)[order(Y0)]])
    if (parm@h <= 0) parm@h <- parm@tau / 20
    if (parm@hN <= 0) parm@hN <- parm@tau / 20    
    if (parm@fsz <= 0) parm@fsz <- round(length(Y0) / 2)
    DF <- m[unlist(lapply(unique(id)[order(Y0)], function(x) which(id == x))),]
    rownames(DF) <- NULL
    Y <- model.response(DF)[,1]
    Status <- model.response(DF)[,2]
    id <- model.extract(DF, id)
    if (is.null(id)) id <- 1:nrow(m)
    Y0 <- unlist(lapply(split(Y, id), max))
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
    vNames <- colnames(X)
    if (parm@disc[1] < 0) parm@disc <- rep(0, p)
    if (length(parm@disc) == 1) parm@disc <- rep(parm@disc, p)
    xlist <- sapply(1:p, function(z) rocTree.Xlist(X[,z], parm@disc[z], Y, id), simplify = FALSE)
    xlist0 <- sapply(1:p, function(z) rocTree.Xlist(X[,z], 1, Y, id), simplify = FALSE) ## for prediction
    Y0 <- unlist(lapply(split(Y, id), max), use.names = FALSE)
    E0 <- unlist(lapply(split(Status, id), max), use.names = FALSE)
    ## Grow
    out <- grow(Y0, E0, xlist, parm)
    ## CV
    if (parm@CV & length(out$beta.seq) > 2) 
        out$con2.seq <- rocTree.cv(out$beta.seq, Y, Status, id, X, parm)
    out <- c(out, rocTree.final(out, Y0, E0, xlist, parm))
    ## out <- c(out, rocTree.haz(out$dfFinal, Y0, ctrl))
    ## names(out$r2Final) <- 
    names(out$dfFinal) <- paste("Node", out$ndFinal, sep = "")
    out$Y0 <- Y0
    out$E0 <- E0
    out$xlist <- xlist
    out$xlist0 <- xlist0
    out$vNames <- vNames
    ## Create Frame for print and plot
    ## Prepare Frame and remove nodes after considering ndFinal
    Frame <- out$treeMat
    Frame <- base::subset(Frame, !is.na(Frame$u))
    if (!is.null(out$ndFinal)) {
        Frame$terminal[which(Frame$nd %in% out$ndFinal)] <- 2
        if (sum((!(Frame$nd %in% out$ndFinal) & Frame$terminal == 2)) > 0) 
            Frame <- Frame[-which(!(Frame$nd %in% out$ndFinal) & Frame$terminal == 2),]
        tmp <- Frame$nd[Frame$terminal == 2]
        tmp2 <- tree.depth(tmp)
        rm1 <- c(sapply(out$nd[which(tmp2 < max(tmp2))],
                        function(a) unlist(sapply(1:(1 + max(tmp2)), function(b) (2^b * a):(2^b * (a + 1) - 1)))))
        rm <- unique(c(rm1, (max(tmp) + 1):(max(Frame$nd) + 1)))
        Frame <- Frame[Frame$nd %in% setdiff(Frame$nd, rm),]
    }
    out$Frame <- Frame
    out$parm <- parm
    out$terms <- attr(m, "terms")
    attr(out$terms, "id") <- Call[[match("id", names(Call))]]
    class(out) <- "rocTree"
    return(out)
}

#'
#' @param fsz forest size; S in the codes; removing this 
#' @noRd
rocTree.control <- function(l) {
    if (missing(l)) l <- NULL
    ## default list
    dl <- list(tau = 0.4, M = 1000, hN = NULL, h = NULL,
               minsp = 20, minsp2 = 5, disc = 0, nflds = 10, CV = FALSE, Trace = FALSE,
               parallel = FALSE, parCluster = detectCores() / 2, B = 500, ghN = .2,
               fsz = function(n) round(n/2))
    naml <- names(l)
    if (!all(naml %in% names(dl)))
        stop("unknown names in control: ", naml[!(naml %in% names(dl))])
    dl[naml] <- l
    if (is.null(dl$hN)) dl$hN <- dl$tau / 20
    if (is.null(dl$h)) dl$h <- dl$tau / 20
    return(dl)
}

#' Prepare xlist in \code{rocTree} and \code{rocForest}
#' 
#' An internal function used to prepare \code{x.list}.
#' This is exported for now to provide transformed covariate path
#' for fitting rpart and ctree.
#'
#' Should this be make public?
#'
#' @param x is a covariate vector.
#' @param disc a binary value indicating whether \code{x} is disctete;
#' \code{disc = 1} if discrete.
#' @param y is the ordered event time observed in the data.
#' @param id subject's id
#' 
#' @export
rocTree.Xlist <- function(x, disc, y, id) {
    yi <- unlist(lapply(split(y, id), max), use.names = FALSE)
    m <- unlist(lapply(split(y, id), length), use.names = FALSE)
    n <- length(unique(id))
    tmp <- unlist(lapply(split(y, id), function(z)
        match(yi, z))) + rep(c(0, cumsum(m)[-length(m)]), each = n)
    xlist <- matrix(x[tmp], n)
    if (!disc) xlist <- t(apply(xlist, 1, rank, ties.method = "max", na.last = "keep")) /
                   rowSums(!is.na(xlist))
    return(xlist)
}

is.rocTree <- function(x) inherits(x, "rocTree")

rocTree.cv <- function(beta.seq, Y, Status, id, X, parm) {
    n <- length(unique(id))
    nflds <- parm@nflds
    di <- unlist(lapply(split(Status, id), max), use.names = FALSE)
    flds <- folds(n, nflds)
    fldstep <- 0
    while(min(unlist(lapply(flds, function(x) sum(di[x])))) == 0) {
        flds <- folds(n, nflds)
        fldstep <- fldstep + 1
        if (fldstep > 50) stop("Not enough events")
    }
    p <- ncol(X)
    x0list <- sapply(1:p, function(z) rocTree.Xlist(X[,z], 1, Y, id), simplify = FALSE)
    if (!parm@parallel) con2.seq <- sapply(1:nflds, function(x) con.cv(Y, id, flds[[x]], Status, x0list, beta.seq, parm))
    if (parm@parallel) {
        cl <- makeCluster(parm@parCluster)
        clusterExport(cl = cl, 
                      varlist = c("Y", "id", "flds", "Status", "x0list", "beta.seq", "parm"),
                      envir = environment())
        con2.seq <- parSapply(cl, 1:nflds, function(x) con.cv(Y, id, flds[[x]], Status, x0list, beta.seq, parm))
        stopCluster(cl)
    }
    colnames(con2.seq) <- paste("result.", 1:nflds, sep = "")
    rownames(con2.seq) <- NULL
    con2.seq
}
