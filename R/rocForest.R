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
#' @param control a list of control parameters. See 'details' for important special
#' features of control parameters.
#' @export
#'
#' @examples
#' library(survival)
#' set.seed(123)
#' dat <- simu(40, 0, 1.1)
#' fit <- rocForest(Surv(Y, death) ~ z1 + z2, id = id, data = dat, control = list(minsp = 3, minsp2 = 1))
#' fit
#' 
#' @return An object of S3 class "\code{rocForest}" representing the fit, with the following components:
rocForest <- function(formula, data, id, subset, control = list()) {
    ctrl <- rocTree.control(control)
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
    if (is.null(ctrl$disc)) ctrl$disc <- rep(0, p)
    if (length(ctrl$disc) == 1) ctrl$disc <- rep(ctrl$disc, p)
    if (is.numeric(ctrl$fsz) && ctrl$fsz > length(Y)) stop("Invalid split size in forest.")
    xlist <- sapply(1:p, function(z) rocTree.Xlist(X[,z], ctrl$disc[z], Y, id), simplify = FALSE)
    xlist0 <- sapply(1:p, function(z) rocTree.Xlist(X[,z], 1, Y, id), simplify = FALSE) ## for prediction
    Y0 <- unlist(lapply(split(Y, id), max), use.names = FALSE)
    E0 <- unlist(lapply(split(Status, id), max), use.names = FALSE)
    out <- NULL
    out$Y0 <- Y0
    out$E0 <- E0
    out$xlist <- xlist
    out$xlist0 <- xlist0
    out$vNames <- vNames
    out$ctrl <- ctrl
    out$terms <- attr(m, "terms")
    attr(out$terms, "id") <- Call[[match("id", names(Call))]]
    if (!ctrl$parallel) out$forest <- lapply(1:ctrl$B, function(x) forest1(Y0, E0, xlist, ctrl))
    if (ctrl$parallel) {
        cl <- makeCluster(ctrl$parCluster)
        clusterExport(cl = cl, 
                      varlist = c("Y0", "E0", "xlist", "ctrl"),
                      envir = environment())
        out$forest <- parLapply(cl, 1:ctrl$B, function(x) forest1(Y0, E0, xlist, ctrl))
        stopCluster(cl)
    }
    class(out) <- "rocForest"
    return(out)
}

#' Function used to grow tree in forest
#'
#' This is an internal function, called by \code{rocForest}.
#' The function assumes a subsample is divided into two evenly-sized haves,
#' L1 and L2 (N1 ~ N / 4, N2 ~ N / 4). 
#' This function gives a tree constructed from L1, but with restrict by L2. 
#'
#' @param Y1 is a size N1 vector specifying the ordered follow-up time for L1.
#' @param E1 is the censoring indicator for the follow-up time \code{Y1} (\code{length(E) = length(Y)}).
#' @param X1.list covaraite path for the subjects associated with \code{Y1}.
#' @param Y2 is a size N2 vector specifying the ordered follow-up time for L2.
#' @param X2.list covaraite path for the subjects associated with \code{Y2}.
#' @param Y is the size N vector specifying the ordered follow-up time.
#' @param control see rocTree.control for default values.
#'
#' @return The function returns the \code{treeMat} in data frame format.
#' @keywords internal
#' @noRd
grow2 <- function(Y1, E1, X1.list, Y2, X2.list, Y, control) {
    tau <- control$tau
    M <- control$M
    hN <- control$hN
    minsp <- control$minsp
    minsp2 <- control$minsp2
    N <- length(Y)
    N1 <- length(Y1)
    N2 <- length(Y2)
    fit <- survfit(Surv(Y1, E1) ~ 1)
    fit$surv <- c(1, fit$surv)
    ss <- fit$surv[findInterval(Y1, fit$time)]
    Smat <- matk2 <- outer(ifelse(E1, Y1, NA), Y1, "<=") 
    fmat <- matk <- t(E1 * sapply(ifelse(E1, Y1, NA), K2, vec = Y1, h = hN) / hN)
    if (min(rowMeans(fmat)[E1 == 1]) < 0) 
        fmat[E1 == 1 & rowMeans(fmat) < 0, ] <- fmat[which(rowMeans(fmat) > 0)[1], ]
    Stau <- fit$surv[findInterval(tau, fit$time) + 1]
    EE1 <- Y1 <= tau * E1
    ss <- ss[EE1]
    fall <- rowMeans(fmat)
    const <- -diff(c(ss, Stau)) / (rowMeans(fmat)[EE1])
    ## Initialization
    fTree <- STree <- matrix(NA, M, sum(Y1 <= E1 * tau))
    fTree[1, ] <- rowMeans(fmat)[Y1 <= tau * E1]
    STree[1, ] <- rowMeans(Smat)[Y1 <= tau * E1]
    ## Define a tree object
    treeMat <- data.frame(nd = 1:M, terminal = 0, u = NA, u2 = NA, p = NA, cut = NA)
    treeMat[1, ] <- c(1, 1, 1, 1, NA, NA)
    ## node number of each observation
    ndInd <- matrix(1, N1, N1)
    ndInd[lower.tri(ndInd)] <- 0
    ndInd2 <- matrix(1, N, N2)
    ndInd2[outer(Y, Y2, FUN = ">")] <- 0
    while (sum(treeMat[, 2] == 1) > 0) {
        sp <- splitTree(X1.list, Y1, E1, fmat, Smat, treeMat, ndInd, const, fTree, STree, control, ceiling(sqrt(length(X1.list))))
        if (sp[1] * 2 < M & !is.na(sp[2])) {
            ndInd[ndInd == sp[1] & X1.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd[ndInd == sp[1] & X1.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            ## update L2 simutaneously
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            treeMat[sp[1], 2] <- 0
            treeMat[sp[1], 5:6] <- sp[2:3]
            treeMat[sp[1] * 2, ] <- c(sp[1] * 2, 1,
                                      mean(diag(ndInd) == sp[1] * 2 & EE1),
                                      min(rowMeans(ndInd[Y1 <= tau, ] == sp[1] * 2)), NA, NA)
            treeMat[sp[1] * 2 + 1, ] <- c(sp[1] * 2 + 1, 1,
                                          mean(diag(ndInd) == sp[1] * 2 + 1 & EE1),
                                          min(rowMeans(ndInd[Y1 <= tau, ] == sp[1] * 2 + 1)), NA, NA)
            treeMat$terminal[which(treeMat$u < minsp / N1 & treeMat$u2 < minsp2 / N1)] <- 2
            if(sum(diag(ndInd) == 2 * sp[1]) > 1) {
                fTree[sp[1] * 2, ] <- rowSums(fmat[EE1, diag(ndInd) == 2 * sp[1], drop = FALSE]) / N1
            } else {
                fTree[sp[1] * 2, ] <- sum(fmat[EE1, diag(ndInd) == 2 * sp[1]]) / N1
            }
            if(sum(diag(ndInd) == 2 * sp[1] + 1)>1) {
                fTree[sp[1] * 2 + 1, ] <- rowSums(fmat[EE1, diag(ndInd) == 2 * sp[1] + 1, drop = FALSE]) / N1
            } else {
                fTree[sp[1] * 2 + 1, ] <- sum(fmat[EE1, diag(ndInd) == 2 * sp[1] + 1]) / N1
            }
            STree[sp[1] * 2, ] <- rowSums(Smat * (ndInd == 2 * sp[1]))[EE1] / N1
            STree[sp[1] * 2 + 1, ] <- rowSums(Smat * (ndInd == 2 * sp[1] + 1))[EE1] / N1
        } else {
            ## if (sum(!is.na(treeMat$u2)) > 1) break
            treeMat$terminal[treeMat$nd == sp[1]] <- 2
            treeMat$terminal[treeMat$nd == 1] <- 0
            break
        }
    }
    treeMat <- treeMat[!is.na(treeMat$u),]
    ndTerm <- treeMat$nd[treeMat$terminal >= 1]
    if (nrow(treeMat) == 1) ndTerm <- 1
    szL2 <- sapply(ndTerm, function(x) rowSums(ndInd2 == x))
    list(treeMat = treeMat, szL2 = szL2, ndInd2 = ndInd2)
}

is.rocForest <- function(x) inherits(x, "rocForest")

#' This function provides one tree in \code{rocForest}
#'
#' @keywords internal
#' @noRd
forest1 <- function(Y, E, xlist, control) {
    n <- length(Y)
    ## S <- 2 * (1 + (n - 1) %/% 4)
    if (is.numeric(control$fsz)) S <- control$fsz
    else S <- control$fsz(n)
    idB <- sample(1:n, S)
    idB1 <- sort(idB[1:(S/2)])
    idB2 <- sort(idB[-(1:(S/2))])
    Y1 <- Y[idB1]
    E1 <- E[idB1]
    Y2 <- Y[idB2]
    X1.list <- lapply(xlist, function(x) x[idB1, idB1])
    X2.list <- lapply(xlist, function(x) x[,idB2])
    out <- grow2(Y1, E1, X1.list, Y2, X2.list, Y, control)
    out$idB2 <- idB2
    return(out)
}

