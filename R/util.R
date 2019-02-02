############################################################################################ 
## This file consists of utility functions used in the rocTree package
############################################################################################

#' Creates summaries for the final nodes.
#'
#' This function is called after split, prune, and CV (if applicable).
#'
#' @param tree is the tree-matrix; this matrix carries important information about the tree.
#' There are 6 columns in the tree. ("nd", "terminal", "u", "u2", "p", "cut").
#' @param Y0 is the ordered follow-up times.
#' @param E is the censoring indicator.
#' @param xlist is the covariate path.
#' @param control a list of control parameters.
#'
#' @return Returns a list of 2.
#' @param ndFinal node numbers after prune and cross-validation.
#' @param dfFinal is an n by nd matrix, where nd is the number of nodes in the final node numbers.
#' Each column gives the estimated unsmoothed \eqn{\hat{f}^{uc}_\tau = d\hat{F}^{uc}_\tau}. 
#' 
#' @keywords internal
#' @noRd
rocTree.final <- function(tree, Y0, E0, xlist, parm) {
    hN <- parm@hN
    n <- length(Y0)
    con2.seq <- tree$con2.seq
    tree2 <- data.frame(tree$treeMat)
    tree2$terminal[!is.na(tree2$p)] <- 3
    if (!is.null(con2.seq)) {
        tmp <- which.max(rank(rowMeans(con2.seq), ties.method = "first"))
        tree2$terminal[tree$optTree.seq[[tmp]]] <- 4
    } else {
        tree2$terminal[tree$optTree.seq[[1]]] <- 4
    }
    tree2 <- tree2[tree2$terminal >= 3,]
    ndInd <- 1 * !is.na(xlist[[1]])
    for (i in 1:NROW(tree2)) {
        sp <- tree2[i,]
        if (sp$terminal == 3) {
            ndInd[ndInd == sp$nd & xlist[[sp$p]] <= sp$cut] <- sp$nd * 2
            ndInd[ndInd == sp$nd & xlist[[sp$p]] > sp$cut] <- sp$nd * 2 + 1
        } 
    }
    fitc <- survfit(Surv(Y0, 1 - E0) ~ 1)
    sc <- with(summary(fitc), sapply(Y0, function(x) c(1, surv)[sum(x > c(0, time))]))
    tmp <- outer(Y0, Y0, "-")
    matk <- E0 * K1(tmp / hN) / hN
    matk2 <- tmp >= 0
    matk3 <- E0 * (tmp == 0)
    ndFinal <- tree2$nd[tree2$terminal == 4]
    dfFinal <- do.call(cbind, lapply(ndFinal, function(x)
        rowSums(t(matk3) * (ndInd == x)) / rowSums(t(matk2) * (ndInd == x))))
    dfFinal[!E0,] <- 0
    dfFinal[is.na(dfFinal)] <- 0
    list(ndFinal = ndFinal, dfFinal = data.frame(dfFinal))
}

#' Creates hazard estimates for the final nodes.
#'
#' This function is called after split, prune, and CV (if applicable).
#'
#' @param dfFinal is a n by nd matrix with each column consists unsmoothed \eqn{\hat{f}^{uc}_\tau = d\hat{F}^{uc}_\tau}.
#' @param Y0 is the ordered follow-up times.
#' @param parm a list of control parameters.
#' 
#' @keywords internal
#' @noRd
rocTree.haz <- function(dfFinal, Y0, parm) {
    hN <- parm@ghN
    tau <- parm@tau
    tt <- seq(0, tau, length = 500) ## move 500 to control later
    r2Final <- matrix(NA, 500, NCOL(dfFinal))
    for (i in 1:NCOL(dfFinal)) {
        r2Final[,i] <- colSums(sapply(tt, function(x) K3(x, Y0, hN)) * dfFinal[,i]) / hN
    }
    r2Final[r2Final < 0] <- 0
    list(tt = tt, r2Final = data.frame(r2Final))
}

#' Function to create k-fold for cross-validation
#'
#' \code{folds} creates k-fold for cross-validation.
#' This function is similar to "creadFolds" in the caret package.
#'
#' @param k is the number of folds
#' @param n is the length of Y; assuming to fold 1:n to k pieces.
#'
#' @keywords internal
#' @noRd
folds <- function(n, k) {
    if (n <= k) return(split(1:n, 1:n))
    if (k == 1) return(list(1:n))
    if (k > 1) {
        ind <- as.numeric(cut(1:n, breaks = min(5, k), include.lowest = TRUE))
        fold <- unlist(lapply(split(ind, ind), function(x) sample(c(rep(1:k, length(x) %/% k), sample(1:k, length(x) %% k)))))
        return(split(1:n, fold))
    }
}

#' Smoothing kernels used in the package
#'
#' @keywords internal
#' @noRd
K1 <- function(u) {
  0.75 * (1 - u ^ 2) * (abs(u) < 1)
}

#' This doesn't adjust for boundary condition and can cause errors
#' @keywords internal
#' @noRd
K2 <- function(s, vec, h) {
    if (is.na(s)) return(rep(NA, length(vec)))
    if (s < h) return(K1((h - vec) / h))
    ## if (s < h) return(Kq((s - vec) / h, s / h))
    else return(K1((s - vec) / h))
}

Kq <- function(x, q) {
    sigk1 <- sqrt(0.2)
    2 / (q + 1) * K1(2 / (q + 1) * (x - (q - 1) / 2)) *
        (1 + ((q - 1) / (q + 1) / sigk1) ^ 2 + 2 / sigk1 ^ 2 * (1 - q) / (1 + q) ^ 2 * x)
}

#' Smoothing kernels for hazard estimation
#'
#' @keywords internal
#' @noRd
K3 <- function(s, vec, h) {
    if (is.na(s)) return(rep(NA, length(vec)))
    if (s < h) return(Kq((s - vec) / h, s / h))
    else return(K1((s - vec) / h))
}
