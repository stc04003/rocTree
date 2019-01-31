#' Function used to grow tree
#'
#' This is an internal function, called by \code{rocTree}.
#' Inputs are based on an ordered Y.
#'
#' @param Y is a size N vector specifying the ordered follow-up time.
#' @param E is the censoring indicator, \code{length(E) = length(Y)}.
#' @param X.list covaraite path, created by \code{rocTree.Xlist}.
#' @param control see rocTree.control for default values.
#'
#' Variables defined in this function:
#' @param const size = sum(E). Two ways to compute const:
#' 1. -diff(c(ss, Stau)) / (rowMeans(fmat)[Y <= tau * E])
#' 2. 1 / rowMeans(fmat)[Y <= tau * E] / sc[Y <= tau * E] / N
#' This gives the inverse of the denominator in CON; f^{uc}_tau(t) P(Y \ge t, Z(t) \in \tau^\prime).
#' 
#' Prepration and split stage:
#' @param fmat is a N by N matrix, the [i, j]th element gives \Delta_i K_n(Y_j - Y_i). rowMeans(fmat) gives \hat{f}^{UC}(Y_i).
#' @param Smat is a N by N matrix, the [i, j]th element gives (Y_j \le Y_i)\Delta. rowMeans(Smat) gives \hat{S}(Y_i) = \hat{P}(Y \ge Y_i).
#' @param fTree is a M by (E = 1) matrix, estimated density function (\hat{f}^{UC}(Y_i)) at each node.
#' @param STree is a M by (E = 1) matrix, estimated survival function (\hat{S}(Y_i) at each node.
#' @param sc estimated S_c(Y).
#' @param ss estimated S(Y).
#' @param Stau S(tau).
#' @param treeMat tree-matrix; 6 columns. ("nd", "terminal", "u", "u2", "p", "cut")
#' \describe{
#' \item{nd}{is the node number}
#' \item{terminal}{node property: 0 = internal, 1 = terminal node that can be split, 2 = terminal node that can't be split.
#' A node is classified as terminal = 2 is there are too few cases}
#' \item{u}{ is the proportion being split in \tau_L.}
#' \item{u2}{ is the minimum proporiton in \tau_L across time.}
#' \item{p}{ is the covariate being split at the node.}
#' \item{cut}{ is the covariate value being split at the node.}
#' }
#'
#' Prune stage:
#' @param treeMatTerm gives all potential nodes id that is splitted.
#' @param left The left most node in the bottom of the tree (assume tree grows downward with root on top).
#' @param right The right most node in the bottom of the tree.
#' @param treeMatTerm2 ordered treeMatTerm
#' @param ndTerm terminal nodes, constructed from the above arguments.
#'        This gives priority ordering on which node to prune first.
#' @param fuTerm a tempory vector, gives the estimated density function at a given (ndTerm) node.
#' @param SuTerm a tempory vector, gives the estimated survival function at a given (ndTerm) node.
#' @param sizeTree is the number of terminal nodes.
#' Lists used in prune stage: 
#' @param indList terminal code list for prune steps (0 = internal, 1, = terminal).
#' @param conList CON list for prune steps.
#' @param termList terminal node list for prune steps.
#' @param optTreeList optimum tree list. This list consists of all possible tree nodes after prune.
#' It's possible to have multiple rows in each list, indicating different ways to prune on the same node.
#' 
#' @return Return a list. Variables are here:
#' @param alpha is the penalty, defined as (CON[j] - CON[-(1:j)]) / #Tree.
#' @param beta.seq beta
#' @param optTree.seq optimum tree in sequence?
#' @param treeMat tree matrix for the fully grew tree.
#' 
#' @keywords internal
#' @noRd 
growRP <- function(Y, E, X.list, parm) {
    tau <- parm@tau
    M <- parm@M
    hN <- parm@hN
    minsp <- parm@minsp
    minsp2 <- parm@minsp2
    N <- length(Y)
    fitc <- survfit(Surv(Y, 1 - E) ~ 1)
    fit <- survfit(Surv(Y, E) ~ 1)
    fit$surv <- c(1, fit$surv)
    fitc$surv <- c(1, fitc$surv)
    sc <- fitc$surv[findInterval(Y, fitc$time)]
    ss <- fit$surv[findInterval(Y, fit$time)]
    Smat <- matk2 <- outer(ifelse(E, Y, NA), Y, "<=") 
    fmat <- matk <- t(E * sapply(ifelse(E, Y, NA), K2, vec = Y, h = hN) / hN)
    if (min(rowMeans(fmat)[E == 1]) < 0) 
        fmat[E == 1 & rowMeans(fmat) < 0, ] <- fmat[which(rowMeans(fmat) > 0)[1], ]
    Stau <- fit$surv[findInterval(tau, fit$time) + 1]
    EE <- Y <= tau * E
    ss <- ss[EE]
    fall <- rowMeans(fmat)
    const <- -diff(c(ss, Stau)) / (rowMeans(fmat)[EE])
    ## Initialization
    fTree <- STree <- matrix(NA, M, sum(Y <= E * tau))
    fTree[1, ] <- rowMeans(fmat)[EE]
    STree[1, ] <- rowMeans(Smat)[EE]
    ## Define a tree object
    treeMat <- data.frame(nd = 1:M, terminal = 0, u = NA, u2 = NA, p = NA, cut = NA)
    treeMat[1, ] <- c(1, 1, 1, 1, NA, NA)
    ## node number of each observation
    ndInd <- matrix(1, N, N)
    ndInd[lower.tri(ndInd)] <- 0
    ## ndInd is N N matrix
    conTree <- sum(0.5 * const * ss * rowMeans(fmat)[EE])
    while (sum(treeMat$terminal == 1) > 0) {
        sp <- splitTree(X.list, Y, E, fmat, Smat, treeMat, ndInd, const, fTree, STree, parm)
        if (sp[1] * 2 < M & !is.na(sp[2])) {
            ndInd[ndInd == sp[1] & X.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd[ndInd == sp[1] & X.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            fTree[sp[1] * 2, ] <- rowSums(fmat[EE, diag(ndInd) == 2 * sp[1], drop = FALSE]) / N
            fTree[sp[1] * 2 + 1, ] <- rowSums(fmat[EE, diag(ndInd) == 2 * sp[1] + 1, drop = FALSE]) / N
            STree[sp[1] * 2, ] <- rowSums(Smat * (ndInd == 2 * sp[1]))[EE] / N
            STree[sp[1] * 2 + 1, ] <- rowSums(Smat * (ndInd == 2 * sp[1] + 1))[EE] / N
            treeMat[sp[1], 2] <- 0
            treeMat[sp[1], 5:6] <- sp[2:3]
            treeMat[sp[1] * 2, ] <- c(sp[1] * 2, 1,
                                      mean(diag(ndInd) == sp[1] * 2 & EE),
                                      min(rowMeans(ndInd[Y <= tau, ] == sp[1] * 2)), NA, NA)
            treeMat[sp[1] * 2 + 1, ] <- c(sp[1] * 2 + 1, 1,
                                          mean(diag(ndInd) == sp[1] * 2 + 1 & EE),
                                          min(rowMeans(ndInd[Y <= tau, ] == sp[1] * 2 + 1)), NA, NA)
            treeMat$terminal[which(treeMat$u < minsp / N & treeMat$u2 < minsp2 / N)] <- 2
            conTree <- conTree + sp[4]
        } else {
            if (sum(!is.na(treeMat$p)) > 1) break
            else treeMat$terminal[treeMat$nd == sp[1]] <- 2
            break
        }
        if (parm@Trace) print(sp)
    }
    ## prune
    treeMatTerm <- treeMat$nd[treeMat$terminal >= 1 & is.na(treeMat$p)]
    ## sort the node
    left <- 2 ^ round(log(max(treeMatTerm), 2))
    right <- 2 * left - 1
    treeMatTerm2 <- treeMatTerm
    for (i in 1:length(treeMatTerm)) {
        nd <- treeMatTerm[i]
        treeMatTerm2[i] <- nd * 2 ^ ceiling(log(left / nd) / log(2))
    }
    ndTerm <- treeMatTerm <- treeMatTerm[order(treeMatTerm2)]
    ## exhaustive searching of all subtrees
    sizeTree <- length(ndTerm)
    termList <- list()
    indList <- list()
    conList <- list()
    termList[[1]] <- matrix(ndTerm, nrow = 1)
    indList[[1]] <- matrix(1, nrow = 1, ncol = sizeTree)
    fuTerm <- fTree[ndTerm[order(ndTerm)], ]
    SuTerm <- STree[ndTerm[order(ndTerm)], ]
    conList[[1]] <- .C("con",
                       as.integer(length(fuTerm) / sum(EE)),
                       as.integer(sum(EE)),
                       as.double(t(fuTerm)),
                       as.double(t(SuTerm)),
                       as.double(const), 
                       out = double(1), PACKAGE = "rocTree")$out
    for (k in 1:(sizeTree - 1)) {
        sizeNew <- sizeTree - k
        mat1 <- termList[[k]]
        mat2 <- indList[[k]]
        num <- 0
        termMat <- matrix(ncol = sizeNew, nrow = 1)
        indMat <- matrix(ncol = sizeNew, nrow = 1)
        consub <- c()
        for (l in 1:dim(mat1)[1]) {
            ndTerml <- mat1[l, ]
            indl <- mat2[l, ]
            ## change here
            locl <- which(diff(ndTerml) == 1 & ndTerml[-length(ndTerml)] %% 2 == 0 &
                          indl[-length(ndTerml)] == 1)
            for (q in locl) {
                ndNew <- ndTerml
                ndNew[q] <- ndNew[q] / 2
                ndNew <- ndNew[-(q + 1)]
                if (l == 1 & q == locl[1]) indNew <- rep(1, sizeNew)
                else indNew <- rep(0:1, c(q - 1, sizeNew - q + 1))
                termMat <- rbind(termMat, ndNew)
                indMat <- rbind(indMat, indNew)
                num <- num + 1
                fuTerm <- fTree[ndNew[order(ndNew)], ]
                SuTerm <- STree[ndNew[order(ndNew)], ]
                consub <- c(consub, .C("con",
                                       as.integer(length(fuTerm) / sum(EE)),
                                       as.integer(sum(EE)),
                                       as.double(t(fuTerm)),
                                       as.double(t(SuTerm)),
                                       as.double(const), 
                                       out = double(1), PACKAGE = "rocTree")$out)
            }
        } ## end for l
        termMat <- termMat[-1, ]
        indMat <- indMat[-1, ]
        if (is.vector(termMat)) {
            termList[[k + 1]] <- matrix(termMat, nrow = 1)
            indList[[k + 1]] <- matrix(indMat, nrow = 1)
        } else {
            termList[[k + 1]] <- termMat
            indList[[k + 1]] <- indMat
        }
        conList[[k + 1]] <- consub
    } ## end for k
    optTreeList <- list()
    con1 <- rep(NA, length(conList))
    for (k in 1:sizeTree) {
        con1[k] <- max(conList[[k]])
        optTreeList[[k]] <- termList[[k]][which(conList[[k]] == con1[k]), ]
    }
    sz1 <- sizeTree:1
    j <- 1
    res <- matrix(c(1, 0), nrow = 1)
    while (j < sizeTree) {
        alpha <- (con1[j] - con1[-(1:j)]) / (1:(sizeTree - j))
        j <- j + which.min(alpha)
        res <- rbind(res, c(j, min(alpha)))
    }
    ## beta.seq <- sqrt(c(res[, 2] * c(res[-1, 2], Inf)))
    beta.seq <- sqrt(abs(res[,2] * c(res[-1, 2], Inf)))
    beta.seq[length(beta.seq)] <- rev(res[,2])[1]
    list(beta.seq = beta.seq,
         optTree.seq = optTreeList[res[, 1]],
         treeMat = treeMat)
}

growSeq <- function(Y, E, X.list, parm) {
    tau <- parm@tau
    M <- parm@M
    hN <- parm@hN
    minsp <- parm@minsp
    minsp2 <- parm@minsp2
    N <- length(Y)
    fitc <- survfit(Surv(Y, 1 - E) ~ 1)
    fit <- survfit(Surv(Y, E) ~ 1)
    fit$surv <- c(1, fit$surv)
    fitc$surv <- c(1, fitc$surv)
    sc <- fitc$surv[findInterval(Y, fitc$time)]
    ss <- fit$surv[findInterval(Y, fit$time)]
    Smat <- matk2 <- outer(ifelse(E, Y, NA), Y, "<=") 
    fmat <- matk <- t(E * sapply(ifelse(E, Y, NA), K2, vec = Y, h = hN) / hN)
    if (min(rowMeans(fmat)[E == 1]) < 0) 
        fmat[E == 1 & rowMeans(fmat) < 0, ] <- fmat[which(rowMeans(fmat) > 0)[1], ]
    Stau <- fit$surv[findInterval(tau, fit$time) + 1]
    EE <- Y <= tau * E
    ss <- ss[EE]
    fall <- rowMeans(fmat)
    const <- -diff(c(ss, Stau)) / (rowMeans(fmat)[EE])
    ## Initialization
    fTree <- STree <- matrix(NA, M, sum(Y <= E * tau))
    fTree[1, ] <- rowMeans(fmat)[EE]
    STree[1, ] <- rowMeans(Smat)[EE]
    ## Define a tree object
    treeMat <- data.frame(nd = 1:M, terminal = 0, u = NA, u2 = NA, p = NA, cut = NA)
    treeMat[1, ] <- c(1, 1, 1, 1, NA, NA)
    ## node number of each observation
    ndInd <- matrix(1, N, N)
    ndInd[lower.tri(ndInd)] <- 0
    ## ndInd is N N matrix
    conTree <- sum(0.5 * const * ss * rowMeans(fmat)[EE])
    while (sum(treeMat$terminal == 1) > 0) {
        sp <- splitSeq(X.list, Y, E, fmat, Smat, treeMat, ndInd, const, fTree, STree, parm)
        if (sp[1] * 2 < M & !is.na(sp[2])) {
            ndInd[ndInd == sp[1] & X.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd[ndInd == sp[1] & X.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            fTree[sp[1] * 2, ] <- rowSums(fmat[EE, diag(ndInd) == 2 * sp[1], drop = FALSE]) / N
            fTree[sp[1] * 2 + 1, ] <- rowSums(fmat[EE, diag(ndInd) == 2 * sp[1] + 1, drop = FALSE]) / N
            STree[sp[1] * 2, ] <- rowSums(Smat * (ndInd == 2 * sp[1]))[EE] / N
            STree[sp[1] * 2 + 1, ] <- rowSums(Smat * (ndInd == 2 * sp[1] + 1))[EE] / N
            treeMat[sp[1], 2] <- 0
            treeMat[sp[1], 5:6] <- sp[2:3]
            treeMat[sp[1] * 2, ] <- c(sp[1] * 2, 1,
                                      mean(diag(ndInd) == sp[1] * 2 & EE),
                                      min(rowMeans(ndInd[Y <= tau, ] == sp[1] * 2)), NA, NA)
            treeMat[sp[1] * 2 + 1, ] <- c(sp[1] * 2 + 1, 1,
                                          mean(diag(ndInd) == sp[1] * 2 + 1 & EE),
                                          min(rowMeans(ndInd[Y <= tau, ] == sp[1] * 2 + 1)), NA, NA)
            treeMat$terminal[which(treeMat$u < minsp / N & treeMat$u2 < minsp2 / N)] <- 2
            conTree <- conTree + sp[4]
        } else {
            if (sum(!is.na(treeMat$p)) > 1) break
            else treeMat$terminal[treeMat$nd == sp[1]] <- 2
            break
        }
        if (parm@Trace) print(sp)
    }
    ## prune
    treeMatTerm <- treeMat$nd[treeMat$terminal >= 1 & is.na(treeMat$p)]
    ## sort the node
    left <- 2 ^ round(log(max(treeMatTerm), 2))
    right <- 2 * left - 1
    treeMatTerm2 <- treeMatTerm
    for (i in 1:length(treeMatTerm)) {
        nd <- treeMatTerm[i]
        treeMatTerm2[i] <- nd * 2 ^ ceiling(log(left / nd) / log(2))
    }
    ndTerm <- treeMatTerm <- treeMatTerm[order(treeMatTerm2)]
    ## exhaustive searching of all subtrees
    sizeTree <- length(ndTerm)
    termList <- list()
    indList <- list()
    conList <- list()
    termList[[1]] <- matrix(ndTerm, nrow = 1)
    indList[[1]] <- matrix(1, nrow = 1, ncol = sizeTree)
    fuTerm <- fTree[ndTerm[order(ndTerm)], ]
    SuTerm <- STree[ndTerm[order(ndTerm)], ]
    conList[[1]] <- .C("con",
                       as.integer(length(fuTerm) / sum(EE)),
                       as.integer(sum(EE)),
                       as.double(t(fuTerm)),
                       as.double(t(SuTerm)),
                       as.double(const), 
                       out = double(1), PACKAGE = "rocTree")$out
    for (k in 1:(sizeTree - 1)) {
        sizeNew <- sizeTree - k
        mat1 <- termList[[k]]
        mat2 <- indList[[k]]
        num <- 0
        termMat <- matrix(ncol = sizeNew, nrow = 1)
        indMat <- matrix(ncol = sizeNew, nrow = 1)
        consub <- c()
        for (l in 1:dim(mat1)[1]) {
            ndTerml <- mat1[l, ]
            indl <- mat2[l, ]
            ## change here
            locl <- which(diff(ndTerml) == 1 & ndTerml[-length(ndTerml)] %% 2 == 0 &
                          indl[-length(ndTerml)] == 1)
            for (q in locl) {
                ndNew <- ndTerml
                ndNew[q] <- ndNew[q] / 2
                ndNew <- ndNew[-(q + 1)]
                if (l == 1 & q == locl[1]) indNew <- rep(1, sizeNew)
                else indNew <- rep(0:1, c(q - 1, sizeNew - q + 1))
                termMat <- rbind(termMat, ndNew)
                indMat <- rbind(indMat, indNew)
                num <- num + 1
                fuTerm <- fTree[ndNew[order(ndNew)], ]
                SuTerm <- STree[ndNew[order(ndNew)], ]
                consub <- c(consub, .C("con",
                                       as.integer(length(fuTerm) / sum(EE)),
                                       as.integer(sum(EE)),
                                       as.double(t(fuTerm)),
                                       as.double(t(SuTerm)),
                                       as.double(const), 
                                       out = double(1), PACKAGE = "rocTree")$out)
            }
        } ## end for l
        termMat <- termMat[-1, ]
        indMat <- indMat[-1, ]
        if (is.vector(termMat)) {
            termList[[k + 1]] <- matrix(termMat, nrow = 1)
            indList[[k + 1]] <- matrix(indMat, nrow = 1)
        } else {
            termList[[k + 1]] <- termMat
            indList[[k + 1]] <- indMat
        }
        conList[[k + 1]] <- consub
    } ## end for k
    optTreeList <- list()
    con1 <- rep(NA, length(conList))
    for (k in 1:sizeTree) {
        con1[k] <- max(conList[[k]])
        optTreeList[[k]] <- termList[[k]][which(conList[[k]] == con1[k]), ]
    }
    sz1 <- sizeTree:1
    j <- 1
    res <- matrix(c(1, 0), nrow = 1)
    while (j < sizeTree) {
        alpha <- (con1[j] - con1[-(1:j)]) / (1:(sizeTree - j))
        j <- j + which.min(alpha)
        res <- rbind(res, c(j, min(alpha)))
    }
    ## beta.seq <- sqrt(c(res[, 2] * c(res[-1, 2], Inf)))
    beta.seq <- sqrt(abs(res[,2] * c(res[-1, 2], Inf)))
    beta.seq[length(beta.seq)] <- rev(res[,2])[1]
    list(beta.seq = beta.seq,
         optTree.seq = optTreeList[res[, 1]],
         treeMat = treeMat)
}
