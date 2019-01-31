
#' Function used to in the splitting procedure
#'
#' This is an internal function, called by \code{grow.Seq}. Inputs are based on an ordered Y.
#'
#' @section Variables
#' @param X is the covariate path (list). The length of the list is equation to the number of covaraites (p).
#' Each list is of size N by N, row indicates time (Yi), col indicates ID (i).
#' @param Y is the follow-up time, size N, consists of both censored and uncensored times.
#' @param E is the censoring indicator, size N.
#' @param fmat is a N by N matrix, the [i, j]th element gives \Delta_i K_n(Y_j - Y_i). rowMeans(fmat) gives \hat{f}^{UC}(Y_i).
#' @param Smat is a N by N matrix, the [i, j]th element gives (Y_j \le Y_i)\Delta. rowMeans(Smat) gives \hat{S}Y(Y_i) = \hat{P}(Y \ge Y_i).
#' @param treeMat tree-matrix. This is used to track what to split.
#' @param ndInd nodes to split. This is an indicator matrix, at each step, this points to what values to consider splitting.
#' @param const size = sum(E). Two ways to compute const:
#'   1. -diff(c(ss, Stau)) / (rowMeans(fmat)[Y <= tau * E])
#'   2. 1 / rowMeans(fmat)[Y <= tau * E] / sc[Y <= tau * E] / N
#'   This gives inverse the denominator in CON; f^{uc}_tau(t) P(Y \ge t, Z(t) \in \tau^\prime).
#' @param fTree is a M by sum(E) matrix. Use to track \hat{f}^{UC}(Y_i) at each step for i = 1, 2, ..., sum(E).
#' @param STree is a M by sum(E) matrix. Use to track \hat{S}Y(Y_i) at each step for i = 1, 2, ..., sum(E).
#' @param randP is a scalar value <= P. When \code{randP} is provided, only a random subset of p will be considered for splitting.
#' @param parm list of parameters prepared by control; see rocTree.control for default values
#'
#' @section Variables defined in this function:
#' @param nd.terminal potenial nodes to split; nd.terminal <- treeMat[treeMat[, 2] == 1, 1]
#'
#' @section At each nd.terminal: 
#' @param cutAll a vector contains all possible values to cut.
#' @param cutList a size p list. Each list is a \code{cutAll}.
#' @param dconList a size p list. Each list gives a vector of CON, corresponding to each cutAll value.
#' These values are set at -1 if length(cutAll) == 0 (nothing to cut).
#' @param dconmaxP a size p vector. Gives the max CON in each covariate.
#' @param r is the hazard estimates (f / S) for all terminal nodes that are not m in the nd.terminal loop.
#' @param rm is the hazard estimates (fm / Sm) for the mth terminal node that's being considering to split. 
#'
#' @section After checking all nd.terminal (set m in nd.terminal):
#' @param dconopt is a size M vector. Max CON at the mth check.
#' @param indm a scalar; "which covaraite gives the the max(CON) at the mth check (or dconopt[m]).
#' @param cutindm a scalar gives the value to cut at \code{indm}.
#' @param dconindm a scalar gives the CON if to cut at \code{indm}. This and \code{cutindm} are used to compute sopt. Not in outputs.
#' @param sopt is a M by 2 matrix. The first column gives "which covariate to split?"; the second column gives "Cut at where?".
#'
#' @section After all computation, here are the returned variables:
#' @param nd.split a scalar, indicating which node to split.
#' @param sopt2 a size 2 vector, the row in \code{sopt}.
#' @param dconopt2 a scalar, the corresponding CON.
#'
#' @section C codes;
#' In \code{cutAll} c codes, the concordance was computed with the property that
#' S_\tau(t) = S_{\tau_L}(t) + S_{\tau_R}(t) and f_\tau(t) = f_{\tau_L}(t) + f_{\tau_R}(t).
#' This is different then the \code{con} C code.
#'
#' @return returns a size 4 vector, (nd.split, sopt2, dconopt2) gives the information:
#' 1. which node to split
#' 2. which covariate to split
#' 3. given that covariate, where to cut? (split value)
#' 4. the corresponding CON
#' 
#' @keywords internal
#' @noRd 
splitSeq <- function(X, Y, E, fmat, Smat, treeMat, ndInd, const, fTree, STree, parm, randP = NULL)  {
    N <- dim(X[[1]])[1]
    P <- length(X)
    if (is.null(randP)) randP <- P
    disc <- parm@disc
    minsp <- parm@minsp
    minsp2 <- parm@minsp2
    M <- parm@M
    tau <- parm@tau
    ## all the terminal nodes that can be split; not all the terminal nodes
    nd.terminal <- treeMat$nd[treeMat$terminal == 1]
    sopt <- matrix(NA, M, 2)
    dconopt <- rep(0, M)
    ## lnd <- length(nd.terminal)-1
    lt <- sum(Y <= tau * E)
    for (m in nd.terminal) {
        ## need to change with discrete data
        dconList <- dconList2 <- list()
        cutList <- list()
        f <- fTree[(treeMat$terminal >= 1) & (treeMat$nd != m), ]
        S <- STree[(treeMat$terminal >= 1) & (treeMat$nd != m), ]
        ## size of nodes
        fm <- fTree[m, ]
        Sm <- STree[m, ]
        rm <- fm / Sm
        r <- f / S
        r[is.na(r)] <- Inf
        rm[is.na(rm)] <- Inf
        for (p in sample(1:P, randP)) {
            if (disc[p] == 0) {
                cutAll <- sort(unique(X[[p]][1, ndInd[1, ] == m]))
            } else {
                cutAll <- sort(unique(X[[p]][ndInd == m]))
            }
            cutAll <- cutAll[-length(cutAll)]
            ## if there is no potential cut off, skip
            ## need to add check l2 size
            if (length(cutAll) == 0) {
                dconList[[p]] <- -1
                next
            }
            cutList[[p]] <- cutAll
            dconList[[p]] <- .C("cutSearch",
                                as.integer(N), ## n
                                as.integer(length(cutAll)), ## cL
                                as.integer(m), ## m
                                as.integer(which(Y <= tau * E) - 1), ## y
                                as.integer(sum(Y <= tau * E)), ## Ny
                                as.double(minsp), ## minsp
                                as.double(minsp2), ## minsp2
                                as.double(ifelse(is.na(X[[p]]), 0, X[[p]])), ## X
                                as.double(ndInd), ## ndInd
                                as.double(cutAll), ## cut
                                as.double(ifelse(is.na(fmat), 0, fmat)), ## fmat
                                as.double(ifelse(is.na(Smat), 0, Smat)), ## Smat
                                as.double(const), ## const cL by 1
                                as.double(t(f)), ## length(nd.terminal) by Ny
                                as.double(t(S)), ## length(nd.terminal) by Ny
                                as.double(t(ifelse(r == Inf, 99999, r))), ## length(nd.terminal) by Ny
                                as.double(ifelse(rm == Inf , 99999, rm)), ## Ny by 1
                                as.integer(sum(treeMat$terminal >= 1 & treeMat$nd != m)), ## sum(treeMat[,2] >= 1 & treeMat[,1] != m)
                                out = double(length(cutAll)), PACKAGE = "rocTree")$out
        } ## end P
        ## dconmaxP <- unlist(lapply(dconList, max))
        ## print(m)
        ## print(lapply(dconList, summary))
        dconmaxP <- unlist(lapply(dconList, function(x) max(x, -1)))
        if (max(dconmaxP) < 0) {
            sopt[m, ] <- c(P + 1, 999)
            dconopt[m] <- -1
        } else {
            dconopt[m] <- max(dconmaxP)
            indm <- which.max(dconmaxP)
            cutindm <- cutList[[indm[1]]]
            dconindm <- dconList[[indm[1]]]
            sopt[m, ] <- c(indm[1], cutindm[which.max(dconindm)])
        }
    }
    nd.split <- which.max(dconopt)
    sopt2 <- sopt[nd.split, ]
    dconopt2 <- max(dconopt)
    c(nd.split, sopt2, dconopt2)
    ## which node to split, which variable, which cut off, which dcon
}

#' This is an internal function, called by \code{grow.RP}. Inputs are based on an ordered Y.
splitRP <- function(X, Y, E, fmat, Smat, treeMat, ndInd, const, fTree, STree, parm, randP = NULL)  {
    N <- dim(X[[1]])[1]
    P <- length(X)
    if (is.null(randP)) randP <- P
    disc <- parm@disc
    minsp <- parm@minsp
    minsp2 <- parm@minsp2
    M <- parm@M
    tau <- parm@tau
    ## all the terminal nodes that can be split; not all the terminal nodes
    nd.terminal <- treeMat$nd[treeMat$terminal == 1]
    sopt <- matrix(NA, M, 2)
    dconopt <- rep(0, M)
    ## lnd <- length(nd.terminal)-1
    lt <- sum(Y <= tau * E)
    for (m in nd.terminal) {
        ## need to change with discrete data
        dconList <- dconList2 <- list()
        cutList <- list()
        f <- fTree[(treeMat$terminal >= 1) & (treeMat$nd != m), ]
        S <- STree[(treeMat$terminal >= 1) & (treeMat$nd != m), ]
        ## size of nodes
        fm <- fTree[m, ]
        Sm <- STree[m, ]
        rm <- fm / Sm
        r <- f / S
        r[is.na(r)] <- Inf
        rm[is.na(rm)] <- Inf
        if (length(dconList2) < m || is.null(dconList2[[m]])) {
            dconList2[[m]] <- list()
            for (p in sample(1:P, randP)) {
                if (disc[p] == 0) {
                    cutAll <- sort(unique(X[[p]][1, ndInd[1, ] == m]))
                } else {
                    cutAll <- sort(unique(X[[p]][ndInd == m]))
                }
                cutAll <- cutAll[-length(cutAll)]
                ## if there is no potential cut off, skip
                ## need to add check l2 size
                if (length(cutAll) == 0) {
                    dconList[[p]] <- -1
                    next
                }
                cutList[[p]] <- cutAll
                dconList2[[m]][[p]] <- .C("cutSearch2", as.integer(N), as.integer(length(cutAll)), as.integer(m),
                                          as.integer(which(Y <= tau * E) - 1), as.integer(sum(Y <= tau * E)), 
                                          as.double(minsp), as.double(minsp2), 
                                          as.double(ifelse(is.na(X[[p]]), 0, X[[p]])),
                                          as.double(ndInd), as.double(cutAll),
                                          as.double(ifelse(is.na(fmat), 0, fmat)), 
                                          as.double(ifelse(is.na(Smat), 0, Smat)), 
                                          as.double(const), as.double(t(f)), as.double(t(S)), 
                                          as.double(t(ifelse(r == Inf, 99999, r))),
                                          as.double(ifelse(rm == Inf, 99999, rm)), 
                                          as.integer(sum(treeMat$terminal >= 1 & treeMat$nd != m)),
                                          out = double(length(cutAll)), PACKAGE = "rocTree")$out
            }
        }
        dconList <- dconList2[[m]]
        ## dconmaxP <- unlist(lapply(dconList, max))
        ## print(m)
        ## print(lapply(dconList, summary))
        dconmaxP <- unlist(lapply(dconList, function(x) max(x, -1)))
        if (max(dconmaxP) < 0) {
            sopt[m, ] <- c(P + 1, 999)
            dconopt[m] <- -1
        } else {
            dconopt[m] <- max(dconmaxP)
            indm <- which.max(dconmaxP)
            cutindm <- cutList[[indm[1]]]
            dconindm <- dconList[[indm[1]]]
            sopt[m, ] <- c(indm[1], cutindm[which.max(dconindm)])
        }
    }
    nd.split <- which.max(dconopt)
    sopt2 <- sopt[nd.split, ]
    dconopt2 <- max(dconopt)
    c(nd.split, sopt2, dconopt2)
    ## which node to split, which variable, which cut off, which dcon
}
