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
#' @param control a list of control parameters. See 'details' for important special
#' features of control parameters.
#' @export
#' 
#' @return An object of S3 class "\code{rocTree}" representing the fit, with the
#' following components:
#' \describe{
#' \item{Frame}{is a data frame describe the resulting tree.}
#' }
#' @references Sun Y. and Wang, M.C. (2018+). ROC-guided classification and survival trees. \emph{Technical report}.
#' @keywords rocTree
#' @seealso See \code{\link{print.rocTree}} and \code{\link{plot.rocTree}} for printing and plotting an \code{rocTree}, respectively.
#' @examples
#' data(simudat)
#' library(survival)
#' system.time(fit <- rocTree(Surv(Time, Status) ~ X1 + X2 + X3, id = ID,
#' data = simudat, control = list(CV = TRUE, nflds = 10)))
#' fit
rocTree <- function(formula, data, id, subset, control = list()) {
    ctrl <- rocTree.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
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
    xlist <- sapply(1:p, function(z) rocTree.Xlist(X[,z], ctrl$disc[z], Y, id), simplify = FALSE)
    xlist0 <- sapply(1:p, function(z) rocTree.Xlist(X[,z], 1, Y, id), simplify = FALSE) ## for prediction
    Y0 <- unlist(lapply(split(Y, id), max), use.names = FALSE)
    E0 <- unlist(lapply(split(Status, id), max), use.names = FALSE)
    ## Grow
    out <- grow3(Y0, E0, xlist, ctrl)
    ## CV
    if (ctrl$CV & length(out$beta.seq) > 2) 
        out$con2.seq <- rocTree.cv(out$beta.seq, Y, Status, id, X, ctrl)
    out <- c(out, rocTree.final(out, Y0, E0, xlist, ctrl))
    out <- c(out, rocTree.haz(out$dfFinal, Y0, ctrl))
    names(out$dfFinal) <- names(out$r2Final) <- paste("Node", out$ndFinal, sep = "")
    out$Y0 <- Y0
    out$E0 <- E0
    out$xlist <- xlist
    out$xlist0 <- xlist0
    out$vNames <- vNames
    ## Create Frame for print and plot
    ## Prepare Frame and remove nodes after considering ndFinal
    Frame <- data.frame(out$treeMat)
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
    out$ctrl <- ctrl
    out$terms <- attr(m, "terms")
    attr(out$terms, "id") <- Call[[match("id", names(Call))]]
    class(out) <- "rocTree"
    return(out)
}

rocTree.control <- function(tau = 0.4, M = 1000, hN = tau / 20, h = hN,
                            minsp = 20, minsp2 = 5, disc = 0, nflds = 10, CV = FALSE, Trace = FALSE,
                            parallel = FALSE, parCluster = detectCores() / 2, ghN = 0.2,
                            split.method = c("CON", "dCON")) {
    list(tau = tau, M = M, hN = hN, h = h, minsp = minsp, minsp2 = minsp2, disc = disc,
         nflds = nflds, CV = CV, Trace = Trace, parallel = parallel, parCluster = parCluster, ghN = ghN,
         split.method = match.arg(split.method))
}

rocTree.Xlist <- function(x, disc, y, id) {
    yi <- unlist(lapply(split(y, id), max), use.names = FALSE)
    m <- unlist(lapply(split(y, id), length), use.names = FALSE)
    n <- length(unique(id))
    tmp <- unlist(lapply(split(y, id), function(z) match(yi, z))) + rep(c(0, cumsum(m)[-length(m)]), each = n)
    xlist <- matrix(x[tmp], n)
    if (!disc) xlist <- t(apply(xlist, 1, rank, ties.method = "max", na.last = "keep")) / rowSums(!is.na(xlist))
    return(xlist)
}

is.rocTree <- function(x) inherits(x, "rocTree")

## ---------------------------------------------------------------------------------------
## If CV == TURE
## ---------------------------------------------------------------------------------------

CV3 <- function(Y1, E1, X1.list, Y2, E2, X2.list, X12.list, beta.seq, control) {
    ## hN1 <- hN2 <- control$hN
    hN1 <- hN2 <- control$h * 1.2
    M <- control$M
    tau <- control$tau
    minsp <- control$minsp
    minsp2 <- control$minsp2
    N1 <- length(X1.list[[1]][1,])
    N2 <- length(X2.list[[1]][1,])
    fitc <- survfit(Surv(c(Y1, Y2), c(1 - E1, 1 - E2)) ~ 1)
    fit <- survfit(Surv(Y1, E1) ~ 1)
    fit$surv <- c(1, fit$surv)
    fitc$surv <- c(1, fitc$surv)
    sc <- fitc$surv[findInterval(Y1, fitc$time)]
    ss <- fit$surv[findInterval(Y1, fit$time)]
    Smat <- matk2 <- outer(ifelse(E1, Y1, NA), Y1, "<=") / matrix(sc, N1, N1, TRUE)
    fmat <- matk <- t(E1 * sapply(ifelse(E1, Y1, NA), K2, vec = Y1, h = hN1) / sc / hN1)
    if (min(rowMeans(fmat)[E1 == 1]) < 0) 
        fmat[E1 == 1 & rowMeans(fmat) < 0, ] <- fmat[which(rowMeans(fmat) > 0)[1], ]
    const <- 1 / rowMeans(fmat)[Y1 <= tau * E1] / sc[Y1 <= tau * E1] / N1
    sc2 <- fitc$surv[findInterval(Y2, fitc$time)]
    Smat2 <- matk2 <- outer(ifelse(E2, Y2, NA), Y2, "<=") / matrix(sc2, N2, N2, TRUE)
    fmat2 <- matk <- t(E2 * sapply(ifelse(E2, Y2, NA), K2, vec = Y2, h = hN2) / sc2 / hN2)
    if (min(rowMeans(fmat2)[E2 == 1]) < 0) 
        fmat[E2 == 1 & rowMeans(fmat2) < 0, ] <- fmat2[which(rowMeans(fmat2) > 0)[1], ]
    const2 <- 1 / rowMeans(fmat2)[Y2 <= tau * E2] / sc2[Y2 <= tau * E2] / N2
    Smat3 <- outer(ifelse(E2, Y2, NA), Y1, "<=") 
    fmat3 <- t(E1 * sapply(ifelse(E2, Y2, NA), K2, vec = Y1, h = hN1) / hN1)
    fTree <- STree <- matrix(NA, M, sum(Y1 <= E1 * tau))
    fTree2 <- STree2 <- fTree3 <- STree3 <- matrix(NA, M, sum(Y2 <= E2 * tau))
    fTree[1, ] <- rowMeans(fmat)[Y1 <= tau * E1]
    STree[1, ] <- rowMeans(Smat)[Y1 <= tau * E1]
    fTree2[1, ] <- rowMeans(fmat2)[Y2 <= tau * E2]
    STree2[1, ] <- rowMeans(Smat2)[Y2 <= tau * E2]
    fTree3[1, ] <- rowMeans(fmat3)[Y2 <= tau * E2]
    STree3[1, ] <- rowMeans(Smat3)[Y2 <= tau * E2]
    treeMat <- matrix(nrow = M, ncol = 6)
    colnames(treeMat) <- c("nd", "terminal", "u", "u2", "p", "cut")
    treeMat[, 1] <- 1:M
    treeMat[, 2] <- 0
    treeMat[1, ] <- c(1, 1, 1, 1, NA, NA)
    ndInd1 <- matrix(1, N1, N1)
    ndInd2 <- matrix(1, N2, N2)
    ndInd12 <- matrix(1, N2, N1)
    ndInd1[lower.tri(ndInd1)] <- 0
    ndInd2[lower.tri(ndInd2)] <- 0
    ndInd12[is.na(X12.list[[1]])] <- 0
    tlst <- which(Y1 >= tau)[1] - 1
    while (sum(treeMat[, 2] == 1) > 0) {
        sp <- split3(X1.list, Y1, E1, fmat, Smat, treeMat, ndInd1, const, fTree, STree, control)
        if (sp[1] * 2 < M & !is.na(sp[2])) {
            ndInd1[ndInd1 == sp[1] & X1.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd1[ndInd1 == sp[1] & X1.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            ndInd12[ndInd12 == sp[1] & X12.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd12[ndInd12 == sp[1] & X12.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            fTree[sp[1] * 2, ] <- rowSums(as.matrix(fmat[Y1 <= tau * E1, diag(ndInd1) == 2 * sp[1]])) / N1
            fTree[sp[1] * 2 + 1, ] <- rowSums(as.matrix(fmat[Y1 <= tau * E1, diag(ndInd1) == 2 * sp[1] + 1])) / N1
            STree[sp[1] * 2, ] <- rowSums(as.matrix(Smat * (ndInd1 == 2 * sp[1])))[Y1 <= tau * E1] / N1
            STree[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat * (ndInd1 == 2 * sp[1] + 1)))[Y1 <= tau * E1] / N1
            fTree2[sp[1] * 2, ] <- rowSums(as.matrix(fmat2[Y2 <= tau * E2, diag(ndInd2) == 2 * sp[1]])) / N2
            fTree2[sp[1] * 2 + 1, ] <- rowSums(as.matrix(fmat2[Y2 <= tau * E2, diag(ndInd2) == 2 * sp[1] + 1])) / N2
            STree2[sp[1] * 2, ] <- rowSums(as.matrix(Smat2 * (ndInd2 == 2 * sp[1])))[Y2 <= tau * E2] / N2
            STree2[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat2 * (ndInd2 == 2 * sp[1] + 1)))[Y2 <= tau * E2] / N2
            ## rTree2[sp[1]*2,] <- fTree2[sp[1]*2,]/STree2[sp[1]*2,]
            ## rTree2[sp[1]*2+1,] <- fTree2[sp[1]*2+1,]/STree2[sp[1]*2+1,]
            fTree3[sp[1] * 2, ] <- rowSums(as.matrix(fmat3[Y2 <= tau * E2, diag(ndInd1) == 2 * sp[1]])) / N1
            fTree3[sp[1] * 2 + 1, ] <- rowSums(as.matrix(fmat3[Y2 <= tau * E2, diag(ndInd1) == 2 * sp[1] + 1])) / N1
            STree3[sp[1] * 2, ] <- rowSums(as.matrix(Smat3 * (ndInd12 == 2 * sp[1])))[Y2 <= tau * E2] / N1
            STree3[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat3 * (ndInd12 == 2 * sp[1] + 1)))[Y2 <= tau * E2] / N1
            ## rTree3[sp[1]*2,] <- fTree3[sp[1]*2,]/STree3[sp[1]*2,]
            ## rTree3[sp[1]*2+1,] <- fTree3[sp[1]*2+1,]/STree3[sp[1]*2+1,]
            treeMat[sp[1], 2] <- 0
            treeMat[sp[1], 5:6] <- sp[2:3]
            treeMat[sp[1] * 2, ] <- c(sp[1] * 2, 1,
                                      mean(diag(ndInd1) == sp[1] * 2 & Y1 <= tau * E1),
                                      min(rowMeans(ndInd1[Y1 <= tau, ] == sp[1] * 2)), NA, NA)
            treeMat[sp[1] * 2 + 1, ] <- c(sp[1] * 2 + 1, 1,
                                          mean(diag(ndInd1) == (sp[1] * 2 + 1) & Y1 <= tau * E1),
                                          min(rowMeans(ndInd1[Y1 <= tau, ] == sp[1] * 2 + 1)), NA, NA)
            treeMat[treeMat[, 3] <= minsp / N1 & treeMat[, 4] <= minsp2 / N1, 2] <- 2
            ## conTree <- conTree + sp[4]
        } else {
            treeMat[treeMat[, 1] == sp[1], 2] <- 2
            break
        }
    }
    ## prune
    treeMatTerm <- treeMat[treeMat[, 2] >= 1 & is.na(treeMat[, 5]), 1]
    ## sort the node
    left <- 2 ^ round(max(log(treeMatTerm, 2)))
    right <- 2 * left - 1
    treeMatTerm2 <- treeMatTerm
    for (i in 1:length(treeMatTerm)) {
        nd <- treeMatTerm[i]
        treeMatTerm2[i] <- nd * 2 ^ ceiling(log(left / nd) / log(2))
    }
    treeMatTerm <- treeMatTerm[order(treeMatTerm2)]
    ## exhaustive searching of all subtrees
    ndTerm <- treeMatTerm
    sizeTree <- length(ndTerm)
    termList <- list()
    indList <- list()
    conList <- list()
    conList2 <- list()
    termList[[1]] <- matrix(ndTerm, nrow = 1)
    indList[[1]] <- matrix(1, nrow = 1, ncol = sizeTree)
    fuTerm <- fTree[ndTerm[order(ndTerm)], ]
    SuTerm <- STree[ndTerm[order(ndTerm)], ]
    fuTerm2 <- fTree2[ndTerm[order(ndTerm)], ]
    SuTerm2 <- STree2[ndTerm[order(ndTerm)], ]
    fuTerm3 <- fTree3[ndTerm[order(ndTerm)], ]
    SuTerm3 <- STree3[ndTerm[order(ndTerm)], ]
    conList[[1]] <- .C("con", 
                       as.integer(length(fuTerm) / sum(Y1 <= tau * E1)),
                       as.integer(sum(Y1 <= tau * E1)),
                       as.double(t(fuTerm)),
                       as.double(t(SuTerm)),
                       as.double(const), 
                       double(1))[[6]]
    conList2[[1]] <- .C("con3", 
                        as.integer(length(fuTerm2) / sum(Y2 <= tau * E2)),
                        as.integer(sum(Y2 <= tau * E2)),
                        as.double(t(fuTerm2)),
                        as.double(t(SuTerm2)),
                        as.double(t(fuTerm3)),
                        as.double(t(SuTerm3)),
                        as.double(const2),
                        double(1))[[8]]
    ## con3(fuTerm3, SuTerm3, fuTerm2, SuTerm2, const2)
    for (k in 1:(sizeTree - 1)) {
        sizeNew <- sizeTree - k
        mat1 <- termList[[k]]
        mat2 <- indList[[k]]
        num <- 0
        termMat <- matrix(ncol = sizeNew, nrow = 1)
        indMat <- matrix(ncol = sizeNew, nrow = 1)
        consub <- c()
        con2sub <- c()
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
                fuTerm <- fTree[treeMat[, 1] %in% ndNew, ]
                SuTerm <- STree[treeMat[, 1] %in% ndNew, ]
                fuTerm2 <- fTree2[treeMat[, 1] %in% ndNew, ]
                SuTerm2 <- STree2[treeMat[, 1] %in% ndNew, ]
                fuTerm3 <- fTree3[treeMat[, 1] %in% ndNew, ]
                SuTerm3 <- STree3[treeMat[, 1] %in% ndNew, ]
                consub <- c(consub, .C("con", 
                                       as.integer(length(fuTerm) / sum(Y1 <= tau * E1)),
                                       as.integer(sum(Y1 <= tau * E1)),
                                       as.double(t(fuTerm)),
                                       as.double(t(SuTerm)),
                                       as.double(const), 
                                       double(1))[[6]])
                con2sub <- c(con2sub, .C("con3", 
                                         as.integer(length(fuTerm2) / sum(Y2 <= tau * E2)),
                                         as.integer(sum(Y2 <= tau * E2)),
                                         as.double(t(fuTerm2)),
                                         as.double(t(SuTerm2)),
                                         as.double(t(fuTerm3)),
                                         as.double(t(SuTerm3)),
                                         as.double(const2),
                                         double(1))[[8]])
            }
        }
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
        conList2[[k + 1]] <- con2sub
    }
    con1 <- unlist(lapply(conList, max))
    sz1 <- sizeTree:1
    j <- 1
    res <- matrix(c(1, 0), nrow = 1)
    while (j < sizeTree) {
        alpha <- (con1[j] - con1[-(1:j)]) / (1:(sizeTree - j))
        j <- j + tail(which(alpha == min(alpha)), 1)
        res <- rbind(res, c(j, min(alpha)))
    }
    con2opt <- (unlist(lapply(conList2, max)))[res[, 1]]
    con1opt <- con1[res[, 1]]
    res <- cbind(res, con1opt, con2opt)
    res <- rbind(res, c(100, Inf, 0, 0))
    ## beta.seq <- c(0.0001,0.001,0.01,0.1)
    cvconbb <- function(bb) res[which(res[, 2] >= bb)[1] - 1, 4]
    cvcon2 <- sapply(beta.seq[-length(beta.seq)], cvconbb)
    c(unlist(cvcon2), tail(res[, 4], 2)[1])
}

con.cv <- function(y, id, ind, d, x0, beta.seq, control) {
    yi <- unlist(lapply(split(y, id), max), use.names = FALSE)
    di <- unlist(lapply(split(d, id), max), use.names = FALSE)
    n <- length(unique(id))
    p <- length(x0)
    y1 <- yi[-ind]
    y2 <- yi[ind]
    d1 <- di[-ind]
    d2 <- di[ind]
    x1 <- x2 <- x12 <- list()
    n1 <- length(y1)
    n2 <- length(y2)
    for (i in 1:p) {
        x1[[i]] <- x0[[i]][-ind, -ind]
        x2[[i]] <- x0[[i]][ind, ind]
        x12[[i]] <- x0[[i]][ind, -ind]
        if (!control$disc[i]) {
            x1tmp <- x1[[i]]
            x2tmp <- x12[[i]]
            x2tmp[rowSums(!is.na(x2tmp)) == 0,] <- 1
            x1[[i]] <- t(apply(x1tmp, 1, rank, ties.method = 'min', na.last = "keep")) / rowSums(!is.na(x1tmp))
            x2[[i]] <- t(sapply(1:n2, function(x) findInterval(x2[[i]][x,], sort(x2tmp[x,])) / sum(!is.na(x2tmp[x,]))))
            x12[[i]] <-t(sapply(1:n2, function(x) findInterval(x12[[i]][x,], sort(x2tmp[x,])) / sum(!is.na(x2tmp[x,]))))
        }
    }
    CV3(y1, d1, x1, y2, d2, x2, x12, beta.seq, control)
}

rocTree.cv <- function(beta.seq, Y, Status, id, X, control = list()) {
    n <- length(unique(id))
    nflds <- control$nflds
    flds <- folds(n, nflds)
    ## if (!(n %% nflds)) flds <- split(sample(1:n, n), rep(1:nflds, each = floor(n / nflds)))
    ## else flds <- split(sample(1:n, n), c(rep(1:nflds, each = floor(n / nflds)), 1:(n %% nflds)))
    ## flds <- lapply(flds, function(x) (1:n)[-x]) ## if returnTrain = TRUE
    ## flds <- createFolds(1:n, k = nflds, list = TRUE, returnTrain = FALSE)
    p <- ncol(X)
    x0list <- sapply(1:p, function(z) rocTree.Xlist(X[,z], 1, Y, id), simplify = FALSE)
    if (!control$parallel) con2.seq <- sapply(1:nflds, function(x) con.cv(Y, id, flds[[x]], Status, x0list, beta.seq, control))
    if (control$parallel) {
        cl <- makeCluster(control$parCluster)
        clusterExport(cl = cl, 
                      varlist = c("Y", "id", "flds", "Status", "x0list", "beta.seq", "control"),
                      envir = environment())
        con2.seq <- parSapply(cl, 1:nflds, function(x) con.cv(Y, id, flds[[x]], Status, x0list, beta.seq, control))
        stopCluster(cl)
    }
    colnames(con2.seq) <- paste("result.", 1:nflds, sep = "")
    rownames(con2.seq) <- NULL
    con2.seq
}

#' Function used to grow tree
#'
#' This is an internal function, called by \code{rocTree}.
#' Inputs are based on an ordered Y.
#'
#' @param Y is a size N vector specifying the ordered follow-up time.
#' @param E is the censoring indicator, \code{length(E) = length(Y)}.
#' @param X.list covaraite path, created by \code{rocTree.Xlist}.
#' @param control see rocTree.control for default values
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
grow3 <- function(Y, E, X.list, control) {
    tau <- control$tau
    M <- control$M
    hN <- control$hN
    minsp <- control$minsp
    minsp2 <- control$minsp2
    N <- length(Y)
    ## N <- length(X.list[[1]][1, ])
    fitc <- survfit(Surv(Y, 1 - E) ~ 1)
    fit <- survfit(Surv(Y, E) ~ 1)
    fit$surv <- c(1, fit$surv)
    fitc$surv <- c(1, fitc$surv)
    sc <- fitc$surv[findInterval(Y, fitc$time)]
    ss <- fit$surv[findInterval(Y, fit$time)]
    Smat <- matk2 <- outer(ifelse(E, Y, NA), Y, "<=") ## opposite of Y >= Y[i]
    fmat <- matk <- t(E * sapply(ifelse(E, Y, NA), K2, vec = Y, h = hN) / hN)
    if (min(rowMeans(fmat)[E == 1]) < 0) 
        fmat[E == 1 & rowMeans(fmat) < 0, ] <- fmat[which(rowMeans(fmat) > 0)[1], ]
    Stau <- fit$surv[findInterval(tau, fit$time) + 1]
    EE <- Y <= tau * E
    ss <- ss[EE]
    ## ss <- ss[Y <= tau & E == 1]
    fall <- rowMeans(fmat)
    const <- -diff(c(ss, Stau)) / (rowMeans(fmat)[EE])
    ## initialization
    fTree <- STree <- matrix(NA, M, sum(Y <= E * tau))
    fTree[1, ] <- rowMeans(fmat)[EE]
    STree[1, ] <- rowMeans(Smat)[EE]
    ## Define a tree object
    treeMat <- matrix(nrow = M, ncol = 6)
    colnames(treeMat) <- c("nd", "terminal", "u", "u2", "p", "cut")
    ## terminal = 0 - internal
    ## terminal = 1 - terminal can be split
    ## terminal = 2 - terminal cannot be split
    treeMat[, 1] <- 1:M
    treeMat[, 2] <- 0
    treeMat[1, ] <- c(1, 1, 1, 1, NA, NA)
    ## node number of each observation
    ndInd <- matrix(1, N, N)
    ndInd[lower.tri(ndInd)] <- 0
    ## ndInd is N N matrix
    conTree <- sum(0.5 * const * ss * rowMeans(fmat)[EE])
    while (sum(treeMat[, 2] == 1) > 0) {
        sp <- split3(X.list, Y, E, fmat, Smat, treeMat, ndInd, const, fTree, STree, control)
        if (sp[1] * 2 < M & !is.na(sp[2])) {
            ndInd[ndInd == sp[1] & X.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd[ndInd == sp[1] & X.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            fTree[sp[1] * 2, ] <- rowSums(fmat[EE, diag(ndInd) == 2 * sp[1]]) / N
            fTree[sp[1] * 2 + 1, ] <- rowSums(fmat[EE, diag(ndInd) == 2 * sp[1] + 1]) / N
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
            ## mean(ndInd[tlst,] == sp[1]*2)
            treeMat[treeMat[, 3] < minsp / N & treeMat[, 4] < minsp2 / N, 2] <- 2
            conTree <- conTree + sp[4]
        } else {
            treeMat[treeMat[, 1] == sp[1], 2] <- 2
            break
        }
        if (control$Trace) print(sp)
    }
    ## prune
    treeMatTerm <- treeMat[treeMat[, 2] >= 1 & is.na(treeMat[, 5]), 1]
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
    list(beta.seq = beta.seq,
         optTree.seq = optTreeList[res[, 1]],
         treeMat = treeMat)
}

#' Function used to in the splitting procedure
#'
#' This is an internal function, called by \code{grow}. Inputs are based on an ordered Y.
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
#' @param control see rocTree.control for default values
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
split3 <- function(X, Y, E, fmat, Smat, treeMat, ndInd, const, fTree, STree, control)  {
    N <- dim(X[[1]])[1]
    P <- length(X)
    disc <- control$disc
    minsp <- control$minsp
    minsp2 <- control$minsp2
    M <- control$M
    tau <- control$tau
    ## all the terminal nodes that can be split; not all the terminal nodes
    nd.terminal <- treeMat[treeMat[, 2] == 1, 1]
    sopt <- matrix(NA, M, 2)
    dconopt <- rep(0, M)
    ## lnd <- length(nd.terminal)-1
    lt <- sum(Y <= tau * E)
    for (m in nd.terminal) {
        ## need to change with discrete data
        dconList <- list()
        cutList <- list()
        f <- fTree[(treeMat[, 2] >= 1) & (treeMat[, 1] != m), ]
        S <- STree[(treeMat[, 2] >= 1) & (treeMat[, 1] != m), ]
        ## size of nodes
        fm <- fTree[m, ]
        Sm <- STree[m, ]
        rm <- fm / Sm
        r <- f / S
        r[is.na(r)] <- Inf
        rm[is.na(r)] <- Inf
        for (p in 1:P) {
            ## if discrete
            if (disc[p] == 0) {
                cutAll <- sort(unique(X[[p]][1, ndInd[1, ] == m]))
            } else {
                cutAll <- sort(unique(X[[p]][ndInd == m]))
            }
            cutAll <- cutAll[-length(cutAll)]
            ## if there is no potential cut off, skip
            if (length(cutAll) == 0) {
                dconList[[p]] <- -1
                next
            }
            cutList[[p]] <- cutAll
            if (control$split.method == "CON") 
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
                                    as.double(ifelse(rm == Inf | is.na(rm), 99999, rm)), ## Ny by 1
                                    as.integer(sum(treeMat[,2] >= 1 & treeMat[,1] != m)), ## sum(treeMat[,2] >= 1 & treeMat[,1] != m)
                                    out = double(length(cutAll)), PACKAGE = "rocTree")$out
            if (control$split.method == "dCON") 
                dconList[[p]] <- .C("cutSearch2", as.integer(N), as.integer(length(cutAll)), as.integer(m),
                                    as.integer(which(Y <= tau * E) - 1), as.integer(sum(Y <= tau * E)), 
                                    as.double(minsp), as.double(minsp2), 
                                    as.double(ifelse(is.na(X[[p]]), 0, X[[p]])),
                                    as.double(ndInd), as.double(cutAll),
                                    as.double(ifelse(is.na(fmat), 0, fmat)), 
                                    as.double(ifelse(is.na(Smat), 0, Smat)), 
                                    as.double(const), as.double(t(f)), as.double(t(S)), 
                                    as.double(t(ifelse(r == Inf, 99999, r))),
                                    as.double(ifelse(rm == Inf | is.na(rm), 99999, rm)), 
                                    as.integer(sum(treeMat[,2] >= 1 & treeMat[,1] != m)),
                                    out = double(length(cutAll)), PACKAGE = "rocTree")$out
        } ## end P
        dconmaxP <- unlist(lapply(dconList, max))
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

## ## ---------------------------------------------------------------------------------------
## ## concordance
## ## ---------------------------------------------------------------------------------------

## con <- function(fuTerm, SuTerm, const) {
##     if (is.vector(fuTerm)) {
##         fuTerm <- t(as.matrix(fuTerm))
##         SuTerm <- t(as.matrix(SuTerm))
##     }
##     mnd <- dim(fuTerm)[1]
##     ind1 <- matrix(rep(1:mnd, mnd), mnd, mnd)
##     ind2 <- matrix(rep(1:mnd, each = mnd), mnd, mnd)
##     matfs1 <- fuTerm[as.vector(ind1), ] * SuTerm[as.vector(ind2), ]
##     matfs1 <- t(t(matfs1) * const)
##     matfs2 <- fuTerm[as.vector(ind2), ] * SuTerm[as.vector(ind1), ]
##     matfs2 <- t(t(matfs2) * const)
##     sum(matfs1[matfs1 > matfs2]) + 0.5 * sum(matfs1[matfs1 == matfs2])
## }


## con3 <- function(fuTerm3, SuTerm3, fuTerm2, SuTerm2, const2) {
##     if (is.vector(fuTerm3)) {
##         fuTerm3 <- t(as.matrix(fuTerm3))
##         SuTerm3 <- t(as.matrix(SuTerm3))
##         fuTerm2 <- t(as.matrix(fuTerm2))
##         SuTerm2 <- t(as.matrix(SuTerm2))
##     }
##     mnd <- dim(fuTerm3)[1]
##     ind1 <- matrix(rep(1:mnd, mnd), mnd, mnd)
##     ind2 <- matrix(rep(1:mnd, each = mnd), mnd, mnd)
##     matfs1 <- fuTerm2[as.vector(ind1), ] * SuTerm2[as.vector(ind2), ]
##     matfs1 <- t(t(matfs1) * const2)
##     matfs2 <- fuTerm2[as.vector(ind2), ] * SuTerm2[as.vector(ind1), ]
##     matfs2 <- t(t(matfs2) * const2)
##     matfs3 <- fuTerm3[as.vector(ind1), ] * SuTerm3[as.vector(ind2), ]
##     matfs4 <- fuTerm3[as.vector(ind2), ] * SuTerm3[as.vector(ind1), ]
##     sum(matfs1[matfs3 > matfs4]) + 0.5 * sum(matfs1[matfs3 == matfs4])
## }
