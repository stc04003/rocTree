#' Class definition
setClass("splitClass",
         representation(tau = "numeric", M = "numeric", hN = "numeric", h = "numeric",
                        minsp = "numeric", minsp2 = "numeric", disc = "numeric", nflds = "numeric",
                        CV = "logical", Trace = "logical", parallel = "logical", parCluster = "numeric",
                        B = "numeric", ghN = "numeric", fsz = "numeric"),
         prototype(tau = 0.4, M = 1000, hN = 0, h = 0, minsp = 20, minsp2 = 5, disc = -1, nflds = 10,
                   CV = FALSE, Trace = FALSE, parallel = FALSE,
                   parCluster = parallel::detectCores() / 2,
                   B = 500, ghN = 0.2, fsz = 0),
         contains = "VIRTUAL")
setClass("CON", contains = "splitClass")
setClass("dCON", contains = "splitClass")
#' Method dispatch
setGeneric("grow", function(Y, E, X.list, parm) standardGeneric("grow"))

setMethod("grow", signature(parm = "CON"), growSeq)
setMethod("grow", signature(parm = "dCON"), growRP)

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
#' \item{r2Final}{estimated hazards at all terminal nodes.}
#' \item{tt}{time values to estimate hazard. Default length is 500. }
#' }
#' @references Sun Y. and Wang, M.C. (2018+). ROC-guided classification and survival trees. \emph{Technical report}.
#' @keywords rocTree
#' @seealso See \code{\link{print.rocTree}} and \code{\link{plot.rocTree}} for printing and plotting an \code{rocTree}, respectively.
#' @examples
#' library(survival)
#' set.seed(123)
#' dat <- simu(100, 0, 1.1)
#' fit <- rocTree(Surv(Y, death) ~ z1 + z2, id = id, data = dat,
#'        control = list(CV = TRUE, nflds = 5))
#' fit
rocTree <- function(formula, data, id, subset, splitBy = c("CON", "dCON"), control = list()) {
    splitBy  <- match.arg(splitBy)
    parm.control <- control[names(control) %in% names(attr(getClass(splitBy), "slots"))]
    parm <- do.call("new", c(list(Class = splitBy), parm.control))
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
    if (parm@disc[1] < 0) parm@disc <- rep(0, p)
    if (length(parm@disc) == 1) parm@disc <- rep(parm@disc, p)
    xlist <- sapply(1:p, function(z) rocTree.Xlist(X[,z], parm@disc[z], Y, id), simplify = FALSE)
    xlist0 <- sapply(1:p, function(z) rocTree.Xlist(X[,z], 1, Y, id), simplify = FALSE) ## for prediction
    Y0 <- unlist(lapply(split(Y, id), max), use.names = FALSE)
    E0 <- unlist(lapply(split(Status, id), max), use.names = FALSE)
    if (parm@h <= 0) parm@h <- parm@tau / 20
    if (parm@hN <= 0) parm@hN <- parm@tau / 20    
    if (parm@fsz <= 0) parm@fsz <- round(length(Y0) / 2)
    ## Grow
    out <- grow(Y0, E0, xlist, parm)
    ## CV
    if (parm@CV & length(out$beta.seq) > 2) 
        out$con2.seq <- rocTree.cv(out$beta.seq, Y, Status, id, X, ctrl)
    out <- c(out, rocTree.final(out, Y0, E0, xlist, ctrl))
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
    out$ctrl <- ctrl
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

#' Cross-validation, evaluated when CV = TRUE
#'
#' The resulting tree should be smaller than the full grew tree in rocTree.
#' This is the same procedure with splitBy = "dCON" or "CON".
#'
#' @keywords internal
#' @noRd
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
    EE1 <- Y1 <= tau * E1
    EE2 <- Y2 <= tau * E2
    const <- 1 / rowMeans(fmat)[EE1] / sc[EE1] / N1
    sc2 <- fitc$surv[findInterval(Y2, fitc$time)]
    Smat2 <- matk2 <- outer(ifelse(E2, Y2, NA), Y2, "<=") / matrix(sc2, N2, N2, TRUE)
    fmat2 <- matk <- t(E2 * sapply(ifelse(E2, Y2, NA), K2, vec = Y2, h = hN2) / sc2 / hN2)
    if (min(rowMeans(fmat2)[E2 == 1]) < 0) 
        fmat[E2 == 1 & rowMeans(fmat2) < 0, ] <- fmat2[which(rowMeans(fmat2) > 0)[1], ]
    const2 <- 1 / rowMeans(fmat2)[EE2] / sc2[EE2] / N2
    Smat3 <- outer(ifelse(E2, Y2, NA), Y1, "<=") 
    fmat3 <- t(E1 * sapply(ifelse(E2, Y2, NA), K2, vec = Y1, h = hN1) / hN1)
    fTree <- STree <- matrix(NA, M, sum(EE1))
    fTree2 <- STree2 <- fTree3 <- STree3 <- matrix(NA, M, sum(EE2))
    fTree[1, ] <- rowMeans(fmat)[EE1]
    STree[1, ] <- rowMeans(Smat)[EE1]
    fTree2[1, ] <- rowMeans(fmat2)[EE2]
    STree2[1, ] <- rowMeans(Smat2)[EE2]
    fTree3[1, ] <- rowMeans(fmat3)[EE2]
    STree3[1, ] <- rowMeans(Smat3)[EE2]
    treeMat <- data.frame(nd = 1:M, terminal = 0, u = NA, u2 = NA, p = NA, cut = NA)
    treeMat[1, ] <- c(1, 1, 1, 1, NA, NA)
    ndInd1 <- matrix(1, N1, N1)
    ndInd2 <- matrix(1, N2, N2)
    ndInd12 <- matrix(1, N2, N1)
    ndInd1[lower.tri(ndInd1)] <- 0
    ndInd2[lower.tri(ndInd2)] <- 0
    ndInd12[is.na(X12.list[[1]])] <- 0
    tlst <- which(Y1 >= tau)[1] - 1
    while (sum(treeMat$terminal == 1) > 0) {
        sp <- splitTree(X1.list, Y1, E1, fmat, Smat, treeMat, ndInd1, const, fTree, STree, control)
        if (sp[1] * 2 < M & !is.na(sp[2])) {
            ndInd1[ndInd1 == sp[1] & X1.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd1[ndInd1 == sp[1] & X1.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            ndInd12[ndInd12 == sp[1] & X12.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd12[ndInd12 == sp[1] & X12.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            fTree[sp[1] * 2, ] <- rowSums(fmat[EE1, diag(ndInd1) == 2 * sp[1], drop = FALSE]) / N1
            fTree[sp[1] * 2 + 1, ] <- rowSums(fmat[EE1, diag(ndInd1) == 2 * sp[1] + 1, drop = FALSE]) / N1
            STree[sp[1] * 2, ] <- rowSums(as.matrix(Smat * (ndInd1 == 2 * sp[1])))[EE1] / N1
            STree[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat * (ndInd1 == 2 * sp[1] + 1)))[EE1] / N1
            fTree2[sp[1] * 2, ] <- rowSums(fmat2[EE2, diag(ndInd2) == 2 * sp[1], drop = FALSE]) / N2
            fTree2[sp[1] * 2 + 1, ] <- rowSums(fmat2[EE2, diag(ndInd2) == 2 * sp[1] + 1, drop = FALSE]) / N2
            STree2[sp[1] * 2, ] <- rowSums(as.matrix(Smat2 * (ndInd2 == 2 * sp[1])))[EE2] / N2
            STree2[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat2 * (ndInd2 == 2 * sp[1] + 1)))[EE2] / N2
            ## rTree2[sp[1]*2,] <- fTree2[sp[1]*2,]/STree2[sp[1]*2,]
            ## rTree2[sp[1]*2+1,] <- fTree2[sp[1]*2+1,]/STree2[sp[1]*2+1,]
            fTree3[sp[1] * 2, ] <- rowSums(fmat3[EE2, diag(ndInd1) == 2 * sp[1], drop = FALSE]) / N1
            fTree3[sp[1] * 2 + 1, ] <- rowSums(fmat3[EE2, diag(ndInd1) == 2 * sp[1] + 1, drop = FALSE]) / N1
            STree3[sp[1] * 2, ] <- rowSums(as.matrix(Smat3 * (ndInd12 == 2 * sp[1])))[EE2] / N1
            STree3[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat3 * (ndInd12 == 2 * sp[1] + 1)))[EE2] / N1
            ## rTree3[sp[1]*2,] <- fTree3[sp[1]*2,]/STree3[sp[1]*2,]
            ## rTree3[sp[1]*2+1,] <- fTree3[sp[1]*2+1,]/STree3[sp[1]*2+1,]
            treeMat[sp[1], 2] <- 0
            treeMat[sp[1], 5:6] <- sp[2:3]
            treeMat[sp[1] * 2, ] <- c(sp[1] * 2, 1,
                                      mean(diag(ndInd1) == sp[1] * 2 & EE1),
                                      min(rowMeans(ndInd1[Y1 <= tau, ] == sp[1] * 2)), NA, NA)
            treeMat[sp[1] * 2 + 1, ] <- c(sp[1] * 2 + 1, 1,
                                          mean(diag(ndInd1) == (sp[1] * 2 + 1) & EE1),
                                          min(rowMeans(ndInd1[Y1 <= tau, ] == sp[1] * 2 + 1)), NA, NA)
            treeMat$terminal[which(treeMat$u <= minsp / N1 & treeMat$u2 <= minsp2 / N1)] <- 2
            ## conTree <- conTree + sp[4]
        } else {
            if (sum(!is.na(treeMat$p)) > 1) break
            else treeMat$terminal[treeMat$nd == sp[1]] <- 2
            break
        }
    }
    if (all(is.na(treeMat$p))) {
        cat("\nNo splits.\n")
        return(NULL)
    }
    ## prune
    treeMatTerm <- treeMat$nd[treeMat$terminal >= 1 & is.na(treeMat$p)]
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
                       as.integer(length(fuTerm) / sum(EE1)),
                       as.integer(sum(EE1)),
                       as.double(t(fuTerm)),
                       as.double(t(SuTerm)),
                       as.double(const), 
                       double(1))[[6]]
    conList2[[1]] <- .C("con3", 
                        as.integer(length(fuTerm2) / sum(EE2)),
                        as.integer(sum(EE2)),
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
                fuTerm <- fTree[treeMat$nd %in% ndNew, ]
                SuTerm <- STree[treeMat$nd %in% ndNew, ]
                fuTerm2 <- fTree2[treeMat$nd %in% ndNew, ]
                SuTerm2 <- STree2[treeMat$nd %in% ndNew, ]
                fuTerm3 <- fTree3[treeMat$nd %in% ndNew, ]
                SuTerm3 <- STree3[treeMat$nd %in% ndNew, ]
                consub <- c(consub, .C("con", 
                                       as.integer(length(fuTerm) / sum(EE1)),
                                       as.integer(sum(EE1)),
                                       as.double(t(fuTerm)),
                                       as.double(t(SuTerm)),
                                       as.double(const), 
                                       double(1))[[6]])
                con2sub <- c(con2sub, .C("con3", 
                                         as.integer(length(fuTerm2) / sum(EE2)),
                                         as.integer(sum(EE2)),
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
    cvconbb <- function(bb) res[max(1, which(res[, 2] >= bb)[1] - 1), 4]
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
