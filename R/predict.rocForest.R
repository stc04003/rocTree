#' Predicting based on a \code{rocForest} model.
#'
#' The function gives predicted values.
#'
#' @param object is an \code{rocForest} object.
#' @param newdata is an optional data frame in which to look for variables with which to predict.
#' If omitted, the fitted predictors are used.
#' If the covariate observation time is not supplied, covariates will be treated as at baseline.
#' @param type is an optional character string specifying whether to predict the survival probability or the hazard rate.
#' @param ... for future developments.
#'
#' @seealso \code{\link{predict.rocTree}}
#' @export
predict.rocForest <- function(object, newdata, type = c("survival", "hazard", "cumHaz"), ...) {
    if (!is.rocForest(object)) stop("Response must be a \"rocForest\" object")
    type <- match.arg(type)
    ctrl <- object$ctrl
    if (missing(newdata)) {
        xlist <- object$xlist
        Y <- object$Y0
        n <- length(Y)
    } else {
        res <- object$terms[[2]][[2]]
        id <- attr(object$terms, "id")
        if (all(is.na(newdata[,names(newdata) == res]))) res <- NULL
        if (all(is.na(newdata[,names(newdata) == id]))) id <- NULL
        newdata <- model.frame(paste(res, "~", paste(c(object$vNames, id), collapse = "+")), newdata)
        res <- object$terms[[2]][[2]]
        id <- attr(object$terms, "id")
        if (!any(res == names(newdata))) newdata$Y <- min(object$Y0)
        else names(newdata)[which(names(newdata) == res)] <- "Y"
        if (!any(id == names(newdata))) newdata$id <- 1 #:nrow(newdata)
        else names(newdata)[which(names(newdata) == id)] <- "id"
        p <- length(object$vNames)
        Y <- newdata$Y
        n <- length(unique(Y))
        newdata$yind <- findInt(Y, object$Y0)
        newdata <- newdata[order(newdata$yind, Y),] ## edit here
        X <- cbind(yind = newdata$yind, model.frame(delete.response(object$terms), newdata))
        xlist <- rep(list(matrix(NA, n, length(unique(newdata$id)))), p)
        sptdat <- split(X, newdata$yind)
        X.path <- lapply(1:length(sptdat), function(z) {
            tmp <- sptdat[[z]]
            ind <- tmp[1,1]
            sapply(1:p, function(y)
                cbind(0, object$xlist[[y]])[ind, findInt.X(tmp[,y+1], object$xlist0[[y]][ind,])]) ## mimicking ecdf
        })
        X.path <- data.frame(do.call(rbind, X.path))
        for (i in 1:p) {
            xlist[[i]] <- do.call(cbind, lapply(split(X.path[,i], newdata$id), function(z) c(z, rep(NA, n - length(z)))))
        }
    }
    nID <- dim(xlist[[1]])[2]
    ndInd <- matrix(1, n, nID)
    dfPred <- ndInd - 1
    Y0 <- sort(unique(Y))
    t0 <- seq(0, quantile(Y0, 0.8), length = 50)
    if (type %in% c("survival", "cumHaz")) {
        W <- oneW(ndInd, xlist, object$forest[[1]])
        for (i in 2:length(object$forest)) {
            tmp <- oneW(ndInd, xlist, object$forest[[i]])
            W <- lapply(1:nID, function(x) tmp[[x]] + W[[x]])
            rm(tmp)
        }
        W <- lapply(W, matrix, n)
        ## matk <- sapply(t0, function(z) object$E0 * K3(z, Y0, object$ctrl$ghN) / object$ctrl$ghN)
        ## matk2 <- outer(t0, Y0, "<=")
        matk3 <- outer(Y0, Y0, "==") * object$E0
        pred <- list()
        W0 <- upper.tri(matrix(0, n, n), TRUE) * 1:n / n ## for NA's in Wi
        for (i in 1:nID) {
            Vi <- colSums(W[[i]]) / length(object$forest)
            Wi <- W[[i]] / rowSums(W[[i]])
            Wi[rowSums(W[[i]]) == 0,] <- W0[rowSums(W[[i]]) == 0,]
            if (type == "survival") {
                pred[[i]] <- data.frame(Time = Y0, Surv = exp(-cumsum(rowSums(matk3 * Wi))))
            }
            if (type == "cumHaz") {
                pred[[i]] <- data.frame(Time = Y0, cumHaz = cumsum(rowSums(matk3 * Wi)))
            }
        }
    }
    if (type == "hazard") {
        W <- oneV(ndInd, xlist, object$forest[[1]])
        for (i in 2:length(object$forest)) {
            tmp <- oneV(ndInd, xlist, object$forest[[i]])
            W <- lapply(1:nID, function(x) tmp[[x]] + W[[x]])
            rm(tmp)
        }
        ## matk <- sapply(Y0, function(z) object$E0 * K3(z, Y0, object$ctrl$ghN) / object$ctrl$ghN)
        matk <- sapply(t0, function(z) object$E0 * K3(z, Y0, object$ctrl$ghN) / object$ctrl$ghN)
        pred <- list()
        for (i in 1:nID) {
            ## pred[[i]] <- data.frame(Time = Y0, haz = colSums(matk * W[[i]]) / length(object$forest))
            pred[[i]] <- data.frame(Time = t0, haz = colSums(matk * W[[i]][,findInterval(t0, Y0) + 1]) / length(object$forest))
        }
    }
    for (i in 1:length(xlist)) attr(xlist[[i]], "dimnames") <- NULL
    out <- list(xlist = xlist, W = W, pred = pred)
    class(out) <- "predict.rocForest"
    return(out)
}  

#' This function provides one weight matrix using a tree from the forest
#'
#' ##@importFrom memoise memoise
#' @keywords internal
#' @noRd
oneV <- function(ndInd, xlist, tree) {
    szL2 <- tree$szL2
    ndInd2 <- tree$ndInd2
    idB2 <- tree$idB2
    Frame <- tree$treeMat
    ndTerm <- Frame$nd[Frame$terminal >= 1]
    if (nrow(Frame) == 1) ndTerm <- 1
    for (i in 1:dim(Frame)[1]) {
        if (Frame$terminal[i] == 0) {
            ndInd[ndInd == Frame$nd[i] & xlist[[Frame$p[i]]] <= Frame$cut[i]] <- Frame$nd[i] * 2
            ndInd[ndInd == Frame$nd[i] & xlist[[Frame$p[i]]] > Frame$cut[i]] <- Frame$nd[i] * 2 + 1
        }
    }
    lapply(1:dim(xlist[[1]])[2], function(z) giveV(ndInd[,z], idB2, ndInd2, ndTerm, szL2))
}

#' @keywords internal
#' @noRd
oneW <- function(ndInd, xlist, tree) {
    szL2 <- tree$szL2
    ndInd2 <- tree$ndInd2
    idB2 <- tree$idB2
    Frame <- tree$treeMat
    ndTerm <- Frame$nd[Frame$terminal >= 1]
    if (nrow(Frame) == 1) ndTerm <- 1
    for (i in 1:dim(Frame)[1]) {
        if (Frame$terminal[i] == 0) {
            ndInd[ndInd == Frame$nd[i] & xlist[[Frame$p[i]]] <= Frame$cut[i]] <- Frame$nd[i] * 2
            ndInd[ndInd == Frame$nd[i] & xlist[[Frame$p[i]]] > Frame$cut[i]] <- Frame$nd[i] * 2 + 1
        }
    }
    ## n <- dim(xlist[[1]])[1]
    lapply(1:dim(xlist[[1]])[2], function(z) giveW(ndInd[,z], idB2, ndInd2, ndTerm, szL2))
}

#' This function provides W in \code{oneWeight} to avoid \code{for} loops.
#'
#' @keywords internal
#' @noRd
giveW <- function(ndi, idB2, ndInd2, ndTerm, szL2) {
    n <- length(ndi)
    .C("giveWC", as.integer(n), as.integer(length(idB2)), as.integer(length(ndTerm)),
       as.integer(ndi), as.integer(idB2 - 1), as.integer(ndInd2), as.integer(ndTerm),
       as.double(szL2), out = double(n * n), PACKAGE = "rocTree")$out
}

#' @importFrom dplyr dense_rank
giveV <- function(ndi, idB2, ndInd2, ndTerm, szL2) {
    n <- length(ndi)
    tmp <- sapply(sort(unique(ndi)), function(z) 
        .C("giveVC", as.integer(n), as.integer(length(idB2)), as.integer(length(ndTerm)),
           as.integer(rep(z, n)), as.integer(idB2 - 1), as.integer(ndInd2), as.integer(ndTerm),
           as.double(szL2), out = double(n), PACKAGE = "rocTree")$out)
    return(tmp[, dense_rank(ndi)])
}
