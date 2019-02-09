#' Predicting based on a \code{rocTree} model.
#'
#' The function gives predicted values.
#'
#' @param object is an \code{rocTree} object.
#' @param newdata is an optional data frame in which to look for variables with which to predict.
#' If omitted, the fitted predictors are used.
#' If the covariate observation time is not supplied, covariates will be treated as at baseline.
#' @param type is an optional character string specifying whether to predict the survival probability or the cumulative hazard rate.
#' @param ... for future developments.
#'
#' @return Returns a \code{data.frame} of the predicted survival probabilities or cumulative hazard. 
#'
#' @seealso \code{\link{predict.rocForest}}
#' @importFrom stats model.frame
#' @export
predict.rocTree <- function(object, newdata, type = c("survival", "hazard"), ...) {
    if (!is.rocTree(object)) stop("Response must be a \"rocTree\" object")
    type <- match.arg(type)
    parm <- object$parm
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
        if (!any(res == names(newdata))) newdata$Y <- max(object$Y0)
        else names(newdata)[which(names(newdata) == res)] <- "Y"
        if (!any(id == names(newdata))) newdata$id <- 1  ##:nrow(newdata)
        else names(newdata)[which(names(newdata) == id)] <- "id"
        p <- length(object$vNames)
        Y <- newdata$Y
        n <- length(unique(Y))
        newdata$yind <- findInt(Y, object$Y0)
        newdata <- newdata[order(newdata$yind, Y),]
        X <- cbind(yind = newdata$yind, model.frame(delete.response(object$terms), newdata))
        xlist <- rep(list(matrix(NA, n, length(unique(newdata$id)))), p)
        sptdat <- split(X, newdata$yind)
        X.path <- lapply(1:length(sptdat), function(z) {
            tmp <- sptdat[[z]]
            ind <- tmp[1,1]
            sapply(1:p, function(y)
                ## object$xlist[[y]][ind, findInt.X(tmp[,y+1], object$xlist0[[y]][ind,])])
                cbind(0, object$xlist[[y]])[ind, findInt.X(tmp[,y+1], object$xlist0[[y]][ind,])]) ## mimicking ecdf
            })
        X.path <- data.frame(do.call(rbind, X.path))
        for (i in 1:p) {
            xlist[[i]] <- do.call(cbind, lapply(split(X.path[,i], newdata$id), function(z) c(z, rep(NA, n - length(z)))))
        }
    }
    ndInd <- matrix(1, n, dim(xlist[[1]])[2])
    dfPred <- ndInd - 1
    Frame <- object$Frame
    for (i in 1:nrow(Frame)) {
        if (Frame$terminal[i] == 0) {
            ind <- xlist[[Frame$p[i]]] <= Frame$cut[i]
            ndInd[ndInd == Frame$nd[i] & ifelse(is.na(ind), FALSE, ind)] <- Frame$nd[i] * 2
            ndInd[ndInd == Frame$nd[i] & ifelse(is.na(ind), FALSE, !ind)] <- Frame$nd[i] * 2 + 1
        }
    }
    ndInd[is.na(xlist[[1]])] <- NA
    for (i in 1:n) {
        for (k in 1:length(object$ndFinal)) {
            dfPred[i, ndInd[i,] == object$ndFinal[k]] <- object$dfFinal[i, k]
        }
    }
    dfPred[is.na(xlist[[1]])] <- NA
    if (type == "survival") {
        pred <- data.frame(Time = sort(unique(Y)),
                           Surv = rowMeans(apply(dfPred, 2, function(x) exp(-cumsum(x))), na.rm = TRUE))
    }
    if (type == "cumHaz") {
        pred <- data.frame(Time = sort(unique(Y)),
                           cumHaz = rowMeans(apply(dfPred, 2, function(x) cumsum(x)), na.rm = TRUE))
    }
    for (i in 1:length(xlist)) attr(xlist[[i]], "dimnames") <- NULL
    out <- list(pred = pred, type = type, xlist = xlist, dfPred = dfPred)
    class(out) <- "predict.rocTree"
    return(out)
}

is.predict.rocTree <- function(x) inherits(x, "predict.rocTree")

#' findInterval with 0 replaced with 1
#' @keywords internal
#' @noRd
findInt <- function(x, y) {
    pmax(1, findInterval(x, sort(y)))
}

#' findInterval with 0 replaced with 1, works with NA's in y
#' @keywords internal
#' @noRd
findInt.X <- function(x, y) {
    order(c(0, y))[pmax(1, findInterval(x, sort(c(0, y))))]
    ## order(y)[pmax(1, findInterval(x, sort(y)))]
}

## This is a working version for simulation, still need to think about how to make it more general.
## To-do
## (Y, ID) = ~~(yes, yes)~~, (no, yes), (yes, no), (no, no).
## When Y is not given, assume the same covariates is observed through out?
