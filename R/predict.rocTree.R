#' Predicting based on a \code{rocTree} model.
#'
#' The function gives predicted values.
#'
#' @param object is an \code{rocTree} object.
#' @param newdata is an optional data frame in which to look for variables with which to predict.
#' If omitted, the fitted predictors are used.
#' If the covariate observation time is not supplied, covariates will be treated as at baseline.
#' @param type is an optional character string specifying whether to predict the survival probability or the hazard rate.
#' @param ... for future developments.
#' 
#' @importFrom stats model.frame
#' @export
predict.rocTree <- function(object, newdata, type = c("survival", "hazard"), ...) {
    if (!is.rocTree(object)) stop("Response must be a \"rocTree\" object")
    type <- match.arg(type)
    ctrl <- object$ctrl
    if (missing(newdata)) {
        xlist <- object$xlist
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
        if (!any(id == names(newdata))) newdata$id <- 1:nrow(newdata)
        else names(newdata)[which(names(newdata) == id)] <- "id"
        p <- length(object$vNames)
        n <- length(unique(newdata$Y))
        newdata$yind <- myFindInt(newdata$Y, object$Y0)
        newdata <- newdata[order(newdata$yind, newdata$Y),]
        X <- cbind(yind = newdata$yind, model.frame(delete.response(object$terms), newdata))
        xlist <- rep(list(matrix(NA, n, length(unique(newdata$id)))), p)
        sptdat <- split(X, newdata$yind)
        X.path <- lapply(1:length(sptdat), function(z) {
            tmp <- sptdat[[z]]
            ind <- tmp[1,1]
            sapply(1:p, function(y)
                object$xlist[[y]][ind, myFindInt.X(tmp[,y+1], object$xlist0[[y]][ind,])])
        })
        X.path <- data.frame(do.call(rbind, X.path))
        for (i in 1:p) {
            xlist[[i]] <- do.call(cbind, lapply(split(X.path[,i], newdata$id), function(z) c(z, rep(NA, n - length(z)))))
        }
    }
    ndInd <- matrix(1, n, length(unique(newdata$id)))
    Frame <- object$Frame
    for (i in 1:nrow(Frame)) {
        if (Frame$terminal[i] == 0) {
            ind <- object$xlist[[Frame$p[i]]] <= Frame$cut[i]
            ndInd[ndInd == Frame$nd[i] & ifelse(is.na(ind), FALSE, ind)] <- Frame$nd[i] * 2
            ndInd[ndInd == Frame$nd[i] & ifelse(is.na(ind), FALSE, !ind)] <- Frame$nd[i] * 2 + 1
        }
    }
    ndInd[is.na(object$xlist[[1]])] <- NA
    dfPred <- matrix(0, n, length(unique(newdata$id)))
    for (i in 1:n) {
        for (k in 1:length(object$ndFinal)) {
            dfPred[i, ndInd[i,] == object$ndFinal[k]] <- object$dfFinal[i, k]
        }
    }
    dfPred[is.na(object$xlist[[1]])] <- NA
    
}

myFindInt <- function(x, y) {
    pmax(1, findInterval(x, sort(y)))
}

myFindInt.X <- function(x, y) {
    order(y)[pmax(1, findInterval(x, sort(y)))]
}

## This is a working version for simulation, still need to think about how to make it more general.
## To-do
## (Y, ID) = ~~(yes, yes)~~, (no, yes), (yes, no), (no, no).
## When Y is not given, assume the same covariates is observed through out?
