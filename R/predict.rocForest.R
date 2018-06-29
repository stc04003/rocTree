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
predict.rocForest <- function(object, newdata, type = c("survival", "hazard"), ...) {
    if (!is.rocTree(object)) stop("Response must be a \"rocForest\" object")
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
        if (!any(id == names(newdata))) newdata$id <- 1:nrow(newdata)
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
}  
