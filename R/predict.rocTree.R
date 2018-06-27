#' Predicting based on a \code{rocTree} model.
#'
#' The function gives predicted values.
#'
#' @param x an \code{rocTree} object.
#' @param newdata an optional data frame in which to look for variables with which to predict.
#' If omitted, the fitted predictors are used.
#' If the covariate observation time is not supplied, covariates will be treated as at baseline.
#' 
#' @export
predict.rocTree <- function(x, newdata, type = c("survival", "hazard"), ...) {
    if (!is.rocTree(x)) stop("Response must be a \"rocTree\" object")
    type <- match.arg(type)
    ctrl <- x$ctrl
    if (missing(newdata)) {
        xlist <- foo$xlist
    } else {
        res <- x$terms[[2]][[2]]
        id <- attr(x$terms, "id")
        if (all(is.na(newdata[,names(newdata) == res]))) res <- NULL
        if (all(is.na(newdata[,names(newdata) == id]))) id <- NULL
        newdata <- model.frame(paste(res, "~", paste(c(x$vNames, id), collapse = "+")), newdata)
        res <- x$terms[[2]][[2]]
        id <- attr(x$terms, "id")
        if (!any(res == names(newdata))) newdata$Y <- min(x$Y0)
        else names(newdata)[which(names(newdata) == res)] <- "Y"
        if (!any(id == names(newdata))) newdata$id <- 1:nrow(newdata)
        else names(newdata)[which(names(newdata) == id)] <- "id"
        p <- length(x$vNames)
        newdata <- newdata[order(newdata$id, newdata$Y),]
        newdata$yind <- findInterval(newdata$Y, x$Y0, left.open = TRUE, all.inside = TRUE)
        xlist <- rep(list(matrix(NA, length(unique(newdata$Y)), length(unique(newdata$id)))), p)
        for (i in unique(newdata$id)) {
            newdatai <- subset(newdata, id == i)
            Xi <- model.frame(delete.response(x$terms), newdatai)
            for (j in 1:p) {
                tmp <- sapply(1:length(yind), function(z)
                    x$xlist[[j]][yind[z], findInterval(Xi[z, p], sort(x$xlist0[[j]][yind[z],]),
                                                       left.open = TRUE, all.inside = TRUE)])
                xlist[[j]][1:length(tmp), i] <- tmp
            }
        }
    }
    return(xlist)
}

## To-do
## (Y, ID) = ~~(yes, yes)~~, (no, yes), (yes, no), (no, no).
## When Y is not given, assume the same covariates is observed through out?
