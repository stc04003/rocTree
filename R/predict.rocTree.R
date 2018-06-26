#' Predicting based on a \code{rocTree} model.
#'
#' The function gives predicted values.
#'
#' @param x an \code{rocTree} object.
#' @param newdata an optional data frame in which to look for variables with which to predict.
#' If omitted, the fitted predictors are used.
#' @export
predict.rocTree <- function(x, newdata, type = c("survival", "hazard"), ...) {
    if (!is.rocTree(x)) stop("Response must be a \"rocTree\" object")
    type <- match.arg(type)
    ctrl <- x$ctrl
    if (missing(newdata)) {
        xlist <- foo$xlist
    } else {
        Terms <- delete.response(x$terms)
        newdata <- model.frame(Terms, newdata, na.action = na.pass)
        check <- attr(Terms, "dataClasses")
        if (!is.null(check)) .checkMFClasses(check, newdata, TRUE)
        p <- ncol(newdata)
        xlist <- sapply(1:p, function(z) rocTree.Xlist(X[,z], ctrl$disc[z], Y, id), simplify = FALSE)
    }

    ## rocTree.final(x$treeMat, )
}

