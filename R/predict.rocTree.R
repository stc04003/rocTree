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
        if (!(x$terms[[2]][[2]] %in% names(newdata))) newdata$Y <- min(x$Y0)
        if (!(attr(x$terms, "id") %in% names(newdata))) newdata$id <- 1:nrow(newdata)
        else names(newdata)[which(names(newdata) == attr(x$terms, "id"))] <- "id"
    }
}

