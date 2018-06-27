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
        X <- model.frame(delete.response(x$terms), newdata)
        xlist <- list(matrix(NA, length(unique(newdata$Y)), length(unique(newdata$id))),
                      matrix(NA, length(unique(newdata$Y)), length(unique(newdata$id))))
        for (i in unique(newdata$id)) {
            
        }
        
    }
}

findInterval(100, 1:10, left.open = TRUE, rightmost.closed = TRUE)

findInterval(100, 1:10, left.open = TRUE, rightmost.closed = TRUE)


## To-do
## (Y, ID) = ~~(yes, yes)~~, (no, yes), (yes, no), (no, no).
## When Y is not given, assume the same covariates is observed through out?
