#' Plotting an \code{rocForest} object
#'
#' Plots a tree from an \code{rocForest} object.
#'
#' @param x  an \code{rocForest} object.
#' @param tree an integer specifying the \eqn{n^{th}} tree in the forest to print.
#' @param output a string specifying the output type; graph (the default) renders the graph using the \code{grViz} function, and \code{visNetwork} renders the graph using the visnetwork function.
#' @param digits the number of digits to print.
#' @param control a list of control parameters. See "details" for important features of control parameters.
#' @param ... arguments to be passed to or from other methods.
#' 
#' @export
plot.rocForest <- function(x, tree = 1L, output = c("graph", "visNetwork"),
                           digits = 4, control = list(), ...) {
    if (!is.rocForest(x)) stop("Response must be a \"rocForest\" object.")
    if (!is.wholenumber(tree)) stop("Tree number must be an integer.")
    if (tree > length(x$forest)) stop("Tree number exceeded the number of trees in forest.")
    tmp <- x$forest[[tree]]
    names(tmp)[1] <- "Frame"
    class(tmp) <- "rocTree"
    plot.rocTree(tmp, output = output, digits = digits, control = control)
}
