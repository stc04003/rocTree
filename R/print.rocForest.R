#' Printing an \code{rocForest} object
#'
#' Prints a tree from an \code{rocForest} object.
#'
#' @param x an \code{rocForest} object.
#' @param tree an integer specifying the \eqn{n^{th}} tree in the forest to print.
#' The function prints the first tree in the forest by default.
#' @param digits the number of digits of numbers to print.
#' @param dt  an optional logical vector. If TRUE, tree structure based on \strong{\code{data.tree} structure} is printed.
#' @param ... for future development.
#'
#' @export
print.rocForest <- function(x, tree = 1L, digits = 5, dt = TRUE, ...) {
    if (!is.rocForest(x)) stop("Response must be a \"rocForest\" object.")
    if (!is.wholenumber(tree)) stop("Tree number must be an integer.")
    if (tree > length(x$forest)) stop("Tree number exceeded the number of trees in forest.")
    Frame <- x$forest[[tree]]$treeMat
    root <- Node$new("Root", type = "root", decision = "", nd = 1)
    for (i in 2:nrow(Frame)) {
        if (i <= 3) parent <- "root"
        if (i > 3) parent <- paste0("Node", Frame$nd[i] %/% 2)
        if (Frame$terminal[i] == 2) {
            type <- "terminal"
            display <- with(Frame, paste0(nd[i], ") ", tree.split.names(nd[i], nd, p, cut, x$vNames, digits), "*"))
        } else {
            type <- "interior"
            display <- with(Frame, paste0(nd[i], ") ", tree.split.names(nd[i], nd, p, cut, x$vNames, digits)))
        }
        eval(parse(text = paste0("Node", Frame$nd[i], "<-", parent,
                                 "$AddChild(display, type = type, nd = Frame$nd[i])")))
    }
    if (!dt) classic.rocTree(x, ...)
    else{
        toPrint <- ToDataFrameTree(root)[[1]]
        cat(" ROC-guided survival tree #", tree, "\n")
        cat("\n")
        cat(" node), split\n")
        cat("   * denotes terminal node\n")
        cat("  ", toPrint, sep = "\n")
        cat("\n")
    }
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
