#' Printing an \code{rocForest} object
#'
#' Prints a tree from an \code{rocForest} object.
#'
#' @param x an \code{rocForest} object.
#' @param tree an optional integer specifying the \eqn{n^{th}} tree in the forest to print.
#' The function prints the contents of an \code{rocForest} object by default. 
#' @param digits the number of digits of numbers to print.
#' @param dt  an optional logical vector. If TRUE, tree structure based on \strong{\code{data.tree} structure} is printed.
#' @param ... for future development.
#'
#' @seealso \code{\link{rocForest}}, \code{\link{print.rocTree}}
#'
#' @export
print.rocForest <- function(x, tree = NULL, digits = 5, dt = TRUE, ...) {
    if (!is.rocForest(x)) stop("Response must be a \"rocForest\" object.")
    if (!is.null(tree)) {
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
        if (!dt) classic.rocTree(Frame, ...)
        else{
            toPrint <- ToDataFrameTree(root)[[1]]
            cat(" ROC-guided survival forest: tree #", tree, "\n")
            cat("\n")
            cat(" node), split\n")
            cat("   * denotes terminal node\n")
            cat("  ", toPrint, sep = "\n")
            cat("\n")
        }
    }
    if (is.null(tree)) {
        cat("rocForest result\n\n")
        cat("Call:\n", deparse(x$call), "\n\n")
        cat("Sample size:                                 ", ncol(x$xlist[[1]]), "\n")
        cat("Number of independent variables:             ", length(x$xlist), "\n")
        cat("Number of trees:                             ", x$parm@B, "\n")
        cat("Split rule:                                  ", x$parm@splitBy, "\n")
        cat("Number of variables tried at each split:     ", ceiling(sqrt(length(x$xlist))), "\n")
        cat("Size of subsample:                           ", x$parm@fsz, "\n")
        cat("Minimum number of failure in a node:         ", x$parm@minsp, "\n")
        cat("Minimum number of failure in a teminal node: ", x$parm@minsp2, "\n")
    }
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#' @export
print.predict.rocForest <- function(x, tree = 1L, ...) {
    if (!is.wholenumber(tree)) stop("Tree number must be an integer.")
    if (tree > length(x$pred)) stop("Tree number exceeded the number of subjects in newdata.")
    if (names(x$pred[[tree]])[[2]] == "Surv") {
        cat(" Fitted survival probabilities for subject #", tree, ":\n")
    }
    if (names(x$pred[[tree]])[[2]] == "cumHaz") {
        cat(" Fitted cumulative hazard for subject #", tree, ":\n")
    }
    print(head(x$pred[[tree]], 5))
    cat("\n")
}
