#' Plotting an \code{rocForest} object
#'
#' Plots a tree from an \code{rocForest} object. See \code{\link{plot.rocTree}} for more details.
#'
#' @param x an \code{rocForest} object.
#' @param tree an integer specifying the \eqn{n^{th}} tree in the forest to print.
#' @param output a string specifying the output type; graph (the default) renders the graph using the \code{grViz} function, and \code{visNetwork} renders the graph using the visnetwork function.
#' @param digits the number of digits to print.
#' @param rankdir is a character string specifying the direction of the tree flow. The available options are top-to-bottom ("TB"), bottom-to-top ("BT"), left-to-right ("LR"),
#' and right-to-left ("RL"); the default value is "TB".
#' @param shape is a character string specifying the shape style.
#' Some of the available options are "ellipse", "oval", "rectangle", "square", "egg", "plaintext", "diamond", and "triangle". The default value is "ellipse".
#' @param nodeOnly is a logical value indicating whether to display only the node number; the default value is "TRUE".
#' @param savePlot is a logical value indicating whether the plot will be saved (exported); the default value is "FALSE".
#' @param file_name is a character string specifying the name of the plot when "savePlot = TRUE". The file name should include its extension. The default value is "pic.pdf"
#' @param file_type is a character string specifying the type of file to be exported. Options for graph files are: "png", "pdf", "svg", and "ps". The default value is "pdf".
#'
#' @export
plot.rocForest <- function(x, tree = 1L, output = c("graph", "visNetwork"),
                           digits = 4,
                           rankdir = c("TB", "BT", "LR", "RL"),
                           shape = "ellipse",
                           nodeOnly = FALSE,
                           savePlot = FALSE, 
                           file_name = "pic.pdf",
                           file_type = "pdf", ...) {
    if (!is.rocForest(x)) stop("Response must be a \"rocForest\" object.")
    if (!is.wholenumber(tree)) stop("Tree number must be an integer.")
    if (tree > length(x$forest)) stop("Tree number exceeded the number of trees in forest.")
    output <- match.arg(output)
    rankdir <- match.arg(rankdir)
    tmp <- x$forest[[tree]]
    names(tmp)[1] <- "Frame"
    class(tmp) <- "rocTree"
    plot.rocTree(tmp, output = output, digits = digits, rankdir = rankdir,
                 shape = shape, nodeOnly = nodeOnly, savePlot = savePlot,
                 file_name = file_name, file_type = file_type)
}

#' @noRd
#' @export
plot.predict.rocForest <- function(x, ..., control = list()) {
    if (!is.predict.rocForest(x)) stop("Response must be a \"predict.rocForest\" object")
    tmp <- data.frame(x = x$pred[[1]]$Time, y = x$pred[[1]]$Surv)
    gg <- ggplot(tmp, aes(x = x, y = y)) + geom_step(lwd = I(1.1)) +
        xlab("Time") + ylab("Survival probabilities")
    gg
}

