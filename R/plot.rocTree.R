#' Plotting an \code{rocTree} object
#'
#' Plots an \code{rocTree} object. The function returns a \code{dgr_graph} object and is rendered in the RStudio Viewer.
#'   
#' @param x an object of class "\code{rocTree}", usually returned by the rocTree function.
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
#' @param ... arguments to be passed to or from other methods.
#'
#' @seealso See \code{\link{rocTree}} for creating \code{rocTree} objects.
#' @export
#' @examples
#' set.seed(1)
#' dat <- simu(100, 0, 1.3)
#' library(survival)
#' system.time(fit <- rocTree(Surv(Time, death) ~ z1 + z2, id = id,
#' data = dat, control = list(prune = TRUE, nflds = 10)))
#' plot(fit)
#' plot(fit, rankdir = "LR", nodeOnly = TRUE)
#' plot(fit, output = "visNetwork", nodeOnly = TRUE, shape = "box")
plot.rocTree <- function(x, output = c("graph", "visNetwork"),
                         digits = 4,
                         rankdir = c("TB", "BT", "LR", "RL"),
                         shape = "ellipse",
                         nodeOnly = FALSE,
                         savePlot = FALSE, 
                         file_name = "pic.pdf",
                         file_type = "pdf", ...) {
    if (!is.rocTree(x)) stop("Response must be a \"rocTree\" object")
    output <- match.arg(output)
    rankdir <- match.arg(rankdir)
    Frame <- x$Frame
    ## create data.tree
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
    SetGraphStyle(root, rankdir = rankdir)
    if (nodeOnly) {
        GetNodeLabel <- function(node) {
            switch(node$type,
                   root = "Root", 
                   terminal = paste0("Node ", node$nd),
                   interior = paste0("Node ", node$nd))
        }
        SetNodeStyle(root, fontname = 'helvetica', label = GetNodeLabel, shape = shape)
    } else {
        SetNodeStyle(root, fontname = 'helvetica', shape = shape)
    }
    ## x$graph <- ToDiagrammeRGraph(x, direction = "climb", pruneFun = NULL)
    plot.Node(root, output = output, control = list(savePlot = savePlot, file_name = file_name, file_type = file_type))
}

#' Plotting the estimated hazard function from an rocTree object
#'
#' Plot the estimated hazard function from rocTree objects.
#'
#' @param x an object of class "rocTree", usually returned by the rocTree function.
#' @param ghN an optional smoothing parameter used in smoothing hazard functions;
#' the default value is 0.2.
#' @export
#'
#' @importFrom ggplot2 ggplot geom_line xlab ylab aes
#' 
#' @examples
#' set.seed(1)
#' dat <- simu(100, 0, 1.3)
#' library(survival)
#' system.time(fit <- rocTree(Surv(Time, death) ~ z1 + z2, id = id,
#' data = dat, control = list(prune = TRUE, nflds = 10)))
#' plotTreeHaz(fit)
plotTreeHaz <- function(x, ghN = NULL) {
    if (!is.rocTree(x)) stop("Response must be a \"rocTree\" object")
    if (!is.null(ghN)) x$parm@ghN <- ghN
    x <- c(x, rocTree.haz(x$dfFinal, x$Y0, x$parm))
    tmp <- data.frame(x = x$tt, y = unlist(x$r2Final),
                      Node = rep(x$ndFinal, each = length(x$tt)), row.names = NULL)
    tmp$Node <- factor(tmp$Node)
    gg <- ggplot(tmp, aes(x = x, y = y, group = Node)) +
        geom_line(aes(linetype = Node, color = Node), lwd = I(1.1)) +
        xlab("Time") + ylab("Hazard")
    gg
}

## ---------------------------------------------------------------------------------------
## From data.tree package that are not exported by them
## ---------------------------------------------------------------------------------------

plot.Node <- function(x, ..., output = "graph", control) {
    graph <- ToDiagrammeRGraph(x, direction = "climb", pruneFun = NULL)
    if (control$savePlot) {
        export_graph(graph, file_name = control$file_name, file_type = control$file_type)
        render_graph(graph, output = output)
    } else {
        render_graph(graph, output = output)
    }
}

#' Plotting the predicted survival function from an rocTree object
#'
#' Plot the predicted survival function from rocTree objects.
#'
#' @importFrom ggplot2 ggplot geom_step
#'
#' @noRd
#' @export
plot.predict.rocTree <- function(x, ..., control = list()) {
    if (!is.null(control$tau)) x$parm@tau <- control$tau
    if (!is.null(control$ghN)) x$parm@ghN <- control$ghN 
    if (!is.predict.rocTree(x)) stop("Response must be a \"predict.rocTree\" object")
    tmp <- data.frame(x = x$pred$Time, y = x$pred$Surv)
    gg <- ggplot(tmp, aes(x = x, y = y)) + geom_step(lwd = I(1.1)) +
        xlab("Time") + ylab("Survival probabilities")
    gg
}
