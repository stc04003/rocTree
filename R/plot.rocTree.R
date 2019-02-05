#' Plotting an \code{rocTree} object
#'
#' Plots an \code{rocTree} object. The function uses a \code{dgr_graph} object, render the graph in the RStudio Viewer.
#'
#' The argument \code{control} consists of some options to be passed down to \code{render_graph} of the \code{DiagrammeR} package.
#' The argument defaults to a list with the following values:
#' \describe{
#'   \item{rankdir}{a character vector specifying the drawing direction.
#' The default value is \code{TB} (top-to-bottom).
#' Other options include \code{BT} (bottom-to-top), \code{LR} (left-to-right), and \code{RL} (right-to-left).}
#'   \item{shape}{a character vector specifying the shape of nodes. The default value is "ellipse". See the polygon-based nodes in \code{graphviz} for more options.}
#'   \item{nodeOnly}{a logical vector specifying whether to display only the node number.}
#' }
#'   
#' @param x an object of class "\code{rocTree}", usually returned by the rocTree function.
#' @param output a string specifying the output type; graph (the default) renders the graph using the \code{grViz} function, and \code{visNetwork} renders the graph using the visnetwork function.
#' @param digits the number of digits to print.
#' @param control a list of control parameters. See "details" for important features of control parameters.
#' @param ... arguments to be passed to or from other methods.
#'
#' @seealso See \code{\link{rocTree}} for creating \code{rocTree} objects.
#' @export
#' @examples
#' data(simudat)
#' library(survival)
#' system.time(fit <- rocTree(Surv(Time, Status) ~ X1 + X2 + X3, id = ID,
#' data = simudat, control = list(CV = TRUE, nflds = 10)))
#' plot(fit)
#' plot(fit, control = list(rankdir = "LR", nodeOnly = TRUE))
#' plot(fit, output = "visNetwork", control = list(nodeOnly = TRUE, shape = "box"))
plot.rocTree <- function(x, output = c("graph", "visNetwork"),
                         digits = 4, control = list(), ...) {
    if (!is.rocTree(x)) stop("Response must be a \"rocTree\" object")
    output <- match.arg(output)
    ctrl <- plot.rocTree.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
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
    SetGraphStyle(root, rankdir = ctrl$rankdir)
    if (ctrl$nodeOnly) {
        GetNodeLabel <- function(node) {
            switch(node$type,
                   root = "Root", 
                   terminal = paste0("Node ", node$nd),
                   interior = paste0("Node ", node$nd))
        }
        SetNodeStyle(root, fontname = 'helvetica', label = GetNodeLabel, shape = ctrl$shape)
    } else {
        SetNodeStyle(root, fontname = 'helvetica', shape = ctrl$shape)
    }
    ## x$graph <- ToDiagrammeRGraph(x, direction = "climb", pruneFun = NULL)
    plot.Node(root, output = output, control = ctrl)
}

plot.rocTree.control <- function(rankdir = c("TB", "BT", "LR", "RL", "TD"),
                                 shape = "ellipse",
                                 nodeOnly = FALSE,
                                 savePlot = FALSE, 
                                 file_name = "pic.pdf",
                                 file_type = "pdf") {
    list(rankdir = match.arg(rankdir), shape = shape, nodeOnly = nodeOnly, savePlot = savePlot,
         file_name = file_name, file_type = file_type, hN = 0.2)
}

#' Plotting the estimated hazard function from an rocTree object
#'
#' Plot the estimated hazard function from rocTree objects.
#'
#' The argument "control" defaults to a list with the following values:
#' \describe{
#'   \item{\code{ghN}}{smoothing parameter used in smoothing hazard functions; the default value is 0.2.}
#' }
#'
#' @param x an object of class "rocTree", usually returned by the rocTree function.
#' @param control A list of control parameters. See "details" for important special features of control parameters.
#' @export
#'
#' @importFrom ggplot2 ggplot geom_line xlab ylab aes
#' 
#' @examples
#' data(simudat)
#' library(survival)
#' system.time(fit <- rocTree(Surv(Time, Status) ~ X1 + X2 + X3, id = ID,
#' data = simudat, control = list(CV = TRUE, nflds = 10)))
#' plotTreeHaz(fit)
plotTreeHaz <- function(x, control = list()) {
    if (is.null(control$tau)) control$tau <- x$parm@tau
    if (is.null(control$ghN)) control$ghN <- x$parm@ghN    
    if (!is.rocTree(x)) stop("Response must be a \"rocTree\" object")
    x <- c(x, rocTree.haz(x$dfFinal, x$Y0, control))
    tmp <- data.frame(x = x$tt, y = unlist(x$r2Final), Node = rep(x$ndFinal, each = length(x$tt)), row.names = NULL)
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
