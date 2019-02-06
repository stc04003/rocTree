#' Printing an \code{rocTree} object
#'
#' The function prints an \code{rocTree} object. It is a method for the generic function print of class "\code{rocTree}".
#'
#' @param x an \code{rocTree} object.
#' @param digits the number of digits of numbers to print.
#' @param dt  an optional logical vector. If TRUE, tree structure based on \strong{\code{data.tree} structure} is printed.
#' @param ... for future development.
#'
#' @export
#' @examples
#' set.seed(1)
#' dat <- simu(100, 0, 1.3)
#' library(survival)
#' system.time(fit <- rocTree(Surv(Y, death) ~ z1 + z2, id = id,
#' data = dat, control = list(CV = TRUE, nflds = 10)))
#' fit
#' print(fit, tree = FALSE)
print.rocTree <- function(x, digits = 5, dt = TRUE, ...) {
    if (!is.rocTree(x)) stop("Response must be a \"rocTree\" object.")
    ## digits = getOption("digits")
    Frame <- x$Frame
    ## create data.tree
    root <- Node$new("Root", type = "root", decision = "", nd = 1)
    if (nrow(Frame) > 1) {
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
    }
    if (!dt) classic.rocTree(x, ...)
    else{
        toPrint <- ToDataFrameTree(root)[[1]]
        cat(" ROC-guided survival tree\n")
        if (nrow(Frame) > 1) {        
            cat("\n")
            cat(" node), split\n")
            cat("   * denotes terminal node\n")
            cat("  ", toPrint, sep = "\n")
        } else {
            cat(" Decision tree found no splits.")
        }
        cat("\n")
    }
}

classic.rocTree <- function(x, spaces = 2L, digits = 5, top2bottom = FALSE, ...) {
    if (!is.rocTree(x)) stop("Response must be a \"rocTree\" object")
    ## digits = getOption("digits")
    ## Prepare Frame and remove nodes after considering ndFinal
    Frame <- data.frame(x$treeMat)
    Frame <- subset(Frame, !is.na(Frame$u))
    Frame$terminal[Frame$nd == 1] <- 0 ## root can't be a terminal node (keep?)
    if (!is.null(x$ndFinal)) {
        Frame$terminal[which(Frame$nd %in% x$ndFinal)] <- 2
        if (sum((!(Frame$nd %in% x$ndFinal) & Frame$terminal == 2)) > 0) 
            Frame <- Frame[-which(!(Frame$nd %in% x$ndFinal) & Frame$terminal == 2),]
        tmp <- Frame$nd[Frame$terminal == 2]
        tmp2 <- tree.depth(tmp)
        rm1 <- c(sapply(x$nd[which(tmp2 < max(tmp2))],
                        function(a) unlist(sapply(1:(1 + max(tmp2)), function(b) (2^b * a):(2^b * (a + 1) - 1)))))
        rm <- unique(c(rm1, (max(tmp) + 1):(max(Frame$nd) + 1)))
        Frame <- Frame[Frame$nd %in% setdiff(Frame$nd, rm),]
    }   
    node <- Frame$nd
    depth <- tree.depth(node)
    if (length(node) == 1) {indent <- paste0(format(node), ")")
    } else {
        indent <- paste(rep(" ", spaces * 32L), collapse = "")
        indent <- substring(indent, 1L, spaces * seq(depth))
        indent <- paste0(c("", indent[depth]), format(node), ")")
    }
    termNd <- rep(" ", length(depth))
    termNd[Frame$terminal == 2] <- "*"
    slab <- with(Frame, sapply(nd, function(z) tree.split.names(z, nd, p, cut, x$vNames, digits)))
    toPrint <- paste(indent, slab, termNd)
    if (!top2bottom) {
        k1 <- max(depth)
        k2 <- 2^ceiling(log(max(node), base = 2) + 1e-7) - 1
        tmpM <- matrix(rep(1:k2, 2^(k1 - tree.depth(1:k2))), k1 + 1, byrow = TRUE)
        ind <- sapply(node, function(z) which.min(z - t(tmpM) > 0))
        rk <- ind %% 2^k1
        rk <- ifelse(rk > 0, rk, 2^k1)
        rk <- rank(rk, ties.method = "first")
        toPrint <- toPrint[order(rk)]
    }
    cat(" ROC-guided survival tree\n")
    cat("\n")
    cat(" node), split\n")
    cat("   * denotes terminal node\n")
    cat(" ", toPrint, sep = "\n")
    cat("\n")
}

tree.depth <- function(nodes) {
    depth <- floor(log(nodes, base = 2) + 1e-7)
    depth - min(depth)
}

tree.split.names <- function(nd0, nd, p, cut, xname, digits = getOption("digits")) {
    if (nd0 == 1) return("root")
    ind <- which(nd == nd0 %/% 2)
    if (nd0 %% 2 == 0) {
        return(paste(xname[p[ind]], "<=", formatC(cut[ind], digits = digits, flag = "#")))
    } else {
        return(paste(xname[p[ind]], ">", formatC(cut[ind], digits = digits, flag = "#")))
    }
}

#' @export
print.predict.rocTree <- function(x, tree = 1L, ...) {
    if (!is.wholenumber(tree)) stop("Tree number must be an integer.")
    if (names(x$pred)[[2]] == "Surv") {
        cat(" Fitted survival probabilities:\n")
    }
    if (names(x$pred)[[2]] == "cumHaz") {
        cat(" Fitted cumulative hazard:\n")
    }
    print(head(x$pred, 5))
    cat("\n")
}
