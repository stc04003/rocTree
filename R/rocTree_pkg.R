#' rocTree:Receiver Operating Characteristic (ROC)-Guided Classification and Survival Tree
#' 
#' \code{rocTree} uses a ROC-guided classification algorithm to grow and prune survival trees.
#' Hazards at terminal notes are also computed.
#'
#' @aliases rocTree-package
#' @section Introduction:
#' Existing classification, regression, and survival trees are typically implemented via a
#' greedy algorithm that maximizes the within node homogeneity or between node heterogeneity.
#' Sun and Wang (2018+) re-define the ROC curves for tree classifiers using the idea of
#' randomized test and show that the target function yields the highest ROC curve.
#' Then a criterion based on ROC is used in the tree-building algorithm. 
#' 
#' @section Methods:
#' The package contains functions to construct ROC-guided survival trees (\code{\link{rocTree}}).
#'
#' @examples
#' library(survival)
#' data(simudat)
#' system.time(fit <- rocTree(Surv(Time, Status) ~ X1 + X2, id = ID, data = simudat,
#' control = list(CV = TRUE, nflds = 10)))
#' fit
#' 
#' @seealso \code{\link{rocTree}}
#' @seealso For more details, see the \code{rocTree} vignette by running: \code{vignette("rocTree")}
#' @docType package
#' 
#' @importFrom stats model.extract model.matrix model.response
#' @importFrom utils tail
#' @importFrom survival survfit Surv
#' @importFrom parallel detectCores makeCluster setDefaultCluster clusterExport stopCluster
#' @importFrom parallel parSapply parLapply
#' @importFrom graphics legend lines plot
#' @importFrom data.tree Node ToDataFrameTree ToDiagrammeRGraph SetGraphStyle SetNodeStyle
#' @importFrom DiagrammeR render_graph %>% export_graph
#'
#' @useDynLib rocTree, .registration = TRUE
"_PACKAGE"
NULL

## Use @name rocTree for package's name if there is no function that share the same name
## #' @name rocTree
