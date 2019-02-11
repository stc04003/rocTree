#' rocTree:Receiver Operating Characteristic (ROC)-Guided Classification Survival Tree and Forest.
#'
#' The \code{rocTree} package uses a Receiver Operating Characteristic guided classification
#' algorithm to grow and prune survival trees.
#' The \code{rocTree} package also provides implementation to grow random forest.
#' 
#'
#' @aliases rocTree-package
#' @section Introduction:
#' The \code{rocTree} package provides implementations to a unified framework for
#' tree-structured analysis
#' with censored survival outcomes.
#' Different from many existing tree building algorithms,
#' the \code{rocTree} package incorporates time-dependent covariates by constructing
#' a time-invariant partition scheme on the survivor population.
#' The partition-based risk prediction function is constructed using an algorithm guided by
#' the Receiver Operating Characteristic (ROC) curve.
#' Specifically, the generalized time-dependent ROC curves for survival trees show that the
#' target hazard function yields the highest ROC curve.
#' The optimality of the target hazard function motivates us to use a weighted average of the
#' time-dependent area under the curve on a set of time points to evaluate the prediction
#' performance of survival trees and to guide splitting and pruning.
#' Moreover, the \code{rocTree} package also offers a novel risk prediction forest algorithm,
#' where the ensemble is on unbiased martingale estimating equations.
#' 
#' @section Methods:
#' The package contains functions to construct ROC-guided survival trees (\code{\link{rocTree}})
#' and random forest (\code{\link{rocForest}}).
#' 
#' @seealso \code{\link{rocTree}}, \code{\link{rocForest}}
#' @seealso For more details, see the vignettes at \link{https://www.sychiou.com/rocTree/index.html}.
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
