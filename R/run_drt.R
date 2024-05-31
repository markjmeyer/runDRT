#' Run Doubly Ranked Tests
#'
#' Performs two (or more) sample doubly ranked tests on pre-processed functional data,
#' formatted as either a matrix (for functions) or an array (for surfaces).
#'
#' @param X an n by T matrix or an S by T by n array containing the functions (or surfaces) to analyze.
#' @param G a vector of length n containing the grouping variable.
#' @param method statistic for summarizing the ranks: 'suff.rank' for sufficient statistic (the default) or
#' 'avg.rank' for arithmetic average.
#' @param data.names a vector of length two containing names that describe `X` and `G`.
#' @param formula a formula of the form `X ~ G`.
#' @param ... additional arguments to pass to [run_drt.default()], e.g. `method`.
#' @details
#' Doubly ranked tests are non-parametric tests that first rank functions (or surfaces) by time (or location). Next, the procedure
#' summarizes the observed ranks using a statistic. The summarized ranks are then analyzed using either a Wilcoxon rank sum
#' test or a Kruskal-Wallis test. To perform a doubly ranked test, realizations of functions must be stored in an n by T matrix where
#' n is the total number of observed functions and T is the number of realizations per function (commonly time points or locations). Surface
#' data in an S by T by n array can be analyzed as well, although currently this feature has under gone only limited testing.
#'
#' By default, `run_drt()` implements a sufficient statistic when summarizing the ranks of each observed function across T, i.e.
#' the argument `method` defaults to `method = suff.rank`. This statistic has the form
#' \deqn{t(z) = \frac{1}{T}\sum_{t=1}^T\log\left[ \left(\frac{z_t}{n}- \frac{1}{2n}\right)\bigg/\left(1-\frac{z_t}{n} + \frac{1}{2n}\right) \right],}
#' where \eqn{z_t} is the observed rank at time \eqn{t}. See Meyer (2024) for additional details. The average rank may also be
#' used by setting `method = 'avg.rank'`, although this summary has not undergone testing in the doubly ranked context.
#'
#' Regardless of the statistic used, the summarized ranks are the analyzed using either [wilcox.test()] or [kruskal.test()],
#' depending on the number of groups in `G`.
#'
#' For functional data, Meyer (2024) suggests using [refund::fpca.face()] for pre-processing the data, but `X` can be pre-processed using any functional
#' data approach or it can just be the raw data. `run_drt()` itself performs no pre-processing and takes `X` as inputted.
#'
#' @return A list with class "`htest`" containing the following components:\tabular{ll}{
#'    \code{statistic} \tab the value of the test statistic with a name describing it. \cr
#'    \tab \cr
#'    \code{parameter} \tab  the parameter(s) for the exact distribution of the test statistic. \cr
#'    \tab \cr
#'    \code{p.value} \tab the p-value for the test. \cr
#'    \tab \cr
#'    \code{null.value} \tab the location parameter. \cr
#'    \tab \cr
#'    \code{alternative} \tab a character string describing the alternative hypothesis. \cr
#'    \tab \cr
#'    \code{data.name} \tab a character string giving the names of the data. \cr
#'    \tab \cr
#'    \code{test_details} \tab the output from the internally run Wilcoxon rank sum or Kruskal-Wallis test. \cr
#'    \tab \cr
#'    \code{method} \tab character string giving the type of doubly ranked test performed. \cr
#'    \tab \cr
#'    \code{ranks} \tab a list containing the ranks by column (if `X` is a matrix) and the summarized ranks. \cr
#'    \tab \cr
#'    \code{data} \tab a list containing `X` and `G`. \cr
#' }
#'
#' @references Meyer, MJ (2024). Doubly ranked tests for grouped functional data. \emph{Available on arXiv} at \url{https://arxiv.org/abs/2306.14761}.
#'
#' @examples
#' #### Two Sample Problem: Resin Viscosity ####
#' library(FDboost)
#' data("viscosity")
#'
#' Xv    <- matrix(viscosity$visAll, nrow = nrow(viscosity$visAll), ncol = ncol(viscosity$visAll))
#' fXv   <- refund::fpca.face(Xv)
#' Yvis  <- fXv$Yhat
#' TR  <- viscosity$T_A
#'
#' run_drt(Yvis ~ TR)
#'
#' #### Four Sample Problem: Canadian Weather ####
#' R     <- fda::CanadianWeather$region
#' XT    <- t(fda::CanadianWeather$dailyAv[,,'Temperature.C'])
#' fXT   <- refund::fpca.face(XT)
#' YT    <- fXT$Yhat
#'
#' run_drt(YT ~ R)
#'
#' @export
run_drt <- function(X, G,
                    method = c('suff.rank','avg.rank'),
                    data.names = NULL){
  UseMethod("run_drt")
}
#' @rdname run_drt
#' @export
run_drt.default <- function(X, G,
                            method = c('suff.rank','avg.rank'),
                            data.names = NULL){

  if(length(method) > 1){
    method  <- 'suff.rank'
  }

  if(method != 'suff.rank' & method != 'avg.rank'){
    stop("method must be either suff.rank or avg.rank")
  }

  if(length(dim(X)) == 2){
    array_flag    <- FALSE
    N             <- nrow(X)
  } else if(length(dim(X)) == 3){
    array_flag    <- TRUE
    N             <- dim(X)[3]
  } else{
    stop('X must be either an N x S matrix or a T x S x N array')
  }

  if(!is.null(dim(G)) | length(G) != N){
    stop('G must be a vector of length N')
  }

  if(array_flag){
    ##### construct ranks for surfaces #####
    R       <- rank_array(X, na.last = 'keep')
    Y       <- switch(method,
                      suff.rank = apply(R, 3, calc_suff_stat, n = dim(R)[3]),
                      avg.rank = apply(R, 3, mean))
  } else{
    ##### construct ranks for curves #####
    R       <- apply(X, 2, rank, na.last = 'keep')
    Y       <- switch(method,
                      suff.rank = apply(R, 1, calc_suff_stat, n = nrow(R), na.rm = TRUE),
                      avg.rank = rowMeans(R, na.rm = TRUE))
  }

  number_bins_in_G <- length(unique(G))

  if(number_bins_in_G == 2){
    test_details       <- suppressWarnings(stats::wilcox.test(Y ~ G))
    stat_method <- switch(method,
                          suff.rank = 'Doubly Ranked Mann-Whitney-Wilcoxon Test with Sufficient Statistic',
                          avg.rank = 'Doubly Ranked Mann-Whitney-Wilcoxon Test with Arithmetic Mean'
    )
  } else{
    test_details   <- suppressWarnings(stats::kruskal.test(Y ~ G))
    stat_method <- switch(method,
                          suff.rank = 'Doubly Ranked Kruskal-Wallis Test with Sufficient Statistic',
                          avg.rank = 'Doubly Ranked Kruskal-Wallis Test with Arithmetic Mean'
    )

  }

  data_list  <- list(X = X, G = G)
  if(is.null(data.names)){
    data.names <- names(data_list)
  }
  out       <- list(statistic = test_details$statistic,
                    parameter = NULL,
                    p.value = test_details$p.value,
                    null.value = test_details$null.value,
                    alternative = test_details$alternative,
                    data.name = paste(data.names, collapse = ' by '),
                    test_details = test_details,
                    method = stat_method,
                    ranks = list(Y = Y, R = R),
                    data = data_list)

  class(out) <- 'htest'

  return(out)

}
#' @rdname run_drt
#' @export
run_drt.formula <- function(formula, ...){
  X           <- eval(formula[[2]])
  G           <- eval(formula[[3]])
  data.names  <- c(formula[[2]], formula[[3]])

  out   <- run_drt(X, G, data.names = data.names, ...)
  return(out)
}
