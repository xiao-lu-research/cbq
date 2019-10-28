#' Plot cbq object
#'
#' General plot function for \code{cbq} objects, which dispatches the chosen
#' type of plotting to the corresponding function.
#'
#' @param object A \code{cbq} object to be plotted.
#' @param type Character string giving the type of plotting. The options are
#'   \code{"trace"} for trace plots, \code{"coef"} for coefficient plots.
#' @param ... Additional arguments to be passed to subsequent plot functions.
#'
#' @export
#'
#'
#' @examples
#'
#'
plot.cbq <- function(object, type = "trace", ...) {
  printFunName <- paste0("plot_", type, ".cbq")
  do.call(printFunName, args = c(list(object = object), list(...)))
}



#' Make traceplots for cbq
#'
#' Plot traceplots from a \code{cbq} object. 
#'
#' @param object A \code{cbq} object.
#' @param ... Additional parameters to be passed.
#'
#' @export
#'
#' @examples
#'
#'
plot_trace.cbq <- function(object, ...) {
  rstan::traceplot(object$stanfit,...)
}

#' Make coefficient plots for cbq
#'
#' Plot traceplots from a \code{cbq} object. 
#'
#' @param object A \code{cbq} object.
#' @param ... Additional parameters to be passed.
#'
#' @export
#'
#' @examples
#'
#'
plot_coef.cbq <- function(object, ...) {
  rstan::plot(object$stanfit,...)
}
