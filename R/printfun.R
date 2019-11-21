#' Extract CBQ Coefficients
#'
#' Create a table of coefficient results from a \code{cbq} object.
#'
#' @param object A \code{cbq} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A table of coefficients with their corresponding lower and upper bounds.
#' @export
#'
#' @method coef cbq
#'
coef.cbq <- coefficients.cbq <- function(object, ...) {
  coefmat <-
    cbind(matrix(object$means, nrow = object$npars), t(object$ulbs))
  row.names(coefmat) <- c(object$xnames)
  colnames(coefmat) <- c("Estimate", "LB", "UB")
  return(coefmat)
}


#' Print cbq object
#'
#' General print function for \code{cbq} objects, which dispatches the chosen type
#' of printing to the corresponding function.
#'
#' @param x A \code{cbq} object to be printed.
#' @param type Character string giving the type of printing, such as
#'   \code{"text"}, \code{"mcmc"}, \code{"coef"}.
#' @param ... Additional arguments to be passed to print functions.
#'
#' @export
#'
#'
#'
#'
print.cbq <- function(x, type = "text", ...) {
  printFunName <- paste0("print_", type, ".cbq")
  do.call(printFunName, args = c(list(object = x), list(...)))
}


#' Print the main results from a \code{cbq} object.
#'
#' @param object A \code{cbq} object.
#' @param digits Number of digits to display.
#'
#' @export
#'
#'
#'
print_text.cbq <- function(object, digits = 3) {
  cat("Conditional binary quantile regression \n")
  cat("\nCall:\n",
      paste(deparse(object$Call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = "")
  cat(
    "MCMC run for",
    object$nsim,
    "iterations, with",
    object$stanfit@sim$warmup2,
    "used. \n\n"
  )
  cat("Coefficients:\n")
  print(round(coef(object), digits))
  cat("\n")
}


#' Print the mcmc results from a cbq object
#'
#' This prints a number of diagnostics about the results of a \code{cbq} objects
#'
#'
#' @param object A \code{cbq} object.
#' @param ... Additional arguments to be passed.
#'
#' @export
#'
#'
#'
print_mcmc.cbq <- function(object, ...) {
  print(object$stanfit, ...)
}



#' Print cbq coefficients
#'
#' @param object A \code{cbq} object.
#' @param digits Number of digits to display.
#'
#' @export
#'
#'
#'
print_coef.cbq <- function(object, digits = 3) {
  print(round(coef(object), digits))
}
