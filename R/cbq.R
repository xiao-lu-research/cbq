#' Check if a predictor is dichotomous, adopted from package \code{circGLM}
#'
#' @param x A character or numerical vector to be tested.
#'
#' @return A logical, \code{TRUE} if the \code{x} has dummy coding (0, 1),
#'   \code{FALSE} otherwise.
#'
is.dichotomous <- function(x) {
  n_unique_x <- length(unique(x))
  if (n_unique_x == 2) {
    if (all(x == 0 | x == 1)) {
      return(TRUE)
    } else {
      warning("A predictor might be dichotomous but not 0|1.")
    }
  } else if (n_unique_x > 2 & n_unique_x < 8) {
    warning(
      paste(
        "A predictor has between 3 and 7 unique values.",
        "It might be categorical with multiple categories",
        "but without dummy coding."
      )
    )
  } else if (n_unique_x == 1) {
    stop("A predictor had only a single unique value.")
  }
  FALSE
}

#' Fitting conditional binary quantile models
#'
#' The main function for running the conditional binary quantile model. The function returns a cbq \code{cbq} object that can be further investigated using standard functions such as \code{plot}, \code{print}, \code{coef}, and \code{predict}.
#'
#' The model can be passed either as a combination of a \code{formula} and a data frame \code{data}, as in \code{lm()}.
#'
#' Convergence diagnotics can be performed using either \code{print(object, "mcmc")} or \code{plot(object, "mcmc")}.
#'
#' @param formula An object of class "Formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param q The quantile value.
#' @param nsim The number of iterations.
#' @param burnin The number of burnin iterations.
#' @param thin Thinning parameter.
#' @param CIsize The size of confidence interval.
#' @param nchain The number of parallel chains.
#' @param seeds Random seeds to replicate the results.
#' @param inverse_distr If FALSE, the ALD will not be reversed. The default is FALSE.
#' @param offset Offset values to enhance sampling stability. The default value is 1e-20.
#'
#' @return A \code{cbq} object, which can be further analyzed with its associated \code{\link{plot.cbq}}, \code{\link{coef.cbq}} and \code{\link{print.cbq}} functions.
#'
#' An object of class \code{cbq} contains the following elements
#'
#'   \describe{
#'
#'   \item{\code{Call}}{The matched call.}
#'   \item{\code{formula}}{Symbolic representation of the model.}
#'   \item{\code{q}}{The quantile value.}
#'   \item{\code{nsim}}{The number of MCMC iterations.}
#'   \item{\code{burnin}}{The number of burnin periods.}
#'   \item{\code{thin}}{Thinning.}
#'   \item{\code{seeds}}{Random seeds.}
#'   \item{\code{CIsize}}{The size of confidence interval.}
#'   \item{\code{data}}{Data used.}
#'   \item{\code{x}}{Covaraites used.}
#'   \item{\code{y}}{The dependent variable.}
#'   \item{\code{xnames}}{Names of the covariates.}
#'   \item{\code{stanfit}}{Outputs from stan.}
#'   \item{\code{sampledf}}{A matrix of posterior samples.}
#'   \item{\code{summaryout}}{A summary based on posterior samples.}
#'   \item{\code{npars}}{Number of covariates.}
#'   \item{\code{ulbs}}{Lower and upper confidence bounds.}
#'   \item{\code{means}}{Estimates at the mean.}
#'
#'
#' }
#'
#' @importFrom stats aggregate coef model.frame model.matrix quantile
#'
#' @export
#'
#' @author
#' Xiao Lu
#'
#' @references
#' Lu, Xiao (forthcoming). Discrete Choice Data with Unobserved Heterogeneity: A Conditional Binary Quantile Model. Political Analysis. https://doi.org/10.1017/pan.2019.29
#'
#'
#' @examples
#' # Simulate the data
#' x <- rnorm(50)
#' y <- ifelse(x > 0, 1, 0)
#' dat <- as.data.frame(cbind(y, x))
#'
#' # Estimate the CBQ model
#' model <- cbq(y ~ x, dat, 0.5, inverse_distr = FALSE, nsim = 1000)
#'
#' # Show the results
#' print(model)
#' coef(model)
#' plot(model)
#'
#'
cbq <- function(formula,
                data,
                q = NULL,
                # quantile
                nsim = 1000,
                burnin = NULL,
                thin = 1,
                CIsize = .95,
                nchain = 1,
                seeds = 12345,
                inverse_distr = FALSE,
                offset = 1e-20) {
  if (is.null(burnin))
    burnin <- floor(nsim / 2)
  if (burnin < 0)
    stop("Burn-in must be non-negative.")
  if (thin < 1)
    stop("Thinning factor must be positive.")
  if (CIsize <= 0)
    stop("Confidence interval size 'CIsize' must be positive.")
  if (CIsize > 1)
    stop(paste0("Confidence interval size 'CIsize' ",
                "can not be larger than 1."))
  if (missing(formula) | missing(data)) {
    stop(paste0("Formula and data should be given."))
  }
  if (is.null(q)) {
    warning("Quantile is not specified. The default quantile 0.5 is used.")
    q <- 0.5
  }
  if (length(q) > 1)
    stop("Only a single quantile value is allowed.")
  if (q >= 1 |
      q <= 0)
    stop("The specified quantile is out of range. The value must be in (0,1).")

  f <- Formula::Formula(formula)
  data <- stats::model.frame(f, data)
  y <- c(as.matrix(stats::model.frame(f, data)[, 1]))

  x = stats::model.matrix(f, data)
  if (grepl("\\|", deparse(f))) {
    xq <- as.matrix(stats::model.matrix(f, data, rhs = 2)[, -1])
    nq <- dim(xq)[2]
  } else {
    xq <- NULL
    nq <- 0
  }

  n_covariate <- dim(x)[2]
  N <- length(y)

  if (nq > 1)
    stop("The number of index variables must be equal to or less than one.")
  if (!is.dichotomous(y))
    stop("The dependent variable must be binary with values in {0,1}.")

  if (nq == 1) {
    indx = as.integer(as.factor(c(xq)))
    if (!all(stats::aggregate(y, list(indx), sum)[, 2] == 1)) {
      stop(
        "In each choice set, there must be only one chosen observation. Multiple 1s or all 0s are not allowed in any choice set."
      )
    }

  }

  if (nq == 0) {
    if (inverse_distr == FALSE) {
      stanmodel <- stanmodels$cbqbv
      datlist <- list(
        N = N,
        Y = c(y),
        D = n_covariate,
        X = x,
        p = q,
        offset = offset
      )

    } else {
      stanmodel <- stanmodels$cbqb
      datlist <- list(
        N = N,
        Y = c(y),
        D = n_covariate,
        X = x,
        p = q,
        offset = offset
      )

    }
  } else {
    if (inverse_distr == FALSE) {
      x <- x[order(indx, y), ]
      y <- y[order(indx, y)]
      stanmodel <- stanmodels$cbqdv
      datlist <- list(
        N = N,
        Y = c(y),
        D_common = n_covariate,
        X_common = x,
        N_indx = length(unique(indx)),
        ind = indx,
        q = q,
        offset = offset
      )

    } else {
      x <- x[order(indx, y), ]
      y <- y[order(indx, y)]
      stanmodel <- stanmodels$cbqd
      datlist <- list(
        N = N,
        Y = c(y),
        D_common = n_covariate,
        X_common = x,
        N_indx = length(unique(indx)),
        ind = indx,
        q = q,
        offset = offset
      )

    }
  }

  pars <- "beta"
  stanout <- sampling(
    stanmodel,
    data = datlist,
    pars = pars,
    seed = seeds,
    iter = nsim,
    thin = thin,
    warmup = burnin,
    chains = nchain
  )

  summaryout <- rstan::summary(stanout)$summary
  sampledf <- as.data.frame(stanout)[, 1:n_covariate]

  out = list()
  class(out) <- c("cbq", class(out))
  out$Call <- match.call()
  out$formula <- formula
  out$q <- q
  out$nsim <- nsim
  out$burnin <- burnin
  out$thin <- thin
  out$seeds <- seeds
  out$CIsize  <- CIsize
  out$data   <- data
  out$x    <- x
  out$y <- y
  out$xnames <- colnames(x)
  out$stanfit <- stanout
  out$sampledf <- sampledf
  out$summaryout <- summaryout
  out$npars <- n_covariate
  out$ulbs <-
    apply(sampledf, 2, stats::quantile, probs = c((1 - CIsize) / 2, 1 - (1 -
                                                                           CIsize) / 2))
  out$means <- summaryout[1:n_covariate, 1]

  return(out)


}
