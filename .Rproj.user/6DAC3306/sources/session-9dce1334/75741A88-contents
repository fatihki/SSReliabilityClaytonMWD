#' Modified Weibull Distribution
#' 
#' @title The Modified Weibull Distribution (MWD)
#' 
#' @description
#' Density, distribution function, quantile function, random generation,
#' and hazard function for the Modified Weibull distribution (MWD)
#' introduced by Lai et al. (2003).
#' 
#' @name ModifiedWeibull
#' 
#' @import stats
#' 
#' @param x Numeric vector of observations.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities in \eqn{[0,1]}.
#' @param n Integer; number of observations to be generated.
#' 
#' @param a Positive scale parameter (\eqn{a > 0}).
#' @param b Non-negative shape parameter (\eqn{b \ge 0}).
#' @param lambda Non-negative parameter (\eqn{\lambda \ge 0}) controlling
#' the growth rate of the hazard function.
#'
#' @param log Logical. If \code{TRUE}, returns log-density,
#' log-distribution, or log-hazard values where applicable.
#'
#' @param lower.tail Logical. If \code{FALSE}, returns \eqn{1 - F(x)}
#' and computes quantiles for \eqn{1 - p}.
#' 
#' #' The Modified Weibull distribution with parameters \eqn{a},
#' \eqn{b}, and \eqn{\lambda} has cumulative distribution function (CDF),
#' probability density function (PDF), and hazard function given by
#' 
#' @details 
#' The Modified Weibull distribution with parameters \eqn{a}, 
#' \eqn{b} and \eqn{\lambda} has cumulative distribution function (CDF),
#' probability density function (PDF), and hazard function given by
#' 
#' \deqn{F(x) = 1 - \exp\left(-a x^b \exp(\lambda x)\right),}
#' \deqn{f(x) = a (b + \lambda x) x^{b - 1} \exp(\lambda x)
#' \exp\left(-a x^b \exp(\lambda x)\right),}
#' \deqn{h(x) = a (b + \lambda x) x^{b - 1} \exp(\lambda x),}
#' 
#' where \eqn{x > 0}, \eqn{a > 0} is the scale parameter, \eqn{b \ge 0} is a shape parameter,
#' and \eqn{\lambda \ge 0} is an acceleration or flexibility parameter that controls
#' how quickly the hazard grows over time.
#' 
#' Special cases:
#' \itemize{
#'   \item If \eqn{\lambda = 0}, the MWD reduces to the Weibull distribution
#'   \eqn{F(x) = 1 - \exp(-a x^b)}.
#'
#'   \item If \eqn{b = 0}, the MWD reduces to a type I extreme-value
#'   (log-gamma) distribution
#'   \eqn{F(x) = 1 - \exp(-a \exp(\lambda x))}.
#' }
#'
#' @return
#' \itemize{
#'   \item \code{dMweibull}: Density values.
#'   \item \code{pMweibull}: Distribution function values.
#'   \item \code{qMweibull}: Quantiles.
#'   \item \code{rMweibull}: Random deviates.
#'   \item \code{hMweibull}: Hazard function values.
#' }
#' 
#' 
#' @references
#' Lai, C. D., Xie, M., and Murthy, D. N. P. (2003).
#' A modified Weibull distribution.
#' \emph{IEEE Transactions on Reliability}, \strong{52}(1), 33--37.
#'
#' @examples
#' n <- 25
#' a <- 0.75; b <- 1.25; lambda <- 0.60
#' set.seed(123)
#' x <- rlnorm(n)
#' ff <- dMweibull(x, a, b, lambda)
#' FF <- pMweibull(x, a, b, lambda)
#' qq <- qMweibull(runif(n), a, b, lambda)
#' dat <- rMweibull(n, a, b, lambda)
#' hf <- hMweibull(x, a, b, lambda )
#' 
#' @export
dMweibull <- function(x, a, b, lambda, log=FALSE) {
  if (any(x <= 0)) stop("x must be > 0")
  if (a <= 0) stop("Parameter 'a' must be > 0")
  if (b < 0) stop("Parameter 'b' must be >= 0")
  if (lambda < 0) stop("Parameter 'lambda' must be >= 0")
  
  # log-density
  logpdf <- log(a) + (b - 1) * log(x) +
    lambda * x + log(b + lambda * x) -
    a * x^b * exp(lambda * x)

  if (log) return(logpdf)

  pdf <- exp(logpdf)

  return(pdf)
}
#' 
#' @rdname ModifiedWeibull
#' @export
pMweibull <- function(q, a, b, lambda, lower.tail=TRUE, log=FALSE){
  if(any(q < 0)) stop("q must be >=0")
  if (a <= 0) stop("Parameter 'a' must be > 0")
  if (b < 0) stop("Parameter 'b' must be >= 0")
  if (lambda < 0) stop("Parameter 'lambda' must be >= 0")
  
  F <- 1 - exp(-a * q^b * exp(lambda*q))
  if(!lower.tail) F <- 1 - F
  if(log) F <- log(F)
  return(F)
}
#'
#' @rdname ModifiedWeibull
#' @export
qMweibull <- function(p, a, b, lambda, lower.tail=TRUE){
  if(!lower.tail) p <- 1 - p
  if(any(p < 0 | p > 1)) stop("p must be in [0,1]")
  if (a <= 0) stop("Parameter 'a' must be > 0")
  if (b < 0) stop("Parameter 'b' must be >= 0")
  if (lambda < 0) stop("Parameter 'lambda' must be >= 0")
  
  # Handle p = 0 and p = 1 directly
  result <- numeric(length(p))
  result[p == 0] <- 0
  result[p == 1] <- Inf
  
  # Indices for which we need actual computation
  idx <- which(p > 0 & p < 1)
  pp <- p[idx]
  
  if(lambda == 0 && b != 0){ 
    # Two-parameter Weibull quantile
    result[idx] <- (-log(1 - pp)/a)^(1/b)
  } else if(b == 0 && lambda != 0){ 
    # Type I extreme-value (or log-gamma) quantile
    result[idx] <- (1/lambda) * log(-log(1 - pp)/a)
  } else {
    # General case: use uniroot for each p
    result[idx] <- sapply(pp, function(pp_i){
      uniroot(function(x) 1 - exp(-a*x^b*exp(lambda*x)) - pp_i,
              interval = c(0, 1e5), tol = .Machine$double.eps^0.5)$root
    })
  }
  return(result)
}
#'
#' @rdname ModifiedWeibull
#' @export
rMweibull <- function(n, a, b, lambda){
  u <- runif(n)
  qMweibull(u, a, b, lambda)
}
#' @rdname ModifiedWeibull
#' @export
hMweibull <- function(x, a, b, lambda, log=FALSE){
  if (any(x <= 0)) stop("x must be > 0")
  if (a <= 0) stop("Parameter 'a' must be > 0")
  if (b < 0) stop("Parameter 'b' must be >= 0")
  if (lambda < 0) stop("Parameter 'lambda' must be >= 0")

  # log-hazard
  loghz <- log(a) + (b - 1) * log(x) +
    lambda * x + log(b + lambda * x) 
  
  if (log) return(loghz)
  
  hz <- exp(loghz)
  
  return(hz)
}
