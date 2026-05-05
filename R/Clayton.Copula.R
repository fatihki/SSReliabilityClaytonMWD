#' 
#' @title Two-dimensional Clayton Copula
#' 
#' @description
#' Computes the joint cumulative distribution function (CDF) and probability
#' density function (PDF) of the two-dimensional Clayton copula.
#' 
#' @name Clayton_Copula 
#' 
#' @param u Numeric vector of values in \eqn{[0,1]}. First marginal (uniform).
#' @param v Numeric vector of values in \eqn{[0,1]}. Second marginal (uniform).
#' @param theta Positive numeric scalar. Dependence parameter 
#'   \eqn{\theta > 0}.
#'
#' @details
#' The joint distribution function of the two-dimensional Clayton copula is
#' \deqn{
#' C(u,v;\theta) = \left(u^{-\theta} + v^{-\theta} - 1\right)^{-1/\theta},
#' }
#' where \eqn{\theta > 0}.
#'
#' The corresponding joint density is given by
#' \deqn{
#' c(u,v;\theta) = (\theta + 1) u^{-(\theta + 1)} v^{-(\theta + 1)}
#' \left(u^{-\theta} + v^{-\theta} - 1\right)^{-\left(1/\theta + 2\right)}.
#' }
#' 
#' 
#' @return
#' \itemize{
#'   \item \code{Clayton_Copula}: Numeric vector of CDF values.
#'   \item \code{Clayton_Copula_pdf}: Numeric vector of PDF values.
#' }
#'
#' @examples
#' u <- c(0.2, 0.5, 0.8)
#' v <- c(0.3, 0.6, 0.9)
#'
#' Clayton_Copula(u, v, theta = 2)
#' Clayton_Copula_pdf(u, v, theta = 2)
#'
#' @references
#' Nelsen, R. B. (2006). \emph{An Introduction to Copulas}. Springer.
#' 
#' @export
Clayton_Copula <- function(u, v, theta) {
  
  if (any(u < 0 | u > 1)) stop("u must be in [0,1]")
  if (any(v < 0 | v > 1)) stop("v must be in [0,1]")
  if (theta <= 0) stop("theta must be > 0")
  
  if(theta == 0) return(u * v) # Independent case
  return( pmax( u^(-theta) + v^(-theta) - 1, 0)^(-1/theta) )
}
#' 
#' @rdname Clayton_Copula
#' @export
Clayton_Copula_pdf <- function(u, v, theta) {
  
  if (any(u < 0 | u > 1)) stop("u must be in [0,1]")
  if (any(v < 0 | v > 1)) stop("v must be in [0,1]")
  if (theta <= 0) stop("theta must be > 0")
  
  if(theta <= 0) return(rep(0, length(u)))
  
  eps <- 1e-10
  u <- pmin(pmax(u, eps), 1 - eps)
  v <- pmin(pmax(v, eps), 1 - eps)
  
  p <- (theta + 1) * (u * v)^(-(theta + 1)) *
    (u^(-theta) + v^(-theta) - 1)^(-2 - 1/theta)
  return(p)
}
#'
#' The first derivative of Clayton copula function wrt theta, which is used in the one-step LSE and WLSE.
#' @noRd
dClayton_Copula_theta <- function(u, v, theta) {
  S <- u^(-theta) + v^(-theta) - 1
  term1 <- log(S)/theta^2
  term2 <- (u^(-theta)*log(u) + v^(-theta)*log(v)) / (theta*S)
  return( Clayton_Copula(u,v,theta) * (term1 + term2) )
}
#'
#' The second derivative of Clayton copula function wrt theta, which is used in the one-step LSE and WLSE.
#' @noRd
d2Clayton_Copula_theta <- function(u, v, theta) {
  S <- u^(-theta) + v^(-theta) - 1
  C <- Clayton_Copula(u,v,theta)
  dC <- dClayton_Copula_theta(u,v,theta)
  dS <- -u^(-theta)*log(u) - v^(-theta)*log(v)
  d2S <- u^(-theta)*(log(u))^2 + v^(-theta)*(log(v))^2
  
  d1 <- log(S)/(theta^2) + dS/(theta*S)
  # derivative of d1
  d1_prime <- - 2*log(S)/theta^3 - (dS^2)/(theta*S^2) + d2S/(theta*S)
  # the second derivative of C
  d2 <- C * (d1^2 + d1_prime)
  return(d2)
}
