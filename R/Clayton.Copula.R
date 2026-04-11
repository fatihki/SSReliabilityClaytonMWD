#' 
#' @title Two-dimensional Clayton Copula
#' 
#' @description
#' This is the joint cumulative distribution and probability density functions
#' for the two-dimensional Clayton copula.
#' 
#' @name Clayton_Copula 
#' 
#' @param u is a vector of points in \eqn{[0,1]} representing the first coordinate of the copula.
#' @param v is a vector of points in \eqn{[0,1]} representing the first coordinate of the copula.
#' @param theta defines \eqn{\theta} dependency parameter between the pairs.
#' 
#' @details  
#' The joint distribution function of the two-dimensional Clayton copula, along with its joint
#' probability density function, are given by
#' \deqn{C(u,v;\theta) = \left(u^{-\theta} + v^{-\theta} - 1\right)^{-1/\theta}, \; \theta \in (0, \infty).}
#' \deqn{c_{\theta}(u,v) = (\theta +1)  u^{-(\theta + 1)} v^{-(\theta + 1)} \left( u^{-\theta} + v^{-\theta}-1 \right)^{-\left (\frac{1}{\theta} + 2 \right)}}
#' where \eqn{\theta > 0}.
#' 
#' @return \code{Clayton_Copula} gives the joint cumulative distribution function, 
#' \code{Clayton_Copula_pdf} gives the joint probability density function values for the Clayton copula evaluated at points \eqn{(u,v)} 
#' with given \eqn{\theta} parameter.
#' 
#' @export
Clayton_Copula <- function(u, v, theta) {
  if(theta == 0) return(u * v) # Independent case
  return( pmax( u^(-theta) + v^(-theta) - 1, 0)^(-1/theta) )
}
#' 
#' @rdname Clayton_Copula
#' @export
Clayton_Copula_pdf <- function(u, v, theta) {

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
