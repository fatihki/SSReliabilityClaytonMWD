#'
#' @title Bivariate Random Data Generation for MWD Marginals via Clayton Copula
#'
#' @description
#' Generates bivariate random samples from a dependent stress–strength model where
#' both marginals follow the Modified Weibull Distribution (MWD), and the dependence
#' structure between the variables is modeled using a Clayton copula.
#' 
#' @name rMweibull_Clayton 
#'
#' @param n Integer; number of observations to be generated.
#' @param a1,b1,lambda1 MWD parameters for the strength variable \eqn{X},
#' with \eqn{a_1 > 0}, \eqn{b_1 \ge 0}, and \eqn{\lambda_1 \ge 0}.
#' @param a2,b2,lambda2 MWD parameters for the stress variable \eqn{Y},
#' with \eqn{a_2 > 0}, \eqn{b_2 \ge 0}, and \eqn{\lambda_2 \ge 0}.
#' @param theta Clayton copula parameter with \eqn{\theta >0}.
#' 
#' @details
#' Further details are provided in Kızılaslan (2026).
#' 
#' @return A list containing:
#' \item{(U, V)}{Independent generated uniform numbers for \eqn{(X,Y)}.}
#' \item{(X, Y)}{Dependent generated \eqn{n} pairs \eqn{(X,Y)} observations.}
#' 
#' @references
#' Kızılaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress–strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/????}{arXiv:????}
#' 
#' @examples
#' set.seed(123)
#' n <- 50
#' a1 <- 0.75; b1 <- 1.5; lambda1 <- 0.6
#' a2 <- 1.2; b2 <- 0.5; lambda2 <- 0.9
#' theta <- 1 # 2, 3, 4, 5
#' # data generation
#' dat <- SSReliabilityClaytonMWD::rMweibull_Clayton(n, a1, b1, lambda1, a2, b2, lambda2, theta)
#' str(dat)
#' 
#' @export
rMweibull_Clayton <- function(n, a1, b1, lambda1, a2, b2, lambda2, theta){
  # Generate Clayton Copula Uniforms (U, V)
  u <- runif(n)
  w <- runif(n)
  
  # Conditional copula method for Clayton
  v <- ( 1 + u^(-theta) * (w^(-theta/(theta + 1)) - 1) )^(-1/theta)
  
  eps <- 1e-7
  u <- pmax(pmin(u, 1 - eps), eps)
  v <- pmax(pmin(v, 1 - eps), eps)
  
  # Transform Uniforms to MWD Marginals
  x <- qMweibull(u, a1, b1, lambda1)
  y <- qMweibull(v, a2, b2, lambda2)
  
  return( list(U=u, V=v, X=x, Y=y) )
}
