#' Bivariate Random Data Generation under Clayton Copula with MWD Marginals
#'
#' @title Random Generation for MWD Marginals via Clayton Copula
#'
#' @description
#' Generates bivariate random samples from a dependent stress–strength model where
#' both marginals follow the Modified Weibull Distribution (MWD), and the dependence
#' structure between the variables is modeled using a Clayton copula.
#' 
#' Generates bivariate random samples from a dependent stress–strength model where
#' both marginals follow the Modified Weibull Distribution (MWD), and dependence
#' between variables is modeled using a Clayton copula.
#' 
#' @import stats
#' 
#' @name rMweibull_Clayton 
#'
#' @param n Integer. Number of observations to be generated.
#'
#' @param a1,b1,lambda1 Parameters of the strength variable \eqn{X},
#' with \eqn{a_1 > 0}, \eqn{b_1 \ge 0}, and \eqn{\lambda_1 \ge 0}.
#'
#' @param a2,b2,lambda2 Parameters of the stress variable \eqn{Y},
#' with \eqn{a_2 > 0}, \eqn{b_2 \ge 0}, and \eqn{\lambda_2 \ge 0}.
#'
#' @param theta Clayton copula dependence parameter, \eqn{\theta > 0}.
#' 
#' 
#' @details
#' This function generates dependent uniform variables using the Clayton copula,
#' which are then transformed via inverse CDFs of the Modified Weibull marginals
#' to obtain \eqn{(X, Y)}.
#'
#' Further details are provided in Kizilaslan (2026).
#' 
#' @return A list containing:
#' @return A list containing:
#' \item{U}{Uniform samples used in the copula construction.}
#' \item{V}{Dependent uniform samples generated via the Clayton copula.}
#' \item{X}{Simulated observations from \eqn{X \sim \mathrm{MWD}(a_1, b_1, \lambda_1)} obtained by transforming \eqn{U}.}
#' \item{Y}{Simulated observations from \eqn{Y \sim \mathrm{MWD}(a_2, b_2, \lambda_2)} obtained by transforming \eqn{V}.}
#' 
#' @references
#' Kizilaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress–strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/2604.12130}{arXiv:2604.12130}
#' 
#' @examples
#' set.seed(123)
#' n <- 50
#' a1 <- 0.75; b1 <- 1.5; lambda1 <- 0.6
#' a2 <- 1.2; b2 <- 0.5; lambda2 <- 0.9
#' theta <- 1 # 2, 3, 4, 5
#' # data generation
#' dat <- rMweibull_Clayton(n, a1, b1, lambda1, a2, b2, lambda2, theta)
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
