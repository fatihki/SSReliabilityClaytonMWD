#' Stress–Strength Reliability under Clayton Copula with MWD Marginals
#'
#' @title Stress–Strength Reliability for MWD under Clayton Copula
#'
#' @description
#' Computes the stress–strength reliability (SSR)
#' \eqn{R = P(X > Y)}, where
#' \eqn{X \sim \text{MWD}(a_1, b_1, \lambda_1)} (strength) and
#' \eqn{Y \sim \text{MWD}(a_2, b_2, \lambda_2)} (stress),
#' with dependence modeled using a Clayton copula.
#'
#' @name Reliability_Clayton_MWD
#' 
#' @param a1,b1,lambda1 Parameters of the strength variable \eqn{X},
#' with \eqn{a_1 > 0}, \eqn{b_1 \ge 0}, and \eqn{\lambda_1 \ge 0}.
#'
#' @param a2,b2,lambda2 Parameters of the stress variable \eqn{Y},
#' with \eqn{a_2 > 0}, \eqn{b_2 \ge 0}, and \eqn{\lambda_2 \ge 0}.
#'
#' @param theta Clayton copula dependence parameter, \eqn{\theta > 0}.
#' 
#' @details
#' The stress–strength reliability is defined as
#' \eqn{R = P(X > Y)}, which can be expressed using the joint distribution
#' induced by the Clayton copula.
#'
#' In copula form, the reliability is computed as
#' \deqn{
#' R = \int_{0}^{\infty}
#' F_X(x)^{-(\theta + 1)}
#' \left( F_X(x)^{-\theta} + G_Y(x)^{-\theta} - 1 \right)^{-\left(\frac{1}{\theta}+1\right)}
#' f_X(x)\, dx,
#' }
#' which can also be written as
#' \deqn{
#' R = \int_{0}^{1}
#' t^{-(\theta + 1)}
#' \left( t^{-\theta} + G_Y(F_X^{-1}(t))^{-\theta} - 1 \right)^{-\left(\frac{1}{\theta}+1\right)}
#' dt,
#' }
#'
#' where \eqn{F_X(x) = 1 - \exp(-a_1 x^{b_1} e^{\lambda_1 x})}
#' and \eqn{G_Y(y) = 1 - \exp(-a_2 y^{b_2} e^{\lambda_2 y})}.
#' 
#' Further theoretical details can be found in  
#' \code{vignette("ssr-theory", package = "SSReliabilityClaytonMWD")}
#' and Kizilaslan (2026).
#' 
#' @return
#' A list identical to the output of \code{stats::integrate}:
#' \item{value}{Numerical value of the integral (reliability value).}
#' \item{abs.error}{Estimated absolute error of the numerical integration.}
#' 
#' @references
#' Kizilaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress--strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/2604.12130}{arXiv:2604.12130}
#' 
#' @examples
#' a1 <- 0.75; b1 <- 1.5; lambda1 <- 0.6
#' a2 <- 1.2; b2 <- 0.5; lambda2 <- 0.9
#' theta <- 1 # 2, 3, 4, 5
#' R <- Reliability_Clayton_MWD(a1, b1, lambda1, a2, b2, lambda2, theta)
#' R$value
#' # approximated reliability R based on MC method
#' R_MC <- Reliability_Clayton_MWD_MC(a1, b1, lambda1, a2, b2, lambda2, theta, N = 100000)
#' R_MC$R
#' 
#' @export
Reliability_Clayton_MWD <- function(a1, b1, lambda1, a2, b2, lambda2, theta ) {
  integrand <- function(x) {
    Finverse_u <- qMweibull(x, a1, b1, lambda1, TRUE)
    G <- pMweibull(Finverse_u, a2, b2, lambda2)
    term <- (x^(-theta) + G^(-theta) - 1)^(-(1 + theta)/theta)
    return( x^(-(theta + 1)) * term )
  }
  
  I <- stats::integrate(integrand, lower=0, upper=1)
  return(I)
}
#'
#' an alternative version of "Reliability_Clayton_MWD.R" to be able to evaluate numeric gradient of R for the ACI of R.
#' @noRd
R_Clayton_MWD <- function(par) {
  a1 = par[1]; b1 = par[2]; lambda1 = par[3]; a2 = par[4]; b2= par[5]; lambda2 = par[6]; theta = par[7]
  integrand <- function(x) {
    Finverse_u <- qMweibull(x, a1, b1, lambda1, TRUE)
    G <- pMweibull(Finverse_u, a2, b2, lambda2)
    term <- (x^(-theta) + G^(-theta) - 1)^(-(1 + theta)/theta)
    return( x^(-(theta + 1)) * term )
  }
  I <- stats::integrate(integrand, lower=0, upper=1)$value
  return(I)
}
#' alternative calculation with infinite upper limit in the integral
#' @noRd
Reliability_Clayton_MWD_inf <- function(a1, b1, lambda1, a2, b2, lambda2, theta,
                                       lower = 0, upper = Inf) {
  
  integrand <- function(x) {
    u  <- pMweibull(x, a1, b1, lambda1)   # F_X(x)
    v  <- pMweibull(x, a2, b2, lambda2)   # F_Y(x)
    fx <- dMweibull(x, a1, b1, lambda1)   # f_X(x)
    
    return( fx * (u^(-theta-1)) * ( u^(-theta) +  v^(-theta) - 1)^(-(1 + theta)/theta) )
  }
  
  I <- stats::integrate(integrand, lower = lower, upper = upper)
  return(I)
}
#'
#' Monte Carlo Estimation of Stress–Strength Reliability under Clayton Copula
#'
#' @title Monte Carlo Estimation of Stress–Strength Reliability (MWD + Clayton Copula)
#' 
#' @name Reliability_Clayton_MWD_MC
#' 
#' @description
#' Computes an approximate value of the stress–strength reliability
#' \eqn{R = P(X > Y)} using Monte Carlo integration method under the modified Weibull
#' model with dependence induced by a Clayton copula.
#' 
#' It calculates an approximate value of \eqn{R} using the Monte Carlo integration method based
#' on a random generated MWD sample from \eqn{MWD(a_1, b_1, \lambda_1)}.
#' 
#' @param a1,b1,lambda1 Parameters of the strength variable \eqn{X},
#' with \eqn{a_1 > 0}, \eqn{b_1 \ge 0}, and \eqn{\lambda_1 \ge 0}.
#'
#' @param a2,b2,lambda2 Parameters of the stress variable \eqn{Y},
#' with \eqn{a_2 > 0}, \eqn{b_2 \ge 0}, and \eqn{\lambda_2 \ge 0}.
#'
#' @param theta Clayton copula dependence parameter, \eqn{\theta > 0}.
#' 
#' @param N Integer. Number of Monte Carlo samples from \eqn{X \sim MWD(a_1,b_1,\lambda_1)}.
#' 
#' @details
#' The approximate stress–strength reliability \eqn{R} is approximated via Monte
#' Carlo integration:
#'
#' \deqn{
#' R \approx \frac{1}{N} \sum_{i=1}^{N}
#' T(x_i; \boldsymbol{\Omega_1}, \boldsymbol{\Omega_2}, \theta) = \tilde{R},
#' }
#'
#' where \eqn{\boldsymbol{\Omega_1} = (a_1, b_1, \lambda_1)} and
#' \eqn{\boldsymbol{\Omega_2} = (a_2, b_2, \lambda_2)} denote the marginal parameter vectors.
#' The function
#'
#' \deqn{
#' T(x; \boldsymbol{\Omega_1}, \boldsymbol{\Omega_2}, \theta)
#' = F_X(x)^{-(\theta + 1)}
#' \left( F_X(x)^{-\theta} + G_Y(x)^{-\theta} - 1 \right)^{-\left(\frac{1}{\theta}+1\right)},
#' }
#'
#' where \eqn{x_1, \dots, x_N} is a random sample generated from the
#' Modified Weibull distribution \eqn{MWD(a_1, b_1, \lambda_1)}.
#'
#' The accuracy of the approximation improves as sample size \eqn{N} increases.
#' In practice, values around \eqn{N = 10^5}  typically provide stable results.
#' 
#' Further details can be found in Kizilaslan (2026).
#' 
#' @return A list containing:
#' \item{R}{Approximated reliability value based on \eqn{N} generated random samples.}
#' \item{R_MCsample}{Monte Carlo approximate of the reliability \eqn{R} based on \eqn{N} generated samples.}
#' \item{sample}{Generated random sample from \eqn{MWD(a_1, b_1, \lambda_1)} with \eqn{N} size.}
#' 
#' @references
#' Kizilaslan, F. (2026).
#' \emph{Reliability estimation in dependent stress–strength model with Clayton copula and modified Weibull margins}.
#' \href{https://arxiv.org/abs/2604.12130}{arXiv:2604.12130}
#' 
#' @examples
#' # example code
#' a1 <- 0.75; b1 <- 1.5; lambda1 <- 0.6
#' a2 <- 1.2; b2 <- 0.5; lambda2 <- 0.9
#' theta <- 1 # 2, 3, 4, 5
#' R <- Reliability_Clayton_MWD(a1, b1, lambda1, a2, b2, lambda2, theta)
#' R$value
#' # approximated reliability R based on MC method
#' R_MC <- Reliability_Clayton_MWD_MC(a1, b1, lambda1, a2, b2, lambda2, theta, N = 100000)
#' R_MC$R
#' 
#' @export
Reliability_Clayton_MWD_MC <- function(a1, b1, lambda1, a2, b2, lambda2, theta, N = 10000) {
  
  x <- rMweibull(N, a1, b1, lambda1)
  u <- pMweibull(x, a1, b1, lambda1)  
  v <- pMweibull(x, a2, b2, lambda2)   
  MC_sample <- (u^(-theta-1)) * ( u^(-theta) +  v^(-theta) - 1)^(-(1 + theta)/theta) 
  
  return(list( R = mean(MC_sample), R_MCsample = MC_sample, sample = x) )
}
