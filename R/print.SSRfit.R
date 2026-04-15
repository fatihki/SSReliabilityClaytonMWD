#' @title Print SSR.ClaytonMWD Fit Results
#' @description
#' Prints the results of an object of class \code{SSRfit} created by 
#' \code{fit.SSR.ClaytonMWD}, including separate tables of parameter estimates 
#' and interval estimates (MLE, LSE, WLSE, and MPS) for the model parameters, 
#' the Clayton copula parameter \eqn{\theta}, and the reliability \eqn{R}.
#' 
#' @importFrom knitr kable
#' @method print SSRfit
#'
#' @param x An object of class \code{SSRfit} obtained from 
#' \code{fit.SSR.ClaytonMWD}.
#' @param ... Additional arguments passed to the print method.
#' For example, \code{digits} controls the number of decimal places
#' used when printing numerical results.
#'
#' @examples
#' data = list(X = TerkosDam, Y = OmerliDam)
#' fit.SSR = fit.SSR.ClaytonMWD(data, ACI = TRUE, bootstrap = TRUE, B = 10,
#'                              seed = 2026, one.step = TRUE, alpha = 0.05)
#' print(fit.SSR)
#' print(fit.SSR, 3)
#' @export
print.SSRfit <- function(x, ...) {
  
  args <- list(...)
  digits <- args$digits
  if (is.null(digits)) digits <- 5
  
  cat("\nFitting Results of Stress-Strength Reliability Model with MWD Marginals via Clayton Copula \n")
  
  print(knitr::kable(x$all.results, caption = "Parameter Estimates", digits = digits))
  
  if (!is.null(x$ACI.parameters)) {
    print(knitr::kable(x$ACI.parameters,
                       caption = paste0(100*(1-x$alpha), "% Asymptotic Confidence Intervals (MLE)"),
                       digits = digits))
  }
  
  if (!is.null(x$boot.mle)) {
    print(knitr::kable(x$boot.mle,
                       caption = paste0(100*(1-x$alpha), "% Bootstrap CIs (MLE)"),
                       digits = digits))
  }
  
  if (!is.null(x$boot.lse)) {
    print(knitr::kable(x$boot.lse,
                       caption = paste0(100*(1-x$alpha), "% Bootstrap CIs (LSE)"),
                       digits = digits))
  }
  
  if (!is.null(x$boot.wlse)) {
    print(knitr::kable(x$boot.wlse,
                       caption = paste0(100*(1-x$alpha), "% Bootstrap CIs (WLSE)"),
                       digits = digits))
  }
  
  if (!is.null(x$boot.mps)) {
    print(knitr::kable(x$boot.mps,
                       caption = paste0(100*(1-x$alpha), "% Bootstrap CIs (MPS)"),
                       digits = digits))
  }
  
  if (!is.null(x$theta.Ktau)) {
    cat("\nKendall's tau estimate for theta:", round(x$theta.Ktau, digits = digits), "\n")
  }
  
  
  invisible(x)
}
