## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## -----------------------------------------------------------------------------
rm(list = ls())

# Load the package
library(SSReliabilityClaytonMWD)

# generate data from MWD(a, b, lambda) with lambda = 0 reduces two-parameter Weibull distribution
n <- 50
a <- 0.75; b <- 1.25; lambda <- 0
set.seed(123)
dat <- rMweibull(n, a, b, lambda)
# random initial points 
init <- runif(2)

## -----------------------------------------------------------------------------
# Fit WD to `dat`
fit.mle <- fitWD(data = dat, est.method = "mle", opt.method = "L-BFGS-B", starts = init,
                  lower = c(1e-05,1e-05), upper = c(Inf,Inf), hessian = FALSE )
fit.mle$estimates

## -----------------------------------------------------------------------------
fit.lse <- fitWD(data = dat, est.method = "lse", opt.method = "L-BFGS-B", starts = init,
                   lower = c(1e-05,1e-05), upper = c(Inf,Inf), hessian = FALSE )
fit.lse$estimates

## -----------------------------------------------------------------------------
fit.wlse <- fitWD(data = dat, est.method = "wlse", opt.method = "L-BFGS-B", starts = init,
                    lower = c(1e-05,1e-05), upper = c(Inf,Inf), hessian = FALSE )
fit.wlse$estimates

## -----------------------------------------------------------------------------
fit.mps <- fitWD(data = dat, est.method = "mps", opt.method = "L-BFGS-B", starts = init,
                   lower = c(1e-05,1e-05), upper = c(Inf,Inf), hessian = FALSE )
fit.mps$estimates

