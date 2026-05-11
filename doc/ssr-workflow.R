## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## -----------------------------------------------------------------------------
# after released in CRAN
# install.packages("SSReliabilityClaytonMWD")

## -----------------------------------------------------------------------------
# install.packages("remotes")
# remotes::install_github("fatihki/SSReliabilityClaytonMWD")

## -----------------------------------------------------------------------------
rm(list = ls())

# Load the package
library(SSReliabilityClaytonMWD)

# set seed
set.seed(123)
n <- 50
a1 <- 0.75; b1 <- 1.5; lambda1 <- 0.6
a2 <- 1.2; b2 <- 0.5; lambda2 <- 0.9
theta <- 3

# simulate data
dat <- rMweibull_Clayton(n, a1, b1, lambda1, a2, b2, lambda2, theta)

## -----------------------------------------------------------------------------
# true stress-strength reliability value
R_true <- Reliability_Clayton_MWD(a1, b1, lambda1, a2, b2, lambda2, theta)
R_true$value

## -----------------------------------------------------------------------------
fit <- fit.SSR.ClaytonMWD(
  data = dat,
  ACI = TRUE,
  bootstrap = TRUE,
  B = 10,
  seed = 2026,
  one.step = TRUE,
  alpha = 0.05
)

## -----------------------------------------------------------------------------
print(fit)

## -----------------------------------------------------------------------------
# Load example data from the package
data(TerkosDam)
data(OmerliDam)

real_data <- list(X = TerkosDam, Y = OmerliDam)

## -----------------------------------------------------------------------------
fit_ssr <- fit.SSR.ClaytonMWD(
  data = real_data,
  ACI = TRUE,
  bootstrap = TRUE,
  B = 10,
  seed = 2026,
  one.step = TRUE,
  alpha = 0.05
)

## -----------------------------------------------------------------------------
print(fit_ssr)

