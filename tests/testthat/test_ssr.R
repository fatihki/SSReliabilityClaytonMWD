test_that("SSR ClaytonMWD model runs without error", {
  
  data <- list(X = TerkosDam, Y = OmerliDam)
  
  fit <- fit.SSR.ClaytonMWD(
    data,
    ACI = TRUE,
    bootstrap = TRUE,
    B = 5,       
    seed = 2026,
    one.step = TRUE,
    alpha = 0.05
  )
  
  expect_true(!is.null(fit))
  expect_true(is.list(fit))
})