test_that("dgegm model definition works", {

  a0 <- c(5, 0.01, 1)
  R0 <- diag(500, ncol = 3, nrow = 3)
  mod <- dgegm(discount_factors = 0.95, lambda = 1,
               prior = list(a = a0, R = R0),
               variance_law = list(type = "poisson"))

  expect_s3_class(mod, "dgegm")

  mod <- dgegm(discount_factors = 0.95, lambda = 1,
               variance_law = list(type = "poisson"))
  expect_s3_class(mod, "dgegm")
  expect_null(mod$prior$a0)

})
