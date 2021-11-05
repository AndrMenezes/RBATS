test_that("dgegm model definition works", {

  mod <- dgegm(discount_factors = 0.95, lambda = 1,
               variance_law = list(type = "poisson"))
  expect_s3_class(mod, "dgegm")
})
