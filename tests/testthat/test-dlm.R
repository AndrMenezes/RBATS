test_that("polynomial model works", {
  (m1 <- dlm(polynomial = list(order = 2, discount_factor = c(0.95, 0.98))))
  expect_equal(unname(diag(m1$D)), 1/c(0.95, 0.98))

  (m2 <- dlm(polynomial = list(order = 3,
                               discount_factor = c(0.95, 0.98, 0.99))))
  expect_equal(unname(diag(m2$D)), 1/c(0.95, 0.98, 0.99))


})

test_that("seasonal free model works", {

  p <- 3
  (m1 <- dlm(polynomial = list(order = 1, discount_factor = 0.95),
             seasonal = list(type = "free", period = p,
                             discount_factor = 0.98)))
  expect_equal(unname(diag(m1$D)), 1/c(0.95, rep(0.98, p - 1)))

  (m2 <- dlm(polynomial = list(order = 2, discount_factor = 0.95),
             seasonal = list(type = "free", period = p,
                             discount_factor = 0.98)))
  expect_equal(unname(diag(m2$D)), 1/c(0.95, 0.95, rep(0.98, p - 1)))

})

test_that("seasonal fourier model works", {

  (m1 <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 7, harmonics = 1:3,
                    discount_factor = 0.98)
  ))

  expect_equal(nrow(m1$FF), 7L)

  (m2 <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 7, harmonics = 1:3,
                    discount_factor = 0.98)
  ))
  expect_equal(nrow(m2$FF), 8L)

  (m3 <- dlm(
    polynomial = list(order = 3, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 3, harmonics = 1,
                    discount_factor = 0.98)
  ))
  expect_equal(nrow(m3$FF), 5L)

})

test_that("dynamic linear regression model works",{

  (m1 <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    regressor = list(xreg = matrix(runif(10), ncol = 2), discount_factor = 0.99)
  ))
  expect_equal(nrow(m1$GG), 3L)
  expect_equal(m1$i_polynomial, 1L)
  expect_equal(m1$i_regressor, c(2L, 3L))

  (m2 <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 3, harmonics = 1,
                    discount_factor = 0.98),
    regressor = list(xreg = matrix(runif(10), ncol = 2), discount_factor = 0.99)
  ))
  expect_equal(nrow(m2$GG), 6L)
  expect_equal(m2$i_polynomial, c(1L, 2L))
  expect_equal(m2$i_seasonal, c(3L, 4L))
  expect_equal(m2$i_regressor, c(5L, 6L))

})

test_that("dlm with trend + seasonal + regression + autoregressive",{

  (m1 <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 3, harmonics = 1,
                    discount_factor = 0.98),
    regressor = list(xreg = matrix(runif(10), ncol = 2), discount_factor = 0.99),
    autoregressive = list(order = 2, discount_factor = 0.998)
  ))
  expect_equal(rownames(m1$FF), c("level", "period_3__fourier_cos_1",
                                  "period_3__fourier_sin_1", "regressor_1",
                                  "regressor_2", "xi_1", "xi_2", "phi_1",
                                  "phi_2"))
  expect_equal(nrow(m1$GG), 9L)
  expect_equal(m1$i_polynomial, 1L)
  expect_equal(m1$i_seasonal, c(2L, 3L))
  expect_equal(m1$i_regressor, c(4L, 5L))
  expect_equal(m1$i_autoregressive, 6L:9L)

})

test_that("dlm with transfer function", {
  # devtools::load_all()
  m <- dlm(
    transfer_function = list(order = 2, xreg = matrix(rnorm(10), ncol = 1),
                             discount_factor = 0.998)
  )

  expect_equal(m$n_parms, 0L)
  expect_equal(m$tf_order, 2)
  expect_equal(m$i_transfer_function, 1:5)
  expect_equal(unname(m$FF), matrix(c(1, rep(0, 4)), ncol = 1))
  expect_equal(unname(m$GG),
               matrix(c(rep(NA, 5),
                        c(1, rep(0, 4)),
                        c(rep(0, 2), 1, rep(0, 2)),
                        c(rep(0, 3), 1, 0),
                        c(rep(0, 4), 1)), ncol = 5, byrow = T))

  # Several components + transfer function
  (m1 <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 3, harmonics = 1,
                    discount_factor = 0.98),
    autoregressive = list(order = 2, discount_factor = 0.998),
    transfer_function = list(order = 2, xreg = matrix(rnorm(10), ncol = 1),
                             discount_factor = 0.998)
  ))
  expect_equal(m1$i_polynomial, 1L)
  expect_equal(m1$i_seasonal, c(2L, 3L))
  expect_equal(m1$i_autoregressive, 4L:7L)
  expect_equal(m1$i_transfer_function, 8L:12L)
  expect_equal(nrow(m1$FF), 12)
  expect_equal(dim(m1$GG), c(12, 12))


  # Several components + TF(1) with fixed value for lambda
  m <- dlm(
    transfer_function = list(order = 1, xreg = matrix(rnorm(10), ncol = 1),
                             lambda = 0.5, discount_factor = 0.998)
  )
  expect_equal(m$fixed_tf_parm, 1)
  expect_equal(m$n_parms, 0)
  expect_equal(nrow(m$FF), 2)

  # TF(1) with fixed value for lambda
  m <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 3, harmonics = 1,
                    discount_factor = 0.98),
    autoregressive = list(order = 2, discount_factor = 0.998),
    transfer_function = list(order = 1, xreg = matrix(rnorm(10), ncol = 1),
                             lambda = 0.5, discount_factor = 0.998)
  )
  expect_equal(m$i_transfer_function, c(8, 9))
  expect_equal(m$fixed_tf_parm, 1)
  expect_equal(m$n_parms, 3)


})

test_that("dlm with cycle",{

  (m1 <- dlm(cycle = list(freq = 2.5, discount_factor = 0.998)))
  expect_equal(nrow(m1$GG), 2L)
  expect_equal(m1$i_cycle, 1:2)

  (m2 <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    cycle = list(freq = 2.5, discount_factor = 0.998, rho = 0.90)))
  expect_equal(nrow(m2$GG), 4L)
  expect_equal(m2$i_cycle, 3:4)

})
