test_that("Nile data", {
  y <- c(Nile)

  # Define the model
  model <- dlm(
    polynomial_order = 1,
    discount_factors = list(polynomial = 0.75)
  )

  out <- forward_filter_cpp.dlm(model = model, y = y,
                                a = matrix(c(y[1]), ncol = 1),
                                R = diag(100, 1))

  fitNile <- StructTS(Nile, "level")

  # rbenchmark::benchmark(
  #   structs = StructTS(Nile, "level"),
  #
  #   rbats = forward_filter_cpp.dlm(model = model, y = y,
  #                              a = matrix(c(y[1]), ncol = 1),
  #                              R = diag(100, 1)),
  #   R = forward_filter.dlm(object = model, y = y,
  #                          a0 = matrix(c(y[1]), ncol = 1),
  #                          R0 = diag(100, 1)),
  #   replications = 10
  # )

  out$s[length(out$s)]
  out$C[length(out$C)]
  plot(c(Nile))
  # lines(out$f[, 1], col = "red", lwd = 2)
  lines(out$m[1, ], col = "orange", lwd = 2)
  lines(c(fitted(fitNile)), col = "blue", lwd = 2)

  plot(out$s)
  abline(h = fitNile$coef[2L], col = "blue")

  plot(out$C)
  abline(h = fitNile$coef[1L], col = "blue")


})

test_that("Nile data with fake covariates", {
  y <- c(Nile)

  # Fake covariate
  X <- matrix(rnorm(2 * length(y)), ncol = 2)
  # Define the model
  model <- dlm(
    polynomial_order = 1, xreg = X,
    discount_factors = list(polynomial = 0.90, regressors = 0.98)
  )

  out <- forward_filter_dlm_X(y = y, F = model[["FF"]], X = t(X),
                              G = model[["GG"]], D = model[["D"]],
                              a = matrix(c(y[1], 0, 0), ncol = 1),
                              R = diag(100, 3), n = 1, s = 1,
                              df_variance = model[["df_variance"]])

  bck <- backward_smoother_dlm_X(F = model[["FF"]], G = model[["GG"]],
                                 X = t(X),
                                 m_seq = out$m, a_seq = out$a, C_seq = out$C,
                                 R_seq = out$R)

  nrow(out$m)


  plot(c(Nile))
  lines(out$m[1, ], col = "blue", lwd = 2)
  lines(bck$ak[1, ], col = "red", lwd = 2)

  plot(out$m[2, ], type = "l", col =  "grey")
  lines(out$m[3, ], col = "blue")
  abline(h = 0, col = "red")

  dim(out$R)
  out$R[,,54]
})

test_that("Nile data with level+growth", {
  # devtools::load_all()
  y <- c(Nile)

  # Define the model
  model <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = 0.90)
  )

  out <- forward_filter_dlm(y = y, F = model[["FF"]],
                            G = model[["GG"]], D = model[["D"]],
                            a = matrix(c(y[1], 0), ncol = 1),
                            R = diag(100, 2), n = 1, s = 1,
                            df_variance = model[["df_variance"]])
  # devtools::load_all()
  bck <- backward_smoother_dlm(F = model[["FF"]], G = model[["GG"]],
                               m_seq = out$m, a_seq = out$a, C_seq = out$C,
                               R_seq = out$R)
  fwf <- forward_filter.dlm(object = model, y = y,
                            a0 = matrix(c(y[1], 0), ncol = 1),
                            R0 = diag(100, 2))
  bck_R <- backward_smoother.dlm(object = model, state_parameters = fwf$state_parameters)

  # rbenchmark::benchmark(
  #
  #   cpp =  backward_smoother_dlm(F = model[["FF"]], G = model[["GG"]],
  #                                  m_seq = out$m, a_seq = out$a, C_seq = out$C,
  #                                  R_seq = out$R),
  #   R = backward_smoother.dlm(object = model, state_parameters = fwf$state_parameters),
  #   replications = 10
  # )

  x11()
  plot(c(Nile), ylim = c(400, 1500))
  lines(out$f[, 1], col = "red", lwd = 2)
  lines(out$f[, 1] + 1.65*sqrt(out$q[, 1]), col = "red", lwd = 2, lty = 2)
  lines(out$f[, 1] - 1.65*sqrt(out$q[, 1]), col = "red", lwd = 2, lty = 2)
  lines(bck$fk[, 1], col = "blue", lwd = 2)
  lines(bck$fk[, 1] + 1.65*sqrt(bck$qk[, 1]), col = "blue", lwd = 2, lty = 2)
  lines(bck$fk[, 1] - 1.65*sqrt(bck$qk[, 1]), col = "blue", lwd = 2, lty = 2)

  x11()
  plot(c(Nile), ylim = c(400, 1500))
  lines(out$m[1, ], col = "red", lwd = 2)
  lines(out$m[1, ] + 1.65*sqrt(out$C[1,1, ]), col = "red", lwd = 2, lty = 2)
  lines(out$m[1, ] - 1.65*sqrt(out$C[1,1, ]), col = "red", lwd = 2, lty = 2)
  lines(bck$ak[1, ], col = "blue", lwd = 2)
  lines(bck$ak[1, ] + 1.65*sqrt(bck$Rk[1,1, ]), col = "blue", lwd = 2, lty = 2)
  lines(bck$ak[1, ] - 1.65*sqrt(bck$Rk[1,1, ]), col = "blue", lwd = 2, lty = 2)

  x11()
  plot(c(Nile), ylim = c(400, 1500))
  lines(bck_R$data_predictive$mean, col = "red", lwd = 2)
  lines(bck_R$data_predictive$ci_lower, col = "red", lwd = 2, lty = 2)
  lines(bck_R$data_predictive$ci_upper, col = "red", lwd = 2, lty = 2)
  lines(bck$fk[, 1], col = "blue", lwd = 2)
  lines(bck$fk[, 1] + 1.96*sqrt(bck$qk[, 1]), col = "blue", lwd = 2, lty = 2)
  lines(bck$fk[, 1] - 1.96*sqrt(bck$qk[, 1]), col = "blue", lwd = 2, lty = 2)

  # bck_R$data_predictive$mean[1:10]
  # bck$fk[1:10, 1]
  # S_t <- out$s[100]
  # plot(S_t / out$s[-100])



})

