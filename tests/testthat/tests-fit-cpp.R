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

  plot(c(Nile))
  # lines(out$f[, 1], col = "red", lwd = 2)
  lines(out$m[1, ], col = "orange", lwd = 2)

  plot(out$m[2, ], type = "l", col =  "grey")
  lines(out$m[3, ], col = "blue")
  abline(h = 0, col = "red")

  dim(out$R)
  out$R[,,54]
})

