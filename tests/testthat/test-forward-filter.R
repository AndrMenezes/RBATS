test_that("forward filter dlm", {

  rm(list = ls())
  devtools::load_all()
  object <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = c(0.90, 0.95))
  )
  object$prior <- list("a" =  matrix(c(1000, 1), ncol = 1),
                       "R" = diag(100, 2),
                       "W" = (1 - object[["D"]]) * diag(100, 2),
                       "n" = 1, "s" = 1)
  Y <- as.numeric(Nile)
  Y[c(10:12, 45)] <- NA_real_
  out_R <- forward_filter(object = object, y = Y,
                          a0 = matrix(c(1000, 1), ncol = 1),
                          R0 = diag(100, 2),
                          n = 1, s = 1, cpp = FALSE)


  out_cpp <- forward_filter_cpp(model = object, y = Y,
                            a = matrix(c(1000, 1), ncol = 1),
                            R = diag(100, 2),
                            n = 1, s = 1)
  rbenchmark::benchmark(
    cpp_list = forward_filter_cpp(model = object, y = Y,
                                  a = matrix(c(1000, 1), ncol = 1),
                                  R = diag(100, 2),
                                  n = 1, s = 1),

    R = forward_filter(object = object, y = Y,
                       a0 = matrix(c(1000, 1), ncol = 1),
                       R0 = diag(100, 2),
                       n = 1, s = 1, cpp = FALSE),
    replications = 10
  )


  plot(Y)
  lines(out_cpp$m[1, ], col = "blue")
  lines(out_cpp$f, col = "red")

  lines(out_R$data_posterior[out_R$data_posterior$parameter == "level", ]$mean,
        col = "green")
  lines(out_R$data_predictive$mean, col = "yellow")


  head(out_cpp$data_predictive, 14)
  head(out_R$data_predictive, 14)
  head(out_cpp$data_posterior, 5)
  head(out_R$data_posterior, 5)
  out_R$data_posterior |>
    dplyr::filter(parameter != "level") |> head(14)
  out_cpp$data_posterior |>
    dplyr::filter(parameter != "level") |> head(14)

  plot(Y)
  lines(out_cpp$data_predictive$mean, col = "blue")
  lines(out_R$data_predictive$mean, col = "red")

  Y <- rnorm(1e4)
  object <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = c(0.90))
  )
  rbenchmark::benchmark(
    R = forward_filter(object = object, y = Y,
                       a0 = matrix(0, nrow = 2),
                       R0 = diag(100, 2),
                       n = 1, s = 1, cpp = FALSE),
    cpp = forward_filter(object = object, y = Y,
                         a0 = matrix(0, nrow = 2),
                         R0 = diag(100, 2),
                         n = 1, s = 1, cpp = TRUE),
    replications = 10
  )

})

Y <- rnorm(1e4)
object <- dlm(
  polynomial_order = 2,
  discount_factors = list(polynomial = c(0.90))
)
rbenchmark::benchmark(
  cpp_list = forward_filter_dlm(
    y = Y,
    F = object$FF,
    G = object$GG,
    D = object$D,
    a = matrix(0, nrow = 2),
    R = diag(100, 2),
    n = 1, s = 1, df_variance = 1),

  cpp_cube = forward_filter_dlm_cub(
    y = Y,
    F = object$FF,
    G = object$GG,
    D = object$D,
    a = matrix(0, nrow = 2),
    R = diag(100, 2),
    n = 1, s = 1, df_variance = 1),
  R = forward_filter(object = object, y = Y,
                     a0 = matrix(0, nrow = 2),
                     R0 = diag(100, 2),
                     n = 1, s = 1, cpp = FALSE),
  cpp_R = forward_filter(object = object, y = Y,
                     a0 = matrix(0, nrow = 2),
                     R0 = diag(100, 2),
                     n = 1, s = 1, cpp = TRUE),
  replications = 10
)










