# Rcpp::sourceCpp(file = "./src/code.cpp", rebuild = TRUE)
test_that("update dlm", {
  # rm(list = ls())
  # devtools::load_all()

  object <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = c(0.90, 0.95))
  )

  object$prior <- list("a" =  matrix(c(10, 1), ncol = 1),
                       "R" = diag(100, 2),
                       "W" = object[["D"]] * diag(100, 2),
                       "n" = 1, "s" = 1)

  tmp <- update_moments(object = object, y = 20)
  tmp2 <- update_moments_cpp.dlm(object = object, y = 20)

  out_R <- list(a = tmp$prior$a, R = tmp$prior$R,
                m = tmp$posterior$m, C = tmp$posterior$C,
                n = tmp$posterior$n, s = tmp$posterior$s)
  out_cpp <- update_dlm(
    y = 20,
    F = object[["FF"]],
    G = object[["GG"]],
    D = object[["D"]],
    a = object[["prior"]][["a"]],
    R = object[["prior"]][["R"]],
    n = 1, s = 1, df_variance = 1)
  expect_equal(drop(out_cpp$a), unname(out_R$a))
  expect_equal(out_cpp$R, unname(out_R$R))
  expect_equal(drop(out_cpp$m), unname(out_R$m))
  expect_equal(out_cpp$C, unname(out_R$C))
  expect_equal(out_cpp$n, out_R$n)
  expect_equal(out_cpp$s, out_R$s)

  # rbenchmark::benchmark(
  #   R = update_moments(object = object, y = NA_real_),
  #   cpp = update_moments_cpp.dlm(object = object, y = NA_real_),
  #   cpp_2 = update_dlm(
  #     y = 20,
  #     F = object[["FF"]],
  #     G = object[["GG"]],
  #     D = object[["D"]],
  #     a = object[["prior"]][["a"]],
  #     R = object[["prior"]][["R"]],
  #     n = 1, s = 1, df_variance = 1),
  #   replications = 100
  # )


})
















