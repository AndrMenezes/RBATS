devtools::load_all()
set.seed(6669)
y <- c(rnorm(40, mean = 80, sd = 0.5),
       rnorm(40, mean = 90, sd = 0.5),
       0.5 * 1:40 + rnorm(40, mean = 90, sd = 0.3))
y[20] <- 95

plot(y)


model <- dlm(
  polynomial_order = 2,
  discount_factors = list(polynomial = c(0.95, 0.90))
)

ffr <- forward_filter_dlm_monitor(y = y, F = model[["FF"]],
                                  G = model[["GG"]], D = model[["D"]],
                                  a = matrix(c(y[1], 0), ncol = 1),
                                  R = diag(100, 2), n = 1, s = 1,
                                  df_variance = model[["df_variance"]],
                                  bf_threshold = 0.135,
                                  location_shift = 4,
                                  scale_shift = 2,
                                  exception_D = diag(1/c(0.10, 0.20), 2))
plot(y)
lines(ffr$f, col = "blue")
lines(ffr$f + 1.96 * sqrt(ffr$q), col = "blue", lty = 2)
lines(ffr$f - 1.96 * sqrt(ffr$q), col = "blue", lty = 2)
