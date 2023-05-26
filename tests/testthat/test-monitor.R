devtools::load_all()
set.seed(6669)
y <- c(rnorm(40, mean = 80, sd = 0.5),
       rnorm(40, mean = 90, sd = 0.5),
       2 * 1:40 + rnorm(40, mean = 90, sd = 0.25))
y[20:30] <- NA_real_
y[38:40] <- NA_real_

# plot(y)


model <- dlm(
  polynomial_order = 2,
  discount_factors = list(polynomial = c(0.97, 0.97))
)

ffr <- forward_filter_dlm_monitor(y = y, F = model[["FF"]],
                                  G = model[["GG"]], D = model[["D"]],
                                  a = matrix(c(y[1], 0), ncol = 1),
                                  R = diag(100, 2), n = 1, s = 1,
                                  df_variance = model[["df_variance"]],
                                  bf_threshold = 0.135,
                                  location_shift = 4,
                                  scale_shift = 1,
                                  exception_D = diag(1/c(0.10, 0.30), 2),
                                  verbose = TRUE)
plot(y)
lines(ffr$f, col = "blue")
lines(ffr$f + 1.96 * sqrt(ffr$q), col = "blue", lty = 2)
lines(ffr$f - 1.96 * sqrt(ffr$q), col = "blue", lty = 2)

plot(ffr$m[2, ], col = "blue")

par(mfrow = c(3, 1))
plot(log(ffr$H))
plot(ffr$L)
plot(ffr$l)
graphics.off()

###############################################################################
# bilateral
devtools::load_all()
set.seed(6669)
y <- c(rnorm(40, mean = 80, sd = 0.5),
       rnorm(40, mean = 70, sd = 0.5),
       c(1, 2, rep(3, 10)) + rnorm(12, mean = 70, sd = 0.25))
# y[20] <- NA_real_

plot(y)


model <- dlm(
  polynomial_order = 1,
  discount_factors = list(polynomial = c(0.97))
)

ffr <- forward_filter_dlm_monitor_bilateral(
  y = y, F = model[["FF"]],
  G = model[["GG"]], D = model[["D"]],
  a = matrix(y[1], ncol = 1),
  R = diag(1, 1), n = 1, s = 1,
  df_variance = model[["df_variance"]],
  bf_threshold = 0.135,
  location_shift = 4,
  scale_shift = 1,
  exception_D = diag(1/c(0.10), 1),
  verbose = TRUE)

plot(y)
lines(ffr$f, col = "blue")
lines(ffr$f + 1.96 * sqrt(ffr$q), col = "blue", lty = 2)
lines(ffr$f - 1.96 * sqrt(ffr$q), col = "blue", lty = 2)

###############################################################################
# with covariates




