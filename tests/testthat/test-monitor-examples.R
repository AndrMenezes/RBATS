wd <- "../../python_packages/pybats-detection/src/pybats_detection/data"

cp6 <- read.csv(file.path(wd, "cp6__west_harrison.csv"))
usethis::use_data(cp6, overwrite = TRUE)


rm(list = ls())
devtools::load_all()
data("cp6")
plot(cp6$sales)
y <- cp6$sales

model <- dlm(
  polynomial_order = 2,
  discount_factors = list(polynomial = c(0.80))
)
D_exp <- model[["D"]]
diag(D_exp) <- 1/c(0.40, 0.80)
ffr <- forward_filter_dlm_monitor(
  y = y, F = model[["FF"]],
  G = model[["GG"]], D = model[["D"]],
  a = matrix(c(y[1], 0), ncol = 1),
  R = diag(1000, 2), n = 1, s = 1,
  df_variance = model[["df_variance"]],
  bf_threshold = 0.2,
  location_shift = 4,
  scale_shift = 1,
  exception_D = D_exp,
  verbose = TRUE,
  monitor_start = 6)

plot(y)
lines(ffr$f, col = "blue")
lines(ffr$f + 1.96 * sqrt(ffr$q), col = "blue", lty = 2)
lines(ffr$f - 1.96 * sqrt(ffr$q), col = "blue", lty = 2)
plot(ffr$H)

#################

telephone_calls <- read.csv(file.path(wd, "telephone_calls__pankratz.csv"))
usethis::use_data(telephone_calls, overwrite = TRUE)


rm(list = ls())
devtools::load_all()
data("telephone_calls")
y <- telephone_calls$average_daily_calls
plot(y)

model <- dlm(
  polynomial = list(order = 2, discount_factor = 0.90)
)
D_exp <- model[["D"]]
diag(D_exp) <- 1/c(0.20, 0.90)
ffr <- forward_filter_dlm_monitor_bilateral(
  y = y, F = model[["FF"]],
  G = model[["GG"]], D = model[["D"]],
  a = matrix(c(y[1], 0), ncol = 1),
  R = diag(100, 2), n = 1, s = 1,
  df_variance = model[["df_variance"]],
  bf_threshold = 0.135,
  location_shift = 4,
  scale_shift = 1,
  exception_D = D_exp,
  verbose = TRUE,
  monitor_start = 40)

plot(y)
lines(ffr$f, col = "blue")
lines(ffr$f + 1.96 * sqrt(ffr$q), col = "blue", lty = 2)
lines(ffr$f - 1.96 * sqrt(ffr$q), col = "blue", lty = 2)


###################
rm(list = ls())
devtools::load_all()
# market_share <- read.csv(file.path(wd, "pole_market_share.csv"))
# usethis::use_data(market_share, overwrite = TRUE)
data("market_share")
X <- as.matrix(market_share[, c("price", "prom", "cprom")])
# X <- t( t(X) - colMeans(X))
head(X)
y <- market_share$share

data_pybats_detection <- read.csv("~/Documents/paper_pybats_detection/replication/pybats_detection.csv")

model <- dlm(
  polynomial = list(order = 1, discount_factor = 1),
  regressor = list(xreg = X, discount_factor = 0.90)
)
R <- diag(4) * 4
R[1, 1] <- 25
D_exp <- model[["D"]]
diag(D_exp) <- 1/c(0.20, rep(0.90, 3))
ffr <- forward_filter_dlm_monitor_bilateral_X(
  y = y,
  F = model[["FF"]],
  X = t(model[["xreg"]]),
  G = model[["GG"]],
  D = model[["D"]],
  a = matrix(c(42, 0, 0, 0), ncol = 1),
  R = R,
  n = 1, s = 1,
  df_variance = model[["df_variance"]],
  bf_threshold = 0.135,
  location_shift = 4,
  scale_shift = 1,
  exception_D = D_exp,
  verbose = TRUE,
  monitor_start = 5)

data_1 <- data_pybats_detection[, c("t", "f", "q", "e", "H_lower", "H_upper")]

data_2 <- data.frame(f = ffr$f, q = ffr$q, e = (y - ffr$f) / sqrt(ffr$q),
                     H_lower = ffr$H_lower, H_upper = ffr$H_upper)

data_1[30:36, ]
data_2[30:36, ]



plot(y)
lines(ffr$f, col = "blue")
lines(ffr$f + 1.96 * sqrt(ffr$q), col = "blue", lty = 2)
lines(ffr$f - 1.96 * sqrt(ffr$q), col = "blue", lty = 2)



ffr <- forward_filter_dlm_X(
  y = y,
  F = model[["FF"]],
  X = t(model[["xreg"]]),
  G = model[["GG"]],
  D = model[["D"]],
  a = matrix(c(y[1], 0, 0, 0), ncol = 1),
  R = diag(4, 4),
  n = 1, s = 1,
  df_variance = model[["df_variance"]])














