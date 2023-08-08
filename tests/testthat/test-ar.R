rm(list = ls())
devtools::load_all()

y <- c(AirPassengers)

mod <- dlm(
  polynomial = list(order = 2, discount_factor = 0.95),
  seasonal = list(type = "fourier", period = 12, harmonics = 1:2,
                  discount_factor = 0.99)
)

out <- forward_filter_dlm_ar_2(y = y, F = mod$FF, G = mod$GG, D = mod$D,
                               m = matrix(c(y[1L], rep(0, 5)), ncol = 1),
                               C = diag(1000, nrow = 6), n = 1, s = 1,
                               df_variance = 1, ar_order = 0L,
                               n_parms = 6L)

plot(y)
lines(out$f[, 1L], col = "blue")
lines(out$m[1L, ], col = "red")

# AR(2) -------------------------------------------------------------------

# Simulate the data
nobs <- 200
sd_y <- 0.2
true_phi_1 <- 0.9
true_phi_2 <- -0.5

xi_1 <- numeric(nobs)
xi_2 <- numeric(nobs)
y <- numeric(nobs)
xi_1[1L] <- 0

# Random errors
set.seed(121)
nu <- rnorm(n = nobs, sd = sd_y)

# First observation
y[1L] <- xi_1[1L] + nu[1L]

for (t in seq_len(nobs)[-1]) {
  xi_1[t] <- true_phi_1 * xi_1[t - 1] + true_phi_2 * xi_2[t - 1] + nu[t]
  xi_2[t] <- xi_1[t - 1]
  y[t] <- xi_1[t]
}

mod <- dlm(autoregressive = list(order = 2, discount_factor = 0.998))
mod$GG
m0 <- matrix(c(0, 0, 1, 0), ncol = 1)
C0 <- diag(x = c(2, 2, 1, 1), nrow = 4)

out <- forward_filter_dlm(y = y, F = mod$FF, G = mod$GG, D = mod$D,
                          m = m0,
                          C = C0,
                          n = 1, s = 1, df_variance = 1,
                          ar_order = 2L,
                          n_parms = 0L)
plot(out$m[3L, ]); abline(h = true_phi_1)
plot(out$m[4L, ]); abline(h = true_phi_2)

# level + AR(2) -----------------------------------------------------------

# Simulate the data
nobs <- 100
sd_y <- 1
sd_mu <- 0.06
true_phi_1 <- 0.8
true_phi_2 <- -0.5

mu <- numeric(nobs)
xi_1 <- numeric(nobs)
xi_2 <- numeric(nobs)
y <- numeric(nobs)
xi_1[1L] <- 0
mu[1L] <- 5

# Random errors
set.seed(1111)
nu <- rnorm(n = nobs, sd = sd_y)
omega_mu <- rnorm(n = nobs, sd = sd_mu)

# First observation
y[1L] <- mu[1L] + xi_1[1L] + nu[1L]

for (t in seq_len(nobs)[-1]) {
  xi_1[t] <- true_phi_1 * xi_1[t - 1] + true_phi_2 * xi_2[t - 1] + nu[t]
  xi_2[t] <- xi_1[t - 1]
  mu[t] <- mu[t - 1] + omega_mu[t]
  y[t] <- mu[t] + xi_1[t]
}
plot(y)


# Model definition

mod <- dlm(
  polynomial = list(order = 1, discount_factor = 0.95),
  autoregressive = list(order = 2, discount_factor = 0.998)
)
m = matrix(c(y[1L], .5, 0, 1, 0), ncol = 1)

out <- forward_filter_dlm(y = y, F = mod$FF, G = mod$GG, D = mod$D,
                          m = matrix(c(y[1L], rep(0, 4)), ncol = 1),
                          C = diag(c(10, 2, 2, 1, 1), nrow = 5),
                          n = 1, s = 1, df_variance = 1,
                          ar_order = 2L,
                          n_parms = 1L)
plot(out$m[4L, ]); abline(h = true_phi_1)
plot(out$m[5L, ]); abline(h = true_phi_2)
plot(mu); lines(out$m[1L, ], col = "blue")
plot(xi_1); lines(out$m[2L, ], col = "blue")
plot(y); lines(out$f[, 1L], col = "blue")



# level + seas + ar(2) ----------------------------------------------------

# Simulate the data
nobs <- 500
sd_y <- 1
sd_mu <- 0.02
sd_a <- 1e-6
sd_b <- 1e-6
true_phi_1 <- 0.8
true_phi_2 <- -0.5

y <- numeric(nobs)
mu <- numeric(nobs)
a <- numeric(nobs)
b <- numeric(nobs)
xi_1 <- numeric(nobs)
xi_2 <- numeric(nobs)
xi_1[1L] <- 0
mu[1L] <- 5
a[1L] <- 0.1
b[1L] <- 0.1

period <- 12L # annual
r <- 1L # number of harmonics
w <- 2 * pi * r / period

# Random errors
set.seed(1111)
nu <- rnorm(n = nobs, sd = sd_y)
omega_mu <- rnorm(n = nobs, sd = sd_mu)
omega_a <- rnorm(n = nobs, sd = sd_a)
omega_b <- rnorm(n = nobs, sd = sd_b)

# First observation
y[1L] <- mu[1L] + a[1L] + xi_1[1L] + nu[1L]

for (t in seq_len(nobs)[-1]) {
  # AR2
  xi_1[t] <- true_phi_1 * xi_1[t - 1] + true_phi_2 * xi_2[t - 1] + nu[t]
  xi_2[t] <- xi_1[t - 1]
  # Level
  mu[t] <- mu[t - 1] + omega_mu[t]
  # Seasonality
  a[t] <- a[t - 1] * cos(w) + b[t - 1] * sin(w) + omega_a[t]
  b[t] <- b[t - 1] * cos(w) - a[t - 1] * sin(w) + omega_b[t]
  # Response
  y[t] <- mu[t] + a[t] + xi_1[t]
}

mod <- dlm(
  polynomial = list(order = 1, discount_factor = 0.95),
  seasonal = list(type = "fourier", period = 12, harmonics = 1,
                  discount_factor = 0.99),
  autoregressive = list(order = 2, discount_factor = 0.998)
)
m0 <- matrix(c(y[1L], rep(0, 6)), ncol = 1)
C0 <- diag(x = c(4, 2, 2, rep(1, 4)), nrow = 7L)
out <- forward_filter_dlm(y = y, F = mod$FF, G = mod$GG, D = mod$D,
                          m = m0,
                          C = C0,
                          n = 1, s = 1, df_variance = 1,
                          ar_order = 2L,
                          n_parms = 3L)

plot(out$m[6L, ]); abline(h = true_phi_1)
plot(out$m[7L, ]); abline(h = true_phi_2)
plot(mu); lines(out$m[1L, ], col = "blue")
plot(xi_1); lines(out$m[4L, ], col = "blue")
plot(y); lines(out$f[, 1L], col = "blue")

