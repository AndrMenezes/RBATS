rm(list = ls())
library(ggplot2)

# Simulate the data
nobs <- 80
sd_y <- 10
sd_mu <- 0.1
sd_E <- 0.5
sd_beta <- 0.1
true_lambda <- 0.85

E <- numeric(nobs)
beta <- numeric(nobs)
mu <- numeric(nobs)
y <- numeric(nobs)
E[1L] <- 0
beta[1L] <- 10
mu[1L] <- 2

set.seed(66)
y[1L] <- mu[1L] + E[1L] + rnorm(n = 1, sd = sd_y)
x <- numeric(nobs) # rpois(nobs, lambda = 1)
x[c(5, 20, 30, 40, 60)] = c(4, 2.5, 10, 7.5, 5)

for (t in seq_len(nobs)[-1]) {
  beta[t] <- beta[t - 1] + rnorm(n = 1, sd = sd_beta)
  E[t] <- true_lambda * E[t - 1] + beta[t] * x[t] + rnorm(n = 1, sd = sd_E)
  mu[t] <- mu[t - 1] + rnorm(n = 1, sd = sd_mu)
  y[t] <- E[t] + rnorm(n = 1, sd = sd_y)
}
all.equal(length(y), length(E), length(beta), length(mu))
x11()
par(mfrow = c(4, 1))
plot(y, main = expression(y[t]))
plot(mu, main = expression(mu[t]))
plot(E, main = expression(E[t]))
plot(beta, main = expression(beta[t]))
graphics.off()


# Model components
g <- function(theta, x) {
  mu_t <- theta[1L]
  E_t <- theta[2L]
  lambda_t <- theta[3L]
  beta_t <- theta[4L]
  c(mu_t, lambda_t * E_t + beta_t * x, lambda_t, beta_t)
}
GG <- function(theta, x) {
  mu_t <- theta[1L]
  E_t <- theta[2L]
  lambda_t <- theta[3L]
  matrix(data = c(1, 0, 0, 0,
                  0, lambda_t, E_t, x,
                  0, 0, 1, 0,
                  0, 0, 0, 1),
         ncol = 4, byrow = TRUE)
}
FF <- matrix(c(1, 1, 0, 0), ncol = 1)

# Matrix of discount factors
discount_factors <- c(0.95, 0.995, 0.995, 0.995)
D <- diag(x = 1/discount_factors, nrow = 4, ncol = 4)
D[-which(D == diag(D))] <- 1
df_variance <- 1

# Prior for the state parameters
a0 <- c(mu[1L], 0, 0, 0)
R0 <- diag(x = c(20, 100, 0.025, 10), nrow = 4)
# Prior for the observational variance
n0 <- 1
s0 <- 1

# Define list to storage the parameters
list_prior <- list(list())
list_prior[[1L]][["a"]] <- a0
list_prior[[1L]][["R"]] <- R0
list_prior[[1L]][["W"]] <- D * R0
list_prior[[1L]][["n"]] <- n0
list_prior[[1L]][["s"]] <- s0
list_posterior <- list()
list_predictive <- list()
loglik <- 0
parms_names <- c("mu", "E", "lambda", "beta")

# Sequential learning
for (t in seq_len(nobs)) {

  # Predictive t
  f <- drop(crossprod(FF, a0))
  q <- s0 + drop(crossprod(FF, R0 %*% FF))

  # Posterior at t
  e <- y[t] - f

  # Adaptive coefficient vector
  A <- (R0 %*% FF) / q

  # Volatility estimate ratio
  r <- (n0 + e^2 / q) / (n0 + 1)

  # Kalman filter update
  n <- n0 + 1
  s <- r * s0
  m <- drop(a0 + A %*% e)
  C <- r * (R0 - q * tcrossprod(A))

  # Log-predictive likelihood
  lpl <- log(1/sqrt(q) * dt((y - f) / sqrt(q), n0))

  # Prior t + 1
  a <- g(theta = m, x = x[t])
  GG_m <- GG(theta = m, x = x[t])
  P <- tcrossprod(GG_m %*% C, GG_m)

  # Discount information
  R <- D * P
  n <- df_variance * n

  # Set names
  colnames(R) <- rownames(R) <- colnames(C) <- parms_names
  names(a) <- names(m) <- parms_names

  # Storage the estimates
  list_posterior[[t]] <- list(m = m, C = C, n = n, s = s)
  list_predictive[[t]] <- list(f = f, q = q)
  list_prior[[t]] <- list(a = a, R = R, W = R - P, n = n, s = s)

  # Update the prior
  a0 <- a
  R0 <- R
  s0 <- s
  n0 <- n

  # Update the log-predictive likelihood
  loglik <- loglik + lpl

}


# Inspect the results
data_predictive <- data.frame(t = seq_along(list_predictive), y = y)
data_predictive$mean <- sapply(list_predictive, function(z) unname(z[["f"]]))
data_predictive$variance <- sapply(list_predictive, function(z) unname(z[["q"]]))
data_predictive$ci_lower <- data_predictive$mean - 1.96 * sqrt(data_predictive$variance)
data_predictive$ci_upper <- data_predictive$mean + 1.96 * sqrt(data_predictive$variance)
ggplot(data_predictive[-c(1:6), ], aes(x = t, y = y)) +
  geom_point(size = 2) +
  geom_line(aes(y = mean), col = "blue") +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey69",
              alpha = 0.8)


time_index <- seq_len(nobs)
lt_parms <- lapply(parms_names, function(parm) {
  tmp <- lapply(list_posterior, function(z) {
    diag_C <- diag(z[["C"]])
    c("mean" = unname(z[["m"]][which(names(z[["m"]]) == parm)]),
      "variance" = unname(diag_C[names(diag_C) == parm]))
  })
  values <- do.call(rbind, tmp)
  out <- data.frame(t = time_index)
  out$parameter <- parm
  out$mean <- values[, "mean"]
  out$variance <- values[, "variance"]
  out
})
data_posterior <- do.call(rbind, lt_parms)
data_posterior$ci_lower <- data_posterior$mean - 1.96 * sqrt(data_posterior$variance)
data_posterior$ci_upper <- data_posterior$mean + 1.96 * sqrt(data_posterior$variance)

ggplot(data_posterior[data_posterior$t > 6, ], aes(x = t, y = mean)) +
  facet_wrap(~parameter, scales = "free_y") +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey69",
              alpha = 0.8) +
  geom_line()


ggplot(data_posterior[data_posterior$parameter == "lambda", ],
       aes(x = t, y = mean)) +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey69",
              alpha = 0.8) +
  geom_line() +
  geom_hline(yintercept = true_lambda, col = "red")




