rm(list = ls())
library(ggplot2)

# Simulate the data
nobs <- 100
sd_nu <- 10
true_phi <- 0.8

y <- numeric(nobs)
y[1L] <- 0

set.seed(66)
for (t in seq_len(nobs)[-1]) {
  y[t] <- true_phi * y[t - 1] + rnorm(n = 1, sd = sd_nu)
}
x11()
plot(y, main = expression(y[t]))
graphics.off()


# Model components
g <- function(theta) {
  xi_t <- theta[1L]
  phi_t <- theta[2L]
  c(phi_t * xi_t, phi_t)
}
GG <- function(theta) {
  xi_t <- theta[1L]
  phi_t <- theta[2L]
  matrix(data = c(phi_t, xi_t,
                  0, 1),
         ncol = 2, byrow = TRUE)
}
FF <- matrix(c(1, 0), ncol = 1)

# Matrix of discount factors
discount_factors <- c(0.98, 0.998)
D <- diag(x = 1/discount_factors, nrow = 2, ncol = 2)
D[-which(D == diag(D))] <- 1
df_variance <- 1

# Prior for the state parameters
a0 <- c(10, 0)
R0 <- diag(x = c(1, 1), nrow = 2)
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
parms_names <- c("xi", "phi")

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
  a <- g(theta = m)
  GG_m <- GG(theta = m)
  P <- tcrossprod(GG_m %*% C, GG_m)

  # Discount information
  R <- D * P
  cat(
    formatC(c(R[1, 1], s), digits = 4, width = 4),
    "\n")
  n <- df_variance * n

  # Set names
  colnames(R) <- rownames(R) <- colnames(C) <- parms_names
  names(a) <- names(m) <- parms_names

  # Storage the estimates
  list_posterior[[t]] <- list(m = m, C = C, n = n, s = s)
  list_predictive[[t]] <- list(f = f, q = q)
  list_prior[[t + 1]] <- list(a = a, R = R, W = R - P, n = n, s = s)

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
  # geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey69",
  #             alpha = 0.8) +
  geom_line()

x11()
ggplot(data_posterior[data_posterior$parameter == "lambda1", ],
       aes(x = t, y = mean)) +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey69",
              alpha = 0.8) +
  geom_line() +
  geom_hline(yintercept = true_lambda_1, col = "red")

x11()
ggplot(data_posterior[data_posterior$parameter == "lambda2", ],
       aes(x = t, y = mean)) +
  geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), fill = "grey69",
              alpha = 0.8) +
  geom_line() +
  geom_hline(yintercept = true_lambda_2, col = "red")




