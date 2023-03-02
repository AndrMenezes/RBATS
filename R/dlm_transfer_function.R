# Transfer function + level
g <- function(theta, x) {
  mu_t <- theta[1L]
  E_t <- theta[2L]
  beta_t <- theta[3L]
  lambda_t <- theta[4L]
  c(mu_t, lambda_t * E_t + beta_t * x, beta_t, lambda_t)
}
GG <- function(theta, x) {
  mu_t <- theta[1L]
  E_t <- theta[2L]
  beta_t <- theta[3L]
  lambda_t <- theta[4L]
  matrix(data = c(1, 0, 0, 0,
                  0, lambda_t, x, E_t,
                  0, 0, 1, 0,
                  0, 0, 0, 1),
         ncol = 4, byrow = TRUE)
}
# Matrix of discount factors
discount_factors <- c(0.95, 0.98, 0.98, 0.998)
D <- diag(x = 1/discount_factors, nrow = 4, ncol = 4)
D[-which(D == diag(D))] <- 1


