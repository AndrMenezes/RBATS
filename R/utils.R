.variance_law <- function(type, mu, p) {
  switch(type,
        "identity" = 1,
        "poisson" = mu,
        "binomial" = mu * (1 - mu),
        "power" = mu ^ p
  )
}

.diag_one <- function(x) {
  m <- diag(x, nrow = length(x), ncol = length(x))
  m[which(m != diag(m))] <- 1
  m
}

.matrix_power <- function(x, k) {
  x_k <- x
  i <- 2L
  while (i <= k) {
    x_k <- x_k %*% x
    i <- i + 1L
  }
  x_k
}

