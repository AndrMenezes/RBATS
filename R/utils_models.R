.polynomial_model <- function(order = 1L, discount_factors = 0.9) {

  comp_names <- c("level", "growth", "growth_change", if (order > 3) paste0("poly_", 4:order))
  comp_names <- comp_names[1:order]

  FF <- matrix(c(1, rep(0, order - 1)), ncol = 1)
  GG <- diag(order)
  GG[row(GG) == col(GG) - 1] <- 1

  discount_factors_cov <- 1
  if (length(discount_factors) != order) {
    discount_factors <- rep(discount_factors[1L], order)
    discount_factors_cov <- discount_factors[1L]
  }
  D <- diag(1/discount_factors, ncol = order, nrow = order)
  D[which(D != diag(D))] <- 1/discount_factors_cov


  colnames(GG) <- rownames(GG) <- comp_names
  colnames(D) <- rownames(D) <- comp_names
  rownames(FF) <- comp_names

  list(FF = FF, GG = GG, D = D)
}

.seasonal_free_model <- function(period, discount_factors = 0.98) {

  p <- as.integer(period) - 1L
  comp_names <- paste0("period_", period, "__", "seas_", seq_len(p))

  FF <- matrix(c(1, rep(0, p - 1L)), ncol = 1)
  GG <- matrix(0, p, p)
  GG[row(GG) == col(GG) + 1] <- 1
  GG[1, ] <- -1

  discount_factors_cov <- if (length(discount_factors) == p) 1 else discount_factors[1L]
  D <- diag(1/discount_factors, ncol = p, nrow = p)
  D[which(D != diag(D))] <- 1/discount_factors_cov

  colnames(GG) <- rownames(GG) <- comp_names
  colnames(D) <- rownames(D) <- comp_names
  rownames(FF) <- comp_names

  list(FF = FF, GG = GG, D = D)
}

.seasonal_fourier_model <- function(period, harmonics, discount_factors = 0.98) {

  p <- length(harmonics)
  n <- 2*p
  i_min <- seq.int(1, by = 2, length.out = p)
  i_max <- seq.int(2, by = 2, length.out = p)
  comp_names <- paste0(
    "period_", period, "__",
    c("fourier_cos_", "fourier_sin_"), rep(harmonics, each = 2))

  FF <- matrix(rep(c(1, 0), times = p))
  GG <- matrix(data = 0, nrow = n, ncol = n)
  for (j in seq_len(p)) {
    cs <- cos(2 * pi * harmonics[j]/period)
    sn <- sin(2 * pi * harmonics[j]/period)
    i_row <- i_min[j]:i_max[j]
    GG[i_row, i_row] <- matrix(c(cs, -sn, sn, cs), nrow = 2, ncol = 2)
  }

  discount_factors_cov <- if (length(discount_factors) == n) 1 else discount_factors[1L]
  D <- diag(1/discount_factors, ncol = n, nrow = n)
  D[which(D != diag(D))] <- 1/discount_factors_cov

  colnames(GG) <- rownames(GG) <- comp_names
  colnames(D) <- rownames(D) <- comp_names
  rownames(FF) <- comp_names

  list(FF = FF, GG = GG, D = D)
}

.fourier_to_seasonal <- function(period, number_harmonics, FF, GG) {
  # Section 8.6.5 WH, pag 254
  L <- matrix(0, nrow = period, ncol = 2 * number_harmonics)
  L[1L, ] <- FF
  for (i in seq_len(period)[-1L]) {
    L[i, ] <- L[i - 1, ] %*% GG
  }
  L
}

.regression_model <- function(X, discount_factors = 0.98) {
  p <- ncol(X)
  if (length(discount_factors) != p) discount_factors <- rep(discount_factors[1L], p)
  comp_names <- if (is.null(colnames(X))) paste0("regressor_", seq_len(p)) else colnames(X)
  colnames(X) <- comp_names
  GG <- diag(nrow = p)

  # D <- .diag_one(discount_factors)
  D <- diag(1/discount_factors, ncol = p, nrow = p)
  D[which(D != diag(D))] <- 1/discount_factors[1L]

  colnames(D) <- rownames(D) <- colnames(GG) <- rownames(GG) <- comp_names
  list(xreg = X, GG = GG, D = D)
}

.bdiag_one <- function(...) {
  bX <- .bdiag(...)
  # Fill elements 0 to 1
  bX[which(bX == 0)] <- 1
  bX
}

# Function from dlm package. (https://github.com/cran/dlm/blob/master/R/DLM.R)
.bdiag <- function(...) {
  if (nargs() == 1L)
    x <- as.list(...)
  else x <- list(...)
  n <- length(x)
  if (n == 0)
    return(NULL)
  x <- lapply(x, function(y) if (length(y))
    as.matrix(y)
  else stop("Zero-length component in x"))
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1, ]
  cc <- d[2, ]
  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(4, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1, -1] <- rcum[-n]
  ind[2, ] <- rcum
  ind[3, -1] <- ccum[-n]
  ind[4, ] <- ccum
  imat <- array(1:(rsum * csum), c(rsum, csum))
  iuse <- apply(ind, 2, function(y, imat) imat[(y[1] + 1):y[2],
                                               (y[3] + 1):y[4]], imat = imat)
  iuse <- as.vector(unlist(iuse))
  out[iuse] <- unlist(x)
  out
}


