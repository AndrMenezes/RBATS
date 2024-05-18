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
  FF <- matrix(rep(NA_real_, p))

  # D <- .diag_one(discount_factors)
  D <- diag(1/discount_factors, ncol = p, nrow = p)
  D[which(D != diag(D))] <- 1 / discount_factors[1L]

  colnames(D) <- rownames(D) <- colnames(GG) <- rownames(GG) <- rownames(FF) <- comp_names
  list(xreg = X, FF = FF, GG = GG, D = D)
}

.cycle_model <- function(freq, rho = NULL, discount_factors = 0.98) {

  p <- length(freq)
  n <- 2*p
  i_min <- seq.int(1, by = 2, length.out = p)
  i_max <- seq.int(2, by = 2, length.out = p)
  comp_names <- c()
  for (f in freq) {
    comp_names <- c(comp_names, paste0(
      "frequency_{", f, "}__", c("xi_1", "xi_2")) )
  }

  FF <- matrix(rep(c(1, 0), times = p))
  GG <- matrix(data = 0, nrow = n, ncol = n)
  for (j in seq_len(p)) {
    cs <- cos(freq[j])
    sn <- sin(freq[j])
    i_row <- i_min[j]:i_max[j]
    GG[i_row, i_row] <- matrix(c(cs, -sn, sn, cs), nrow = 2, ncol = 2)
  }
  if (!is.null(rho))
    GG <- rho * GG
  # Create the non-linear term for estimate the damping factor (rho)
  if (is.null(rho)) {
    FF <- rbind(FF, 0)
    GG <- cbind(GG, NA_real_)
    GG <- rbind(GG, c(0, 0, 1))
    n <- n + 1L
    comp_names <- c(comp_names, "rho")
  }

  # Discount factor matrix
  discount_factors_cov <- if (length(discount_factors) == n) 1 else discount_factors[1L]
  D <- diag(1/discount_factors, ncol = n, nrow = n)
  D[which(D != diag(D))] <- 1/discount_factors_cov

  colnames(GG) <- rownames(GG) <- comp_names
  colnames(D) <- rownames(D) <- comp_names
  rownames(FF) <- comp_names

  list(FF = FF, GG = GG, D = D)

}

.autoregressive_model <- function(order, discount_factors) {
  FF <- matrix(c(1, rep(0, 2 * order - 1)), ncol = 1L)

  GG <- diag(1, 2 * order)
  if (order > 1L) {
    # Pop up the second main diagonal of the matrix
    diag(GG[1:order, 1:order]) <- 0
    for (i in 1:(order - 1)) {
      GG[i + 1, i] <- 1
    }
  }
  # Fill with NA the first row to represent the posterior mean of past time
  GG[1L, ] <- NA_real_

  # Discount factor matrix
  discount_factors <- if (length(discount_factors) == 1L) rep(discount_factors, 2*order) else discount_factors
  D <- diag(1/discount_factors, ncol = 2 * order, nrow = 2 * order)
  D[which(D != diag(D))] <- 1

  # Components names
  comp_names <- c(paste0("xi_", 1:order), paste0("phi_", 1:order))

  colnames(D) <- rownames(D) <- colnames(GG) <- rownames(GG) <- rownames(FF) <- comp_names

  list(FF = FF, GG = GG, D = D)
}

.transfer_function_model <- function(order, X, discount_factors, lambda = NULL) {

  if ((order == 1) & !is.null(lambda)) {
    FF <- matrix(c(1, 0), ncol = 1L)
    GG <- diag(c(lambda, 1), nrow = 2)

    discount_factors <- if (length(discount_factors) == 1L) rep(discount_factors, 2) else discount_factors
    D <- diag(1/discount_factors, ncol = 2, nrow = 2)
    D[which(D != diag(D))] <- 1

    # Named the components
    comp_names <- c("E_1", "psi")

  }
  else {
    FF <- matrix(c(1, rep(0, 2 * order)), ncol = 1L)

    GG <- diag(1, 2 * order + 1)
    if (order > 1L) {
      # Pop up the second main diagonal of the matrix
      diag(GG[1:order, 1:order]) <- 0
      for (i in 1:(order - 1)) {
        GG[i + 1, i] <- 1
      }
    }
    # Fill with NA the first row to represent the posterior mean of past time
    GG[1L, ] <- NA_real_

    # Discount factor matrix
    discount_factors <- if (length(discount_factors) == 1L) rep(discount_factors, 2 * order + 1) else discount_factors
    D <- diag(1/discount_factors, ncol = 2 * order + 1, nrow = 2 * order + 1)
    D[which(D != diag(D))] <- 1

    # Named the components
    comp_names <- c(paste0("E_", 1:order), paste0("lambda_", 1:order), "psi")

  }

  colnames(D) <- rownames(D) <- colnames(GG) <- rownames(GG) <- rownames(FF) <- comp_names

  list(FF = FF, GG = GG, xreg = X, D = D)
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


