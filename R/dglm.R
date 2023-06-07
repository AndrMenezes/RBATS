update.poisson_dglm <- function(y, F, G, a, R, D) {

  # Prior for linear predictor

  # Conjugate prior for mu_t (VB step)

  # Predictive distribution (parameters)

  # Posterior for mu_t (conjugate)

  # Update the posterior moments of the linear predictor

  # Linear Bayes to update the posterior moments of the state

  # Discount factor for the W_t

}
forward_filter.poisson_dglm <- function(y, F, G, a, R, D) {

}

poisson_dglm <- function(polynomial = list(order = 1L, discount_factor = 0.95),
                         seasonal = list(type = c("none", "free", "fourier"),
                                         period = NULL, harmonics = NULL,
                                         discount_factor = 0.98),
                         regressor = list(xreg = NULL, discount_factor = 0.99)) {

  if (polynomial[["order"]] < 1L) {
    polynomial[["order"]] <- 1L
    polynomial[["discount_factor"]] <- 0.95
    warning("It is not possible to specify negative or zero polynomial order. Using the level model (polynomial_order = 1)")
  }

  # Create model structure
  mod <- .polynomial_model(order = polynomial[["order"]],
                           discount_factors = polynomial[["discount_factor"]])
  mod[["i_polynomial"]] <- 1:polynomial[["order"]]

  # Create model structure for seasonality
  seas_type <- match.arg(seasonal[["type"]], c("none", "free", "fourier"))
  if (seas_type != "none") {

    if (seas_type == "free") {
      mod_seas <- .seasonal_free_model(period = seasonal[["period"]],
                                       discount_factors = seasonal[["discount_factor"]])
    }
    if (seas_type == "fourier") {
      mod_seas <- .seasonal_fourier_model(period = seasonal[["period"]],
                                          harmonics = seasonal[["harmonics"]],
                                          discount_factors = seasonal[["discount_factor"]])
      # Seasonal effect from Fourier components
      mod[["L"]] <- .fourier_to_seasonal(
        period = seasonal[["period"]],
        number_harmonics = length(seasonal[["harmonics"]]),
        FF = mod_seas[["FF"]], GG = mod_seas[["GG"]]
      )
    }

    # Index of the seasonal components
    mod[["i_seasonal"]] <- (nrow(mod[["FF"]]) + 1):(nrow(mod[["FF"]]) + nrow(mod_seas[["FF"]]))

    # Superposition of models
    mod[["FF"]] <- rbind(mod[["FF"]], mod_seas[["FF"]])
    mod[["GG"]] <- .bdiag(mod[["GG"]], mod_seas[["GG"]])
    mod[["D"]] <- .bdiag_one(mod[["D"]], mod_seas[["D"]])
  }

  comp_names <- rownames(mod[["FF"]])
  # Create model structure for regressors
  if (!is.null(regressor[["xreg"]])) {
    mod_reg <- .regression_model(X = regressor[["xreg"]],
                                 discount_factors = regressor[["discount_factor"]])
    mod[["xreg"]] <- mod_reg[["xreg"]]
    mod[["GG"]] <- .bdiag(mod[["GG"]], mod_reg[["GG"]])
    mod[["D"]] <- .bdiag_one(mod[["D"]], mod_reg[["D"]])
    comp_names <- c(comp_names, colnames(mod_reg[["GG"]]))
    # Index of the regressor components
    mod[["i_regressor"]] <- (nrow(mod[["FF"]]) + 1):(nrow(mod[["FF"]]) + ncol(mod_reg[["xreg"]]))
  }

  colnames(mod[["GG"]]) <- rownames(mod[["GG"]]) <- comp_names

  # Additional information
  mod[["polynomial_order"]] <- polynomial[["order"]]
  seasonal[["type"]] <- seas_type
  mod[["seasonal"]] <- seasonal
  mod[["parameters_names"]] <- comp_names
  mod[["loglik"]] <- 0
  mod[["time"]] <- 0L
  mod[["call"]] <- match.call()


  structure(mod, class = "poisson_dglm")

}
