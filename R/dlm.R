#' @title Bayesian Dynamic Linear Model
#'
#' @description Create object of class \code{dlm}.
#'
#' @param polynomial list. Order of polynomial model.
#' @param seasonal list. Components to specify seasonal model.
#' @param regressor list. Components to specify a dynamic regression model.
#' @param autoregressive list. Components to specify a dynamic autoregressive model.
#' @param cycle list. Components to specify a dynamic cycle model.
#' @param transfer_function list. Components to specify a dynamic transfer function model.
#' @param df_variance numeric. Discount factor for observation variance. Use a beta-gamma random walk.
#' @param variance_law list. Variance law \code{type} and \code{power} parameter.
#' The variance law \code{type} are \code{identity}, \code{poisson}, \code{binomial},
#' and \code{power}. The variance law \code{power} should be numeric \eqn{p \geq 1}.
#'
#' @author André F. B. Menezes
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.

#' @rdname dlm
#' @export
dlm <- function(polynomial = list(order = 1L, discount_factor = 0.95),
                seasonal = list(type = c("none", "free", "fourier"),
                                period = NULL, harmonics = NULL,
                                discount_factor = 0.98),
                regressor = list(xreg = NULL, discount_factor = 0.99),
                autoregressive = list(order = NULL, discount_factor = 0.998),
                cycle = list(frequency = NULL, rho = NULL, discount_factor = 0.998),
                transfer_function = list(order = NULL, xreg = NULL, discount_factor = 0.998),
                df_variance = 1,
                variance_law = list(type = "identity", power = 1)) {

  if (missing(polynomial) & missing(seasonal) & missing(regressor) &
      missing(autoregressive) & missing(cycle) & missing(transfer_function)) {
    polynomial <- list(order = 1L, discount_factor = 0.95)
    warning("The default model is the local level with discount factor 0.95")
  }

  if (!missing(polynomial)) {
    # Polynomial component
    mod <- .polynomial_model(order = polynomial[["order"]],
                             discount_factors = polynomial[["discount_factor"]])
    mod[["i_polynomial"]] <- 1:polynomial[["order"]]
    comp_names <- rownames(mod[["FF"]])
  }

  # Seasonality component
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
        FF = mod_seas[["FF"]],
        GG = mod_seas[["GG"]]
      )
    }
    # Index of the seasonal components
    mod[["i_seasonal"]] <- (nrow(mod[["FF"]]) + 1):(nrow(mod[["FF"]]) + nrow(mod_seas[["FF"]]))
    # Superposition of models
    mod[["FF"]] <- rbind(mod[["FF"]], mod_seas[["FF"]])
    mod[["GG"]] <- .bdiag(mod[["GG"]], mod_seas[["GG"]])
    mod[["D"]] <- .bdiag_one(mod[["D"]], mod_seas[["D"]])
    # Component names
    comp_names <- c(comp_names, rownames(mod_seas[["FF"]]))
  }

  # Create regression component
  if (!is.null(regressor[["xreg"]])) {
    mod_reg <- .regression_model(X = regressor[["xreg"]],
                                 discount_factors = regressor[["discount_factor"]])
    if (missing(polynomial) & missing(seasonal)) {
      mod <- mod_reg
      comp_names <- rownames(mod[["FF"]])
      mod[["i_regressor"]] <- 1:nrow(mod[["FF"]])
    } else {
      # Index of the regressor components
      mod[["i_regressor"]] <- (nrow(mod[["FF"]]) + 1):(nrow(mod[["FF"]]) + ncol(mod_reg[["xreg"]]))
      # Superposition
      mod[["xreg"]] <- mod_reg[["xreg"]]
      mod[["FF"]] <- rbind(mod[["FF"]], mod_reg[["FF"]])
      mod[["GG"]] <- .bdiag(mod[["GG"]], mod_reg[["GG"]])
      mod[["D"]] <- .bdiag_one(mod[["D"]], mod_reg[["D"]])
      comp_names <- c(comp_names, colnames(mod_reg[["GG"]]))
    }
  }

  # Create cycle component
  if (!is.null(cycle[["freq"]])) {
    mod_cycle <- .cycle_model(freq = cycle[["freq"]], rho = cycle[["rho"]],
                              discount_factors = cycle[["discount_factor"]])
    if (missing(polynomial) & missing(autoregressive) & missing(regressor)) {
      mod <- mod_cycle
      comp_names <- rownames(mod[["FF"]])
      mod[["i_cycle"]] <- 1:nrow(mod[["FF"]])
    } else {
      # Index of the cycle components
      mod[["i_cycle"]] <- (nrow(mod[["FF"]]) + 1):(nrow(mod[["FF"]]) + nrow(mod_cycle[["FF"]]))
      # Superposition
      mod[["FF"]] <- rbind(mod[["FF"]], mod_cycle[["FF"]])
      mod[["GG"]] <- .bdiag(mod[["GG"]], mod_cycle[["GG"]])
      mod[["D"]] <- .bdiag_one(mod[["D"]], mod_cycle[["D"]])
      comp_names <- c(comp_names, colnames(mod_cycle[["GG"]]))
    }
  }

  # Create autoregressive components
  if (!is.null(autoregressive[["order"]])) {
    mod_autoreg <- .autoregressive_model(order = autoregressive[["order"]],
                                         discount_factors = autoregressive[["discount_factor"]])

    if (missing(polynomial) & missing(regressor)) {
      mod <- mod_autoreg
      comp_names <- rownames(mod[["FF"]])
      mod[["i_autoregressive"]] <- 1:nrow(mod[["FF"]])
    } else {
      # Index of the autoregressive components
      mod[["i_autoregressive"]] <- (nrow(mod[["FF"]]) + 1):(nrow(mod[["FF"]]) + nrow(mod_autoreg[["FF"]]))
      # Superposition
      mod[["FF"]] <- rbind(mod[["FF"]], mod_autoreg[["FF"]])
      mod[["GG"]] <- .bdiag(mod[["GG"]], mod_autoreg[["GG"]])
      mod[["D"]] <- .bdiag_one(mod[["D"]], mod_autoreg[["D"]])
      # Appending the component names
      comp_names <- c(comp_names, colnames(mod_autoreg[["GG"]]))
    }

  }

  # Create transfer function component
  if (!is.null(transfer_function[["order"]])) {
    mod_tf <- .transfer_function_model(order = transfer_function[["order"]],
                                       X = transfer_function[["xreg"]],
                                       discount_factors = transfer_function[["discount_factor"]])

    if (!exists("mod", environment(), inherits = FALSE)) {
      mod <- mod_tf
      comp_names <- rownames(mod[["FF"]])
      mod[["i_transfer_function"]] <- 1:nrow(mod[["FF"]])
    } else {
      # Index of the transfer function components
      mod[["i_transfer_function"]] <- (nrow(mod[["FF"]]) + 1):(nrow(mod[["FF"]]) + nrow(mod_tf[["FF"]]))

      # Superposition
      mod[["FF"]] <- rbind(mod[["FF"]], mod_tf[["FF"]])
      mod[["GG"]] <- .bdiag(mod[["GG"]], mod_tf[["GG"]])
      mod[["D"]] <- .bdiag_one(mod[["D"]], mod_tf[["D"]])
      # Appending the component names
      comp_names <- c(comp_names, colnames(mod_tf[["GG"]]))
    }
    # Covariate for the transfer function
    mod[["xreg_tf"]] <- mod_tf[["xreg"]]

  }

  colnames(mod[["GG"]]) <- rownames(mod[["GG"]]) <- comp_names

  # Discount factor and law for the observational variance
  mod[["df_variance"]] <- df_variance
  mod[["variance_law"]] <- variance_law

  # Additional information
  mod[["polynomial_order"]] <- if (missing(polynomial)) NULL else polynomial[["order"]]
  mod[["seasonal"]] <- if (missing(seasonal)) NULL else seasonal

  mod[["ar_order"]] <- if (missing(autoregressive)) 0 else autoregressive[["order"]]
  mod[["tf_order"]] <- if (missing(transfer_function)) 0 else transfer_function[["order"]]


  # Fill the index of each component
  mod[["i_cycle"]] <- if (is.null(mod[["i_cycle"]])) 0 else mod[["i_cycle"]]
  mod[["i_autoregressive"]] <- if (is.null(mod[["i_autoregressive"]])) 0 else mod[["i_autoregressive"]]
  mod[["i_transfer_function"]] <- if (is.null(mod[["i_transfer_function"]])) 0 else mod[["i_transfer_function"]]
  mod[["i_transfer_function"]] <- if (is.null(mod[["i_transfer_function"]])) 0 else mod[["i_transfer_function"]]


  mod[["parameters_names"]] <- comp_names
  mod[["n_parms"]] <- length(comp_names) - 2*mod[["ar_order"]] - 2*mod[["tf_order"]]
  mod[["loglik"]] <- 0
  mod[["time"]] <- 0L
  mod[["call"]] <- match.call()


  structure(mod, class = "dlm")
}
