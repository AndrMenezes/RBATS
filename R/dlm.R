#' @title Bayesian Dynamic Linear Model
#'
#' @description Create object of class \code{dlm}.
#'
#' @param polynomial list. Order of polynomial model.
#' @param seasonal list. Components to specify seasonal components.
#' @param regressor list.
#' @param df_variance numeric. Discount factor for observation variance. Use a beta-gamma random walk.
#' @param variance_law list. Variance law \code{type} and \code{power} parameter.
#' @param monitor list.
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
                df_variance = 1,
                variance_law = list(type = "identity", power = 2),
                monitor = list(
                  execute = FALSE,
                  verbose = TRUE,
                  start_time = 10L,
                  bilateral = FALSE,
                  bf_threshold = 0.135,
                  location_shift = 4,
                  scale_shift = 1,
                  discount_factors = list(
                    polynomial = 0.20,
                    seasonal = 0.80,
                    regressors = 0.80))) {

  if (polynomial[["order"]] < 1L) {
    polynomial[["order"]] <- 1L
    polynomial[["discount_factor"]] <- 0.95
    warning("It is not possible to specify negative or zero polynomial order. Using the level model (polynomial_order = 1)")
  }

  # Create model structure
  mod <- .polynomial_model(order = polynomial[["order"]],
                           discount_factors = polynomial[["discount_factor"]])

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
    }
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
  }

  colnames(mod[["GG"]]) <- rownames(mod[["GG"]]) <- comp_names

  # Discount factor and law for the observational variance
  mod[["df_variance"]] <- df_variance
  mod[["variance_law"]] <- variance_law

  # Additional information
  mod[["polynomial_order"]] <- polynomial[["order"]]
  seasonal[["type"]] <- seas_type
  mod[["seasonal"]] <- seasonal
  mod[["parameters_names"]] <- comp_names
  mod[["loglik"]] <- 0
  mod[["time"]] <- 0L
  mod[["call"]] <- match.call()

  # Include monitor information if it is available
  if (monitor[["execute"]]) {
    # Polynomial
    df_exception <- 1/monitor[["discount_factors"]][["polynomial"]]
    if (length(df_exception) != polynomial[["order"]]) {
      df_exception <- c(df_exception[1L], diag(mod[["D"]])[2:polynomial[["order"]]])
    }
    # Seasonal
    if (seas_type != "none") {
      df_seas <- monitor[["discount_factors"]][["seasonal"]]
      df_seas <- if (is.null(df_seas)) diag(mod_seas[["D"]]) else rep(1/df_seas, nrow(mod_seas[["FF"]]))
      df_exception <- c(df_exception, df_seas)
    }
    # Regressor
    if (!is.null(regressor[["xreg"]])) {
      df_reg <- 1/monitor[["discount_factors"]][["regressor"]]
      df_reg <- if (length(df_reg) != ncol(regressor[["xreg"]])) rep(df_reg, ncol(regressor[["xreg"]])) else df_reg
      df_reg <- if (is.null(df_reg)) 1/diag(mod_reg[["D"]]) else df_reg
      df_exception <- c(df_exception, df_reg)
    }
    # Exception discount factor matrix
    D_exception <- mod[["D"]]
    diag(D_exception) <- df_exception
    mod[["monitor"]] <- list(
      verbose = monitor[["verbose"]],
      start_time = monitor[["start_time"]],
      bilateral = monitor[["bilateral"]],
      bf_threshold = monitor[["bf_threshold"]],
      location_shift = monitor[["location_shift"]],
      scale_shift = monitor[["scale_shift"]],
      D = D_exception)
  }


  structure(mod, class = "dlm")
}
