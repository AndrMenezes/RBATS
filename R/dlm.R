#' @title Bayesian Dynamic Linear Model
#'
#' @description Create object of class \code{dlm}.
#'
#' @param polynomial_order integer. Order of polynomial model.
#' @param seasonal list. Components to specify seasonal components.
#' @param xreg matrix. Regressors.
#' @param discount_factors list. Discount factors for each model components.
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
dlm <- function(polynomial_order = 1L,
                seasonal = list(type = c("none", "free", "fourier"),
                                period = NULL, harmonics = NULL),
                xreg = NULL,
                discount_factors = list(polynomial = 0.95, seasonal = 0.98,
                                        regressors = 0.99),
                df_variance = 1,
                variance_law = list(type = "identity", power = 2)) {

  if (polynomial_order < 1L) {
    polynomial_order <- 1L
    discount_factors[["polynomial"]] <- 0.95
    warning("It is not possible to specify negative or zero polynomial order. Using the level model (polynomial_order = 1)")
  }

  # Create model structure
  mod <- .polynomial_model(order = polynomial_order,
                           discount_factors = discount_factors[["polynomial"]])

  # Create model structure for seasonality
  seas_type <- match.arg(seasonal[["type"]], c("none", "free", "fourier"))
  if (seas_type != "none") {
    if (seas_type == "free") {
      mod_seas <- .seasonal_free_model(period = seasonal[["period"]],
                                       discount_factors = discount_factors[["seasonal"]])
    }
    if (seas_type == "fourier") {
      mod_seas <- .seasonal_fourier_model(period = seasonal[["period"]],
                                          harmonics = seasonal[["harmonics"]],
                                          discount_factors = discount_factors[["seasonal"]])
    }
    # Superposition of models
    mod[["FF"]] <- rbind(mod[["FF"]], mod_seas[["FF"]])
    mod[["GG"]] <- .bdiag(mod[["GG"]], mod_seas[["GG"]])
    mod[["D"]] <- .bdiag_one(mod[["D"]], mod_seas[["D"]])
  }

  comp_names <- rownames(mod[["FF"]])
  # Create model structure for regressors
  if (!is.null(xreg)) {
    mod_reg <- .regression_model(X = xreg,
                                 discount_factors = discount_factors[["regressors"]])
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
  mod[["polynomial_order"]] <- polynomial_order
  seasonal[["type"]] <- seas_type
  mod[["seasonal"]] <- seasonal
  mod[["parameters_names"]] <- comp_names
  mod[["loglik"]] <- 0
  mod[["time"]] <- 0L
  mod[["call"]] <- match.call()

  structure(mod, class = "dlm")
}
