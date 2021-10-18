#' @title Smoothing moments of state parameters
#'
#' @description Method to smoothing moments of state parameters and the
#' one-step ahead predictive distributions for objects of class  \code{dlm}.
#'
#' @param object an object of \code{dlm}.
#' @param time integer. The time index.
#'
#' @author André F. B. Menezes
#'
#' @references
#' Prado, R.; West, M. Time Series Modeling, Computation, and Inference. CRC Press, 2010.

#' @rdname retrospective_moments
#' @export
retrospective_moments <- function(object, time){
  UseMethod("retrospective_moments", object)
}

#' @rdname retrospective_moments
#' @export
retrospective_moments.dlm <- function(object, time) {

  # t = time and T = end_time
  GG <- object[["GG"]]
  FF <- object[["FF"]]
  if (!is.null(object[["xreg"]]))
    FF <- rbind(FF, t(object[["xreg"]][time, , drop = FALSE]))
  sT <- object[["sT"]]
  st <- object[["predictive"]][[time]][["s"]]

  # a_T(t - T + 1) and R_T(t - T + 1)
  a_T_t_1 <- object[["smooth"]][[time + 1]][["a_T"]]
  R_T_t_1 <- object[["smooth"]][[time + 1]][["R_T"]]

  # Posterior at t and prior at t + 1
  m_t <- object[["posterior"]][[time]][["m"]]
  C_t <- object[["posterior"]][[time]][["C"]]
  a_t_1 <- object[["prior"]][[time + 1]][["a"]]
  R_t_1 <- object[["prior"]][[time + 1]][["R"]]
  B_t <- tcrossprod(C_t, GG) %*% solve(R_t_1)

  # Retrospective the prior and predictive
  # a_T(t - T) and R_T(t - T)
  a_T <- drop(m_t - B_t %*% (a_t_1 - a_T_t_1))
  R_T <- C_t - tcrossprod(B_t %*% (R_t_1 - R_T_t_1), B_t)

  # Change in the scale to reflect the update error variance estimate (t-Student)
  R_T <- sT / st * R_T

  # f_t(-k) e q_t(-k)
  f_T <- drop(crossprod(FF, a_T))
  q_T <- drop(crossprod(FF, R_T %*% FF))

  list(a_T = a_T, R_T = R_T, f_T = f_T, q_T = q_T)

}

