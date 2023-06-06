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
  
}