#' @keywords internal
# one run of the random walk metropolis sampler for an alpha of a gamma distribution

# Random walk MH according to Clark 2003
pois_RW_once <- function(lambdas, alpha_old, beta, a0, b0, alpha_scale = 0.1){

  # Propose new alpha
  alpha_cand <- rgamma(1, shape = alpha_scale, rate = alpha_scale/alpha_old)

  # Posterior for alpha
  posterior_new <- alpha_llk(lambdas = lambdas, alpha = alpha_cand, beta = beta) + dgamma(alpha_cand, shape = a0, rate = b0, log = TRUE)
  posterior_old <- alpha_llk(lambdas = lambdas, alpha = alpha_old, beta = beta) + dgamma(alpha_old, shape = a0, rate = b0, log = TRUE)

  # print(posterior_old)

  # Acceptance ratio
  acc <- min(log(1),(posterior_new + dgamma(alpha_old, shape = alpha_scale, rate = alpha_scale/alpha_cand, log = TRUE)) - (posterior_old + dgamma(alpha_cand, shape = alpha_scale, rate = alpha_scale/alpha_old, log = TRUE)))

  # Draw from runif and return alpha
  if(acc < log(1)) {
    unif         <- log(runif(1))
  } else {
    unif         <- log(1)
  }
  if (unif <= acc) {
    return(list(alpha = alpha_cand, accept = 1))
  } else {
    return(list(alpha = alpha_old, accept = 0))
  }

}


#' @keywords internal
# one run of the random walk metropolis sampler for an alpha of a gamma distribution
alpha_llk <- function(lambdas, alpha, beta){

  n <- length(lambdas)
  return(n*log(exp(alpha*log(beta)-log(base::gamma(alpha)))) + (alpha-1)*sum(log(lambdas)))

}
