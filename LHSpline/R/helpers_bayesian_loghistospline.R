# Internal helper: map low-rank coefficients a into the full gamma vector.
build_gamma_from_a <- function(a, Q2, vecs_k, vals_k, inside_obs = NULL) {
  u <- as.vector(vecs_k %*% (sqrt(vals_k) * a))
  gamma <- as.vector(Q2 %*% u)
  
  # Optional boundary mask. If inside_obs is 1 inside the observed z range
  # and 0 outside, this turns off the GP deviation in the extrapolation
  # region so tail behavior is controlled by the polynomial basis.
  if (!is.null(inside_obs)) {
    gamma <- gamma * inside_obs
  }
  
  list(u = u, gamma = gamma)
}

# Internal helper: log-likelihood as a function of low-rank coefficients a.
# Used inside the gamma/a elliptical slice sampler.
loglik_gamma_lowrank <- function(a, beta, X, y, Q2, vecs_k, vals_k,
                                 n_eff, binWidth, inside_obs = NULL) {
  tmp <- build_gamma_from_a(a, Q2, vecs_k, vals_k, inside_obs)
  loglik_counts(beta, tmp$gamma, X, y, n_eff, binWidth)
}

# Internal helper: normalize eta draws into log-density draws on the z grid.
normalize_log_density_draws <- function(eta_draws, binWidth) {
  if (length(binWidth) == 1L) {
    binWidth <- rep(binWidth, nrow(eta_draws))
  }
  stopifnot(length(binWidth) == nrow(eta_draws))
  
  apply(eta_draws, 2, function(eta) {
    m <- max(eta)
    logZ <- m + log(sum(binWidth * exp(eta - m)))
    eta - logZ
  })
}

