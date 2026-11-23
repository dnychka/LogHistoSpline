# Poisson log-likelihood for binned counts
#
# Model:
#   y_j ~ Poisson(mu_j)
#   eta_j = X_j %*% beta + gamma_j
#   f_j = exp(eta_j) / sum_k binWidth_k * exp(eta_k)
#   mu_j = n_eff * binWidth_j * f_j
#
# Interpretation:
#   f_j is the normalized estimated density height in bin j.
#   n_eff * binWidth_j is the exposure/offset for bin j.
#   The normalization ensures sum_j binWidth_j * f_j = 1.
#
# The constant term sum(log(y_j!)) is omitted because it does not depend
# on beta or gamma and therefore does

loglik_counts <- function(beta, gamma, X, y, n_eff, binWidth) {
  eta <- as.vector(X %*% beta + gamma)
  
  if (length(y) != length(eta)) {
    stop("y must have length nrow(X).")
  }
  if (length(gamma) != length(eta)) {
    stop("gamma must have length nrow(X).")
  }
  if (length(binWidth) == 1L) {
    binWidth <- rep(binWidth, length(eta))
  }
  if (length(binWidth) != length(eta)) {
    stop("binWidth must be scalar or have length nrow(X).")
  }
  if (n_eff <= 0 || any(binWidth <= 0)) {
    stop("n_eff and binWidth must be positive.")
  }
  
  eta <- pmin(eta, 700)
  mu_log <- log(n_eff * binWidth) + eta
  mu <- exp(mu_log)

  sum(y * mu_log - mu)
}






