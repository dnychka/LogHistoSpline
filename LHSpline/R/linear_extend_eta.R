linear_extend_eta <- function(eta_draws, z_grid, z_grid_out, zlim_obs) {
  if (nrow(eta_draws) != length(z_grid)) {
    stop("nrow(eta_draws) must equal length(z_grid).")
  }
  
  if (length(z_grid) < 2) {
    stop("z_grid must contain at least two points.")
  }
  
  if (any(diff(z_grid) <= 0) || any(diff(z_grid_out) <= 0)) {
    stop("z_grid and z_grid_out must be strictly increasing.")
  }
  
  eta_ext <- matrix(NA_real_, nrow = length(z_grid_out), ncol = ncol(eta_draws))
  
  for (j in seq_len(ncol(eta_draws))) {
    eta <- eta_draws[, j]
    
    slope_left <- (eta[2] - eta[1]) / (z_grid[2] - z_grid[1])
    slope_right <- (eta[length(eta)] - eta[length(eta) - 1]) /
      (z_grid[length(z_grid)] - z_grid[length(z_grid) - 1])
    
    # Force right tail to decrease.
    slope_right <- min(slope_right, -1e-8)
    
    eta_ext[, j] <- approx(
      x = z_grid,
      y = eta,
      xout = z_grid_out,
      rule = 1
    )$y
    
    left <- z_grid_out < z_grid[1]
    if (any(left)) {
      eta_ext[left, j] <- eta[1] +
        slope_left * (z_grid_out[left] - z_grid[1])
    }
    
    right <- z_grid_out > z_grid[length(z_grid)]
    if (any(right)) {
      eta_ext[right, j] <- eta[length(eta)] +
        slope_right * (z_grid_out[right] - z_grid[length(z_grid)])
    }
  }
  
  eta_ext
}
