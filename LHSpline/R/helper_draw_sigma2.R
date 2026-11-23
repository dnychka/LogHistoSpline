draw_sigma2 <- function(a, a0, b0) {
  shape_sigma2 <- a0 + length(a) / 2
  rate_sigma2  <- b0 + 0.5 * sum(a^2)
  
  1 / rgamma(1, shape = shape_sigma2, rate = rate_sigma2)
}