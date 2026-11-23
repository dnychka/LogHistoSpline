fitMLE<- function( y, z_grid,X, threshold= c( TRUE, FALSE) ){
  
  n_eff<- sum( y)
  N<- length(z_grid )
  binWidth<- z_grid[2]-z_grid[1]
  glmOffset<- rep(log(n_eff *binWidth), N)
  glmFit<- glm(
    y ~ X - 1 + offset(glmOffset),
    family = poisson()
  )
  etaHat<- predict( glmFit)
  
  predictList<- setupDensityEvaluation( etaHat, 
                                        z_grid,
                                        threshold=threshold
  )
  out<- c(list( etaHat=etaHat, glmFit=glmFit ), predictList)
  return(out)
}