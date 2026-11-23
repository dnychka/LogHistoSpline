### set.seed(123); z<- rnorm( 2e3); look<- createCounts(z, 50)
createCounts<- function( z, N, expandRange=1.0, zlim=NULL){
  # range expanded to add some zero counts
  # 
  zlim_obs <- range(z)
  
  if( is.null(zlim)){
    zlim<- c( zlim_obs[1], zlim_obs[2]*expandRange )
  }
  
  # use simpler binning that should not matter for large N ...
  # also preserves zero  at left. 
  #
  brk<-  seq(zlim[1] ,
             zlim[2] ,
             length.out = N + 1)
  binWidth <- (brk[2] - brk[1]) 
  histObject <- hist(z, breaks = brk, plot = FALSE)
  y <- histObject$counts
  # these are the bin midpoints
  
  z_grid <- (brk[1:N] + brk[ 2: (N+1)])/2
  
  #x_grid <- exp(z_grid)
  n_eff <- sum(y) 
  
  # Step 2: polynomial basis and orthogonal complement
  # z_mean <- mean(z_grid)
  # z_sd <- sd(z_grid)
  # z_center <- (z_grid - z_mean) / z_sd
   zR<- range( z_grid)
  z_center<- (z_grid - zR[1])/( zR[2]- zR[1])
  ######### X matrix 
  delta<- z_center[2] -z_center[1]
  #X <- cbind(1, z_center, log(z_grid))
  # X <- cbind(1, z_center, 
  #            #z_center^2, 
  #            log( z_center + delta/2 )
  #            )
  return( 
    list( 
      binWidth= binWidth,
      z_grid = z_grid,
      y=y,
      z_center=z_center, 
          #X=X,
      brk= brk)
    )
}