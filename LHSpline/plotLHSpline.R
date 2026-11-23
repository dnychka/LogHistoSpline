plotLHSpline<- function( obj, nclass=50){
  hist( obj$z, nclass=nclass, prob=TRUE)
  lines( obj$z_grid, exp( obj$gHat), col="blue")
}