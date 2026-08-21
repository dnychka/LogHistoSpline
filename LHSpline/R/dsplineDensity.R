dsplineDensity <- function(obj, x = NULL) {
  # copy linear coefficients for log density tail 
  # right tail:  aTail[1]*( x- xr[1]) + bTail[1]
  #  left tail:  aTail[2]*( x- xr[2]) + bTail[2]
    y <- exp( splint(obj$logPdf$x, obj$logPdf$y, x) )
    xr<- obj$xr
    if( obj$threshold[1]){
      y[ x< xr[1]]<-0
    }
    if( obj$threshold[2]){
      y[ x > xr[2] ]<-0
    }
  return(y)
}
