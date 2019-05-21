psplineDensity<- function(obj, x, probTail=FALSE) {
  # copy linear coefficients for log density tail 
  # right tail:  aTail[1]*( x- xr[1]) + bTail[1]
  #  left tail:  aTail[2]*( x- xr[2]) + bTail[2]
  a <- obj$aTail
  b <- obj$bTail
  xr<- obj$xr
  
  # probabilities of right tail interval, middle from histospline and left. 
  I <- obj$I
    ind1 <- x < xr[1]
    ind3 <- x > xr[2]
    ind2 <- !(ind3 | ind1)
    y <- rep(NA, length(x))
    if (!obj$threshold[1]) {
      prob<- exp( a[1]*(x[ind1]- xr[1]) + b[1] ) / a[1]
      if( probTail){
        y[ind1]<-1-prob
      }else{
        y[ind1]<-prob
      }  
    }
    else{
      y[ind1] <- 0
    }
    if (!obj$threshold[2]) {
        prob<- -exp( a[2]*(x[ind3]- xr[2]) + b[2] ) / a[2]
      if( probTail){
        y[ind3]<-prob
      }else{
        y[ind3]<-1 - prob
      }
    }
    else{
      y[ind3] <- 1.0
    }
   #prob <- splint(obj$cdfRoot$x, obj$cdfRoot$y, x[ind2])
# piecewise linear interpolation    
    prob <- approx(obj$cdfRoot$x, obj$cdfRoot$y, x[ind2])$y
     if( probTail){
       y[ind2]<-1-prob
     }else{
       y[ind2]<-prob
    }
   
  return(y)
}

