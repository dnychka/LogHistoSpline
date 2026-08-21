
qsplineDensity<- function(obj, prob, probTail=FALSE) {
  # copy linear coefficients for log density tail 
  # right tail:  aTail[1]*( x- xr[1]) + bTail[1]
  #  left tail:  aTail[2]*( x- xr[2]) + bTail[2]
  
  a <- obj$aTail
  b <- obj$bTail
  xr<- obj$xr
  # probabilities of right tail interval, middle from histospline and left. 
  I <- obj$I
  probRange<- range(obj$cdfRoot$y )
  if( probTail){
  ind1 <- prob > 1-probRange[1]
  ind3 <- prob < 1-probRange[2]
  }
  else{
    ind1 <- prob < probRange[1]
    ind3 <- prob > probRange[2]
  }
  ind2 <- !(ind3 | ind1)
   Q <- rep(NA, length(prob))
  if (!obj$threshold[1]) {
    if( probTail){
      Q[ind1] <-  xr[1] + (log((1-prob[ind1])*a[1])- b[1])/ a[1] 
    }
    else{
      Q[ind1] <-  xr[1] + (log((prob[ind1])*a[1])- b[1])/ a[1]
    }
    #  exp( a[1]*(x[ind1]- xr[1]) + b[1] ) / a[1]
  }
  else{
    Q[ind1] <- xr[1]
  }
  if (!obj$threshold[2]) {
    if( probTail){
      Q[ind3] <-  xr[2] + (log( -prob[ind3]*a[2] )- b[2])/ a[2]
    }
    else{
      Q[ind3] <-  xr[2] + (log( -(1-prob[ind3])*a[2] )- b[2])/ a[2]
    }
  }
  else{
    Q[ind3] <- 1.0
  }
   if( probTail){ 
  Q[ind2] <- approx(obj$cdfRoot$y, obj$cdfRoot$x, 1-prob[ind2])$y}
   else{
# piecewise linear interpolation  
  Q[ind2] <- approx(obj$cdfRoot$y, obj$cdfRoot$x, prob[ind2])$y
   }
  return(Q)
}


