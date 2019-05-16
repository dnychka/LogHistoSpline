dsplineDensity <- function(obj, x = NULL) {
  # copy linear coefficients for log density tail 
  # right tail:  aTail[1]*( x- xr[1]) + bTail[1]
  #  left tail:  aTail[2]*( x- xr[2]) + bTail[2]
    y <- exp( splint(obj$logPdf$x, obj$logPdf$y, x) )
    if( obj$threshold[1]){
      y[ x< xr[1]]<-0
    }
    if( obj$threshold[2]){
      y[ x > xr[2] ]<-0
    }
  return(y)
}

rsplineDensity<- function( n, obj){
  qsplineDensity( obj, prob=runif( n))
}

rsplineDensityNEW<- function( n, obj){
  #qsplineDensity( obj, prob=runif( n))
   ind<- runif( n) < .5
   prob<- .5*runif( n)
   y<- rep( NA, n)
   y[ind]<-
         qsplineDensity( obj, prob=prob[ind])
   y[!ind]<- 
          qsplineDensity( obj, prob=prob[!ind], probTail=TRUE)
  return( y)
}

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




plot.splineDensity<- function( obj, 
                               which=c( TRUE, TRUE, TRUE, TRUE),
                               truePdf=NULL,
                               xlim =NULL, ylim=NULL
){
  #  set.panel( 3,1)
  if( which[1]){
    plot( obj$hist, freq=FALSE )
    lines(obj$pdf, co="blue", lwd=2)
    if( !is.null( truePdf)){
      lines( truePdf, col="grey20", lty=2, lwd=2)
    }
    
  }
  if( which[2]){
    plot(obj$logPdf, col="blue", type="l", lwd=2,
         xlim=xlim, ylim =ylim,
         ylab="log density function", xlab=" x ")
    xline( range(obj$Y), col="red")
    if( !is.null( truePdf)){
      lines( truePdf$x,log(truePdf$y), col="grey20", lty=2, lwd=2 )
    }
  }
  ## QQ plot 
  if( which[3]){
    N<- length( obj$Y)
    sortY<- sort(obj$Y)
    pY<- psplineDensity( obj, sortY)
    pGrid<-  ( (1:N) -.5)/N
    QY<- qsplineDensity(obj, pGrid)
    plot( QY, sortY, pch=16,
          xlab= "Fitted Quantiles",
          ylab="Sorted Sample", col="grey30")
    abline( 0,1,col="red", lwd=2) 
  }
  
  if( which[4]& (nrow(obj$info) > 1 ) ){
    whichCV<- obj$whichCV
    CVlabels<- c("-lnLike", "CVGreenSilverman",
                 "CVIntensity", "CVUnweighted")
    plot( obj$info[,1], obj$info[,whichCV], type="b",
          lwd=1, col="grey", pch=16, log="x",
          main= CVlabels[whichCV-2], 
          ylab="CV score",
          xlab="log lambda"
          )
    xline( obj$lambdaEst)
  }
}