testD<- function( x, xlim = range( x), tau=NULL, overlap=2.5,
                  nKnots=5, threshold=c(FALSE, FALSE), N=150, ... ){
  require(mgcv)
  require(fields)
  ind<- x< xlim[1] | x > xlim[2] | is.na(x)
  x<- x[!ind]
  brk<- seq( xlim[1], xlim[2],, (N+1) )
  mid<- brk[1:N] + (brk[2]-brk[1])/2
  histObject<- hist( x, breaks=brk, plot=FALSE)
  dataLocal<- data.frame( x=mid, counts= histObject$counts, tau=tau)
  LKNodes<- seq(0,tau,,nKnots)
  delta<- (LKNodes[2] - LKNodes[1]) *overlap
  D<- rdist( mid, LKNodes)
  X<- Wendland2.2( D/delta)
  fit<- gam(  counts~ X + mid , family=poisson, data=dataLocal,...)
  logf<- log(fit$fitted.values)
  #xr<- xlim
  xr<- range( mid)
  dx<- xr[2] - xr[1]
  
  
  # spline and derivatives at the endpoints
  yEnd<- splint( mid, logf,
                 c( xr[1], xr[2]))
  yD<-  splint( mid, logf,
                c( xr[1], xr[2]),
                derivative=1 )
  I1<-I3<- 0
  if( !threshold[1]){
    # integral before xr[1] assuming exponential tails.
    I1<- exp( yEnd[1])/ yD[1]
  }
  if( !threshold[2]){
    # integral beyond xr[2] assuming exponential tails. 
    I3<-  -exp( yEnd[2])/ yD[2]
  }
  
  # integrate the spline part
  # interval is xr[1], xr[2]
  # I2 is its approximate integral
  NG<- 10^4
  xGrid<-  seq( xr[1], xr[2],, NG)
  dx<- xGrid[2] - xGrid[1]
  yGrid<- splint( mid, logf, xGrid)
  # simple integration rule over the interval
  I2<- .5*sum( (exp( yGrid[-1]) + exp( yGrid[-N]) )*dx )
  print(c(I1,I2,I3))
  nConstant<- I1 + I2 + I3
  logf<- logf - log(nConstant)
  
  # get linear coefficients
  bTail<- yEnd - log(nConstant)
  aTail<- yD
  
  # cdf in the interval xr evaluated  at midpoints of histogram bins
  # normalized so that it is identical with the integral value I2 found
  # above
  cdfTemp<- c(0,cumsum( .5*(exp( yGrid[-1]) + exp( yGrid[-N]) )) )
  cdfTemp<- cdfTemp/ cdfTemp[NG]
  cdfTemp<-  (I1 + I2*cdfTemp)/ nConstant
  # values at  bin midpoints should be stable because xGrid is much larger 
  # than number of bins
  cdf<- splint( xGrid,cdfTemp, mid)
  
  obj<-list(      fit = fit, 
                  hist = histObject,
                  logPdf = list( x = mid, y = logf     ),
                  pdf = list( x = mid, y = exp(logf)), 
                  cdf = list( x = mid, y = cdf),
                  return.level= list( x=  1/((1-cdf)*365), y=mid),
                  xr = xr,
                  mid=mid,
                  logf= logf,
                  I = c(I1,I2,I3),
                  aTail = aTail,
                  bTail = bTail,
                  threshold = threshold,
                  LKNodes=LKNodes,
                  delta= delta)
  
  class(obj)<- "splineDensity"
  return(obj)
}

GPDPart<-  function(u,tau){ifelse(u<tau, 1, (u) )}