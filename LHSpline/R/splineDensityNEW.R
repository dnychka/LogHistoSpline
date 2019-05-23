splineDensity<- function( Y, xlim = NULL,
                  threshold=c(FALSE, FALSE), N=150,
                  lambda=NA, 
                  nGrid=10,
                  verbose=FALSE,debug=FALSE, epsilon=0,
                  gStart=NULL, expandRange = NULL,
                  whichCV=4
                  ){
  if( is.na( lambda[1])){
    lambda<- 10**seq( -2,2,,nGrid)
  }
  require(fields)
  if( is.null( xlim)){
    xlim<- range( Y, na.rm=TRUE)
  }
  if( !is.null( expandRange)){
    xlim[2]<- xlim[2]*expandRange
    N<- round( N*expandRange)
  }
  ind<- Y< xlim[1] | Y > xlim[2] | is.na(Y)
  Y<- Y[!ind]
#
  binWidth<- (xlim[2]-xlim[1])/ (N-1)
  brk<- seq( xlim[1] - binWidth/2, xlim[2]+binWidth/2,
             length.out= (N+1) )

  binMidpoints<- brk[1:N] + (brk[2]-brk[1])/2
  histObject<- hist( Y, breaks=brk, plot=FALSE)
  fitMgcv<- NULL
  if(debug){
# this is just included for comparision of old method that uses mgcv
#    
  require(mgcv)
   dataLocal<- data.frame( x=binMidpoints, counts= histObject$counts)
    fitMgcv<- gam(  counts~ s(x, k= N-3) + binMidpoints
     ,family=poisson, data=dataLocal, ...)
#   print( summary( fitMgcv))
  }
# heavy lifting happens here. NotepoissonSpline
# is a wrapper for for loop over individual lambda fits.  
   fit<- poissonSpline(x=binMidpoints, y=histObject$counts,
                       lambdaGrid= lambda, verbose=verbose,
                       gStart=gStart, whichCV=whichCV)
# info is a  5 column matrix lambda, three criteria for lambda 
# and number of interations for convergence. 
# info[,4] is the usual least squares CV applied to linear 
# problem at convergence. 
   info<- fit$info
# basic idea is to define  the log density based on an interpolating 
# cubic and linear extrapolation beyond the range on histogram bins 
# ( see xr)
# this is either CV minimum from grid search or the result of single
# lambda fit.    
   lambdaEst<- fit$lambdaEst
   gHat<- fit$gHat
#xr<- xlim 
# spline interval may be different than the data range   
  xr<- range( binMidpoints)
  dx<- xr[2] - xr[1]
# spline and derivatives at the endpoints
  yEnd<- splint( binMidpoints, gHat,
               c( xr[1], xr[2]) )
  yD<-  splint( binMidpoints, gHat,
               c( xr[1], xr[2]),
                      derivative=1 )
  I1<-I3<- 0
  if( !threshold[1]){
# integral before xr[1] assuming exponential tails. 
  if( yD[1]< 0){
    warning("left derviative at end point is negative")
  }
    I1<- exp( yEnd[1])/ yD[1]
  }
  if( !threshold[2]){
  # integral beyond xr[2] assuming exponential tails. 
  I3<-  -exp( yEnd[2])/ yD[2]
  if( yD[2] > 0){
    warning("right derviative at end point is positive")
  }
  }

# integrate the spline part
# interval is xr[1], xr[2]
# I2 is its approximate integral
# huge grid for simple integral  
  NG<- 40*N
  xGrid<-  seq( xr[1], xr[2],, NG)
  dx0<- binMidpoints[2]- binMidpoints[1]
  dx<- xGrid[2] - xGrid[1]
  yGrid<- splint( binMidpoints, gHat, xGrid)
# simple integration rule over the interval
  # cdf in the interval xr evaluated  at bin midpoints of histogram bins
  # normalized so that it is identical with the integral value I2 found
  # above
  cdfTemp<- c(0,cumsum( .5*(exp( yGrid[-1]) + exp( yGrid[-NG]) ))*dx )
  I2<-cdfTemp[NG]
  cdfTemp<- cdfTemp/I2
  nConstant<-   I1 + I2 + I3
  #print( sum(exp(gHat)*dx0) )
  #print( sum(exp(yGrid)*dx) )
  #print( I2)
  #print( nConstant)
  probLeft<-    I1 / nConstant
  probMiddle<-  I2 / nConstant 
  probRight<-   I3 / nConstant
  cdfTemp<-  (probLeft + probMiddle*cdfTemp)
#print( c( I1,I2, I3))
  gHat<- gHat - log(nConstant)
# get linear coefficients for log density tail 
# right tail:  aTail[1]*( x- xr[1]) + bTail[1]
#  left tail:  aTail[2]*( x- xr[2]) + bTail[2]
  bTail<- yEnd - log(nConstant)
  aTail<- yD
# values at  bin midpoints should be stable because xGrid is much larger 
# than number of bins
  cdf<- splint( xGrid, cdfTemp, binMidpoints )
  obj<-list(      Y = Y,
         fitPoisson = fit, 
               hist = histObject,
             logPdf = list( x = binMidpoints, y = gHat     ),
                pdf = list( x = binMidpoints, y = exp(gHat)), 
                cdf = list( x = binMidpoints, y = cdf),
# discretise the cdf estimate so that the quantile computation can be accurate.
            cdfRoot = list( x = xGrid, y= cdfTemp),
                 xr = xr,
               probLeft= probLeft,
                probRight= probRight,
                probMiddle= probMiddle,
                  I = c(I1,I2,I3),
              aTail = aTail,
              bTail = bTail,
          threshold = threshold,
             lambda = lambda,
          lambdaEst = fit$lambdaEst,
            fitMgcv = fitMgcv,
               info = info,
         xlim=xlim,
         whichCV= whichCV
          )

  class(obj)<- "splineDensity"
      return(obj)
}

