setupDensityEvaluation<- function( gHat, binMidpoints,
                  threshold=c(TRUE, FALSE)
                  ){
  N<- length( binMidpoints)
  binWidth<- binMidpoints[2]- binMidpoints[1]
#xr<- xlim 
# spline interval may be different than the data range   
  xr<- range( binMidpoints)
  binWidth<- diff( binMidpoints)
  # for future use 
  
  # binBrk<- c(
  #   (binMidpoints[1] - (binMidpoints[2] - binMidpoints[1])/2),
  #   (binMidpoints[2:N] + binMidpoints[1:(N-1)])/2,
  #   (binMidpoints[N] + (binMidpoints[N] - binMidpoints[N-1])/2)
  #         )
  
  # use a linear discretization to evaluate log density
  yEnd<- c(gHat[1], gHat[N])
  yD<- c( (gHat[2]-gHat[1])/binWidth[1],
          (gHat[N]-gHat[N-1])/binWidth[N-1]
  )
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
  yGrid<- splintLinear( binMidpoints, gHat, xGrid)
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
  cdf<- splintLinear( xGrid, cdfTemp, binMidpoints )
  obj<-list(  
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
threshold= threshold
          )
      return(obj)
}

