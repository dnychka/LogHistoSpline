
poissonSpline<- function( x,y, lambdaGrid,  gStart= NULL,
                     nIteration=50, tol=1e-6,
                     verbose=FALSE,
                     whichCV=2){
  
 if( any(diff(x)< 0 )){stop("x needs to be in sorted order")}
 if( any(diff(x) == 0 )){stop("cannot  handle replicated x locattions")}
# this is used as 
# x midpoints of histogram bins
# y the bin counts
# gHat -- cubic spline is the log intensity
   N<- length( x)

#      
   N<- sum(y)
   if( is.null(gStart)){
   require(mgcv)
# 4 knot spline as the starting estimate   
        dataLocal<- data.frame( x=x, counts= y)
        fitMgcv<- gam(  counts~ s(x, k= 4) + x
                        ,family=poisson, data=dataLocal)
        gStart<- predict(fitMgcv)
    }
    nGrid<- length( lambdaGrid)
    info<- NULL
# normalize gStart to counts    
    gStart<- gStart - log( sum( exp( gStart))) + log( N)
    for( k in 1:nGrid ){
      gPrevious<- gStart
      output<- poissonSplineLambda( 
                     x = x,
                     y = y,
                gStart = gPrevious,
                lambda = lambdaGrid[k],
            nIteration = nIteration,
                   tol =  1e-8,
               verbose = verbose
      )
      info<- rbind( info, c( 
                             output$CVlnLike,
                             output$CVap,
                             output$CVap2,
                             output$CVap3,
                             output$nit,
                             output$trA)
                    )
# now set next starting values from this one.      
      gPrevious<- output$gHat
    }
    dimnames( info)<- list( NULL, c( "CVlnLike",
                                   "CVap", "CVap2" ,  "CVap3",
                                            "nit", "trA")
    )
    # basic idea is to define  the log density based on an interpolating 
    # cubic and linear extrapolation beyond the range on histogram bins 
    # ( see xr)
    # simple maximization over a grid of lambda's 
    if( length(lambdaGrid) > 1){
      if( (whichCV < 1)| (whichCV >4)){
        stop("bad whichCV value")
      }
      approxL2CV<- info[,whichCV]
      lr<- range( log( lambdaGrid))
      tempGrid<- seq( lr[1], lr[2], ,1e3)
      tempGridY<- splint( log(lambdaGrid),approxL2CV,
                         tempGrid )
      I<- which.min(
        splint(log( lambdaGrid), approxL2CV, tempGrid ) 
      )
      if( verbose){
        cat("lambda at maximum:", 
            exp(tempGrid[I]), I,  fill=TRUE )
      }
      # refit at this lambda
      fit<- poissonSplineLambda( 
        x = x,
        y = y,
        gStart = gStart,
        lambda = exp(tempGrid[I]),
        nIteration = nIteration,
        tol =  1e-8,
        verbose = verbose
      )
    }
    else{
      fit<- output
    }
# the fit list has the single estimate from CV or just one lambda. 
    ix<- order( x)
    predicted<- list( x= x[ix], y=fit$gHat[ix])
    outObject<- list(x=x,y=y, 
                     predicted= predicted,
                     fitted.values= fit$fitted.values,
                     gHat= fit$gHat,
                     lambdaGrid= lambdaGrid,
                     lambdaEst= fit$lambda,
                     whichCV= whichCV,
                     info= info, 
                     z=fit$z,   
                     weights=fit$weights,
                     gStart=gStart,
                     fit=fit
    )
#     
    class(outObject)<- "PoissonSpline"
 return( outObject )
}
