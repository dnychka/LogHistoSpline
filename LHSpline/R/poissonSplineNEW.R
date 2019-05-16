
poissonSpline<- function( x,y, lambdaGrid,  gStart= NULL,
                     nIteration=50, tol=1e-6,
                     epsilon= 0,
                     verbose=FALSE,
                     whichCV=4){
   require(mgcv)
# this is used as 
# x midpoints of histogram bins
# y the bin counts
# gHat -- cubic spline is the log intensity
   N<- length( x)
# useful decompositions related to prior and penalty matrix 
   qrT<- qr( cbind( 1, x))
   Q<- qr.qy(qrT, diag(N) )
   Q2<- Q[,3: N]
   Q1<- Q[,1: 2]
   K<- Rad.cov( x, x , m=2)
   Sigma<- t(Q2)%*%K%*%Q2 
   SigmaInv<- solve( Sigma)
#      
   N<- sum(y)
   if( is.null(gStart)){
# 4 knot spline as the starting estimate   
        dataLocal<- data.frame( x=x, counts= y)
        fitMgcv<- gam(  counts~ s(x, k= 4) + x
                        ,family=poisson, data=dataLocal)
        gStart<- predict(fitMgcv)
  }
    nGrid<- length( lambdaGrid)
    gHat<- NULL
    info<- NULL
    z<- NULL
# normalize gStart to counts    
    gStart<- gStart - log( sum( exp( gStart))) + log( N)
    hStart<- exp( gStart)
    for( k in 1:nGrid ){
      output<- poissonSplineLambda( 
                     x = x,
                     y = y,
                gStart = gStart,
                lambda = lambdaGrid[k],
              SigmaInv = SigmaInv, 
               epsilon = epsilon, 
                     Q = Q,
            nIteration = nIteration,
                   tol =  1e-8,
               verbose = verbose
                    
      )
      info<- rbind( info, c( output$lambda, 
                             output$lnLike,
                             output$CVlnLike,
                             output$CVap,
                             output$CVap2,
                             output$CVap3,
                             output$lnLikeLinear,
                             output$nit)
                    )
      gHat<- cbind( gHat, output$gHat)
      z<- cbind(z, output$z )
# now set next starting values from this one.      
      gStart<- output$gHat
    }
    # basic idea is to define  the log density based on an interpolating 
    # cubic and linear extrapolation beyond the range on histogram bins 
    # ( see xr)
    # simple maximization over a grid of lambda's 
    
    if( length(lambdaGrid) > 1){
      if( (whichCV <3)| (whichCV >6)){
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
      gStartFinal<-gHat[, which.min(approxL2CV )] 
      # refit at this lambda
      fit<- poissonSplineLambda( 
        x = x,
        y = y,
        gStart = gStartFinal,
        lambda = exp(tempGrid[I]),
        SigmaInv = SigmaInv, 
        epsilon = epsilon, 
        Q = Q,
        nIteration = nIteration,
        tol =  1e-8,
        verbose = FALSE
      )
    }
    else{
      fit<- output
    }
# the fit list has the single estimate from CV or just one lambda.    
 return( list(x=x,y=y, 
         gHat= fit$gHat,
         lambdaGrid= lambdaGrid,
         lambdaEst= fit$lambda,
         info= info, 
         z=z,    # this is just for testing
         fit=fit,
         fitted.values= fit$fitted.values,
         gStart=gStart
         )
       )
}
##################################################
############ fitting engine
#################################################
poissonSplineLambda<- function( x,y, 
                                gStart,
                               lambda, 
                               SigmaInv, 
                               epsilon, 
                               Q,
                        nIteration = 10, 
                               tol = 1e-6,
                           verbose = FALSE,
                            gFloor = -15,
                          gCeiling =  10
                               ) {
  
  gOLD<- gStart 
# set up some matrices that depend on lambda 
  N<- length(y)
  Omega<- matrix( 0, N,N)
  Omega[3:N,3:N ]<-  lambda*SigmaInv
# NOTE:  a diffuse prior on the linear terms:
#  MN( 0, (1/epsilon)I)
  Omega[1,1]<- epsilon
  Omega[2,2]<- epsilon
  Gamma <- Q%*%Omega%*%t(Q)
# begin iteration loop  
  for( iteration in 1:nIteration){
    # notation follows fisher scoring Green and Silveman page 100
    # truncate gOLD to reasonable values
    gOLD<- ifelse( gOLD > gCeiling, gCeiling, gOLD)
    gOLD<- ifelse( gOLD <   gFloor,   gFloor, gOLD)
    bPrime<- exp(gOLD)
    # NOTE for  Poisson exponential  bPrime == bDoublePrime
    z<- (y - bPrime )/ bPrime + gOLD
    # here bPrime == w 
    weights<- bPrime
    # R is Green&Silverman "K"
    if( verbose){
    print( lambda)
    print( stats( cbind(gOLD, weights, z) ) )
    }
    gNEW<- splint( x, z, xgrid=x, wt= weights, lambda= lambda/N)
    #
    #other ways of finding gNEW -- should be identical to above.
    #
    #gNew<- Tps( x, z,weights=weights, lambda= lambda,
    #            scale.type="unscaled", GCV=FALSE)$fitted.values
    #gNEW<- solve( (diag(weights) + Gamma),  c(weights*z) )
    
    # convergence checks.     
    testTol<- sqrt(mean( (gNEW-gOLD)^2)/mean(gOLD^2)  )
    if( verbose){
      cat( "iteration: ", iteration, "  ", 
           sqrt(mean( (gNEW-gOLD)^2)), testTol,
           fill=TRUE) 
    }
    if( testTol< tol) break
    gOLD<- gNEW
  }
  # convert to beta parameters
  gHat<- gNEW
  hHat<-  exp( gNEW)
  beta<- Q%*%gHat
  
  #approx cross validated estimate of g_i omitting y_i 
  bPrime<- exp(gHat)
  weights<- bPrime
  W<- diag(weights )
  # NOTE for  Poisson exponential  bPrime == bDoublePrime
  z<- (y - bPrime)/ bPrime + gHat
  diagAMatrix<- diag( solve( W  + Gamma )%*%W )
  CVres<- (z - gHat)/( 1- diagAMatrix )
  # approximating the estimate of g_i having left out that 
  # observation, y_i
  gHatI<- z - CVres 
  # minus log Poission likelihood of CV estimates.     
  CVlnLike<- -1*sum( y*gHatI - exp(gHatI) ) 
  # CV suggested by Green and Silverman 
  CVap<- mean( bPrime*(z- gHat)^2 / (1-diagAMatrix)^2 )
  # Approx log Likelihood from linear approximation
  CVap2<- mean( bPrime*( y - exp(gHatI) )^2 )
  # leaving out the weights 
  CVap3<- mean(( y - exp(gHatI) )^2 )
  # approximate Gaussian process likelihood 
  # this does not seem to work!
  # lnDetA<- 2*sum( log(diag(
  #                  chol( W + Gamma ))))
  #lnLikeLinear<- -(1/2)* t(z)%*%solve(W + Gamma)%*%z + (1/2)* lnDetA
  # take minus this 
  # lnLikeLinear<- -1* lnLikeLinear
   lnLikeLinear<- NA
# approx posterior density integrating out beta
# Not sure this is working or correct  
#  hessian<-   t(Q) %*% diag( hHat) %*% Q  
# lnDetR<- 2*sum( log(diag(
#   chol( Gamma ))))
#  lnDetH <-    2*sum( log(diag(
#    chol( hessian + Gamma ))))             
#  scriptL<- sum( gHat*y) - sum( hHat ) -
#    (t(gHat)%*%Gamma%*%gHat)/2 +
#    + (1/2)*lnDetR 
#  lnLike<- scriptL - (1/2)*lnDetH 
  lnLike<- NA
  
  if( any( gHat==gFloor)){
    warning("gHat set to the lower bound at some values")
  }
  
  if( any( gHat==gCeiling)){
    warning("gHat set to the upper bound at some values")
  }
  
  return( list( nit= iteration, 
                gHat= gHat,
                fitted.values= exp( gHat),
                lambda=lambda,
                beta = beta,
                lnLike=lnLike,
                CVlnLike= CVlnLike,
                CVap= CVap,
                CVap2= CVap2,
                CVap3= CVap3,
                epsilon= epsilon,
                lnLikeLinear= lnLikeLinear,
                z= z,
                gStart=gStart)
  )
  
}  
