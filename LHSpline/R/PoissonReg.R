

PoissonReg<- function( X, y, Nit=50, tol=1e-6){
  
n<- nrow( X)

# starting value for lambda  is the mean
nuOld<- rep( log(mean( y)), n)
  
for( k in 1:Nit){
  #I 
  lambdaTmp<- exp( nuOld) 
  W<- c(lambdaTmp)
  # pseudo data 
  z<- nuOld  + (1/W)*( y- lambdaTmp)
  # weighted least squares 
  tempObj<- lm(z~X-1,weights=W)
  nuNew <- tempObj$fitted.values
  # test for convergence
  test<- ( mean( abs(nuNew- nuOld)))
  if( test<= tol){
    break
  }
  nuOld<- nuNew
}
return( 
  list( fit=tempObj, 
        test=test, 
        iter = k,
        betaHat = tempObj$coefficients,
        nuHat=nuNew,
        fitted.values = exp(nuNew),
        coefficients= tempObj$coefficients)
)
}


  
  