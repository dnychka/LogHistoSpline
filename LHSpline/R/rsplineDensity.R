rsplineDensity<- function( n, obj){
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
