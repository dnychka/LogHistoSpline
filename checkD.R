library( logPdf)
library( fields)
library( splines)
library( extRemes)

pR<- 1/(c(5,10,25,50,100)*365)
shapePar<- 3
scalePar<- 1/2
N<- 10*365
GPDpart<-  function(u,tau){ifelse(u<tau, 1, u )}
returnTrue<- qgamma( pR, shape=shapePar, scale=scalePar, lower.tail=FALSE)
tailProb<- .05
nReps<- 200
set.seed( 124)
  Y<- rgamma(N, shape=shapePar, scale=scalePar)
  threshold<- quantile( Y, 1- tailProb)
  GPFit<- fevd( Y, threshold=threshold, type="GP",
                method="MLE")
  look<- splineDensity( Y, xlim=c(0,8), threshold=c( TRUE, FALSE))
  xg<- look$logPdf$x
  scaleE<- var( Y)/ mean( Y)
  shapeE<- var(Y)/scaleE^2  
  
# some plots  
  
plot(xg, log(dgamma( xg,shape=shapePar, scale=scalePar) ) ,
       col="red", lwd=2, type="l", xlim=c( 2,7))

#lines( xg, log(dgamma( xg,shape=shapeE, scale=scaleE) ), lty=2, lwd=3)
#xline( max( Y),lwd=4)

print(qgamma( pR, shape=shapePar, scale=scalePar, lower.tail=FALSE))
print(qgamma( pR, shape=shapeE, scale=scaleE, lower.tail=FALSE))
print(return.level(GPFit, c(5,10,25,50,100)))
dGPhat<- devd(xg, scale=GPFit$results$par[1], shape= GPFit$results$par[2],
              threshold= GPFit$threshold, type="GP")*tailProb
lines( xg+GPFit$threshold, log(dGPhat ), col="green4", lwd=4)

lines( look$logPdf, col="green", lwd=2)

tau<- 2

look2<- testD(Y,xlim=c( 0, 8),  tau=tau, threshold= c( TRUE, FALSE))
lines( look2$logPdf, col="blue", lwd=4, lty=2)
xline( c(tau,tau+ look2$delta), col="blue")


#plot( look2$return.level, log="x")

#rt<- max( look2$cdf$x[look2$cdf$y>= pR[1]] )
  
  
