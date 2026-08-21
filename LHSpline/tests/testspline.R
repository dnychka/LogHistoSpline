library( LHSpline)
library( fields)

set.seed(123)
y<- rnorm(50)
fit<- splineDensity(y, N=10)

xr<- fit$xr
xg<- seq(xr[1], xr[2],,1e4)
yg<-dsplineDensity( fit, x=xg)
dx<- xg[2]-xg[1]
test2<- .5*sum( yg[-1]+ yg[-1e4])*dx

NP<- 1e4
xg<- seq( -8, xr[1], , NP)
yg<-dsplineDensity( fit, x=xg)
dx<- xg[2] - xg[1]
yg2<- exp(fit$bTail[1] + fit$aTail[1]*( xg- xr[1]))
test1<-  .5*sum(  yg[-1] + yg[-NP] ) *dx
test1b<- .5*sum( yg2[-1] + yg2[-NP]) *dx
print( exp(fit$bTail[1])/fit$aTail[1] )

yg<-psplineDensity( fit, x=xg)
test2b<- psplineDensity( fit, x= fit$xr)
test.for.zero(test2, diff( test2b), tol=5e-6,
               tag="middle interval" )

test3<- splint( fit$logPdf$x, fit$logPdf$y,c( xr[1], xr[1]-1))
test.for.zero( -diff(test3), fit$aTail[1])
test3b<- splint( fit$logPdf$x, fit$logPdf$y,c( xr[1]-1.5))
test4<- fit$aTail[1]*(-1.5) + fit$bTail[1]
test.for.zero( test3b, test4)





