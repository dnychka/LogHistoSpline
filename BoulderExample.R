library(extRemes)
library(dataWorkshop)
data( BoulderDaily)
setwd("/Users/nychka/Home/Src/logPDF")
y<- BoulderDaily$precip
y<- y*2.54
ind<- y> 0
yr<- BoulderDaily$time
thresh<- 2


par( mar=c( 5,5,1,1))
plot( yr, y,  type="h", ylab="(cm)", xlab="Years", ylim =c(0,15))

ind<- y>= thresh
points( yr[ind], y[ind], col="green4", cex=.5, pch=16)
yline( thresh, col="red4", lwd=.5)
ind2<- which.max( y)
text( yr[ind2], 13, label="23.1 cm ",lty=2, lwd=2, adj=1, col="green4")
arrows(yr[ind2],0,yr[ind2], 15, col="green4", lwd=3, angle=20, length=.1)
#dev.copy2pdf(file="/Users/nychka/Home/Current_talks/AGUExtremes/pix/Boulder1.pdf")

ind<- !is.na(y)
BoulderFit<- fevd( y[ind], threshold= thresh, type="GP")
ind<- !is.na(y) & y < 20
BoulderFit<- fevd( y[ind], threshold= thresh, type="GP")
hold<- ci.fevd( BoulderFit, return.period= 25)


par( mar=c( 5,5,1,1))
fields.style()
ind2<- (y >= thresh)& (y< 20)
brk<- seq( 1e-4, 15,,40)
hist( y[ind2&ind], breaks=brk, xlim=c(0,15), prob=TRUE, col="green4", border="green2", xlab="Daily Precip (cm)", ylab="Density", main="", ylim=c(0,.8) )
#xline( c( 8.14, 9.55, 10.9), lty=c(2,1,2), col="grey30", 
#lwd=3)
xline( c( 8.15, 9.8, 11.5), lty=c(2,1,2), col="grey",  lwd=3)
xline( 2, col="red4", lwd=2)

xgrid<- seq( 0, 13,,300)
ydenGP<- devd( xgrid, threshold= thresh, type="GP", scale = BoulderFit$results$par[1],       shape = BoulderFit$results$par[2]) 

lines( xgrid+thresh, ydenGP, col="black", lwd=3)
xline( 2, col="red4", lwd=2)

#dev.copy2pdf(file="/Users/nychka/Home/Current_talks/AGUExtremes/pix/Boulder2.pdf")

###################

library(dataWorkshop)
data( BoulderDaily)

y<- BoulderDaily$precip
y<- y*2.54
thresh2<- .1
ind3<- !is.na(y) & y>= thresh2
yr<- BoulderDaily$time
#

Sfit<- splineDensity(log(y[ind3]), lambda=.1)
xgrid1<- seq( thresh, 13,,300)
delta<- (xgrid1[2] - xgrid1[1])
xraw<- exp(Sfit$predicted[,1])
yden<- splint( xraw, Sfit$predicted[,2], xgrid)
yden<- exp( yden)/xgrid1
yden<- yden/ sum( yden*delta)

Sfit2<- splineDensity(log(y[ind3]), sp=1)
xgrid1<- seq( thresh, 13,,300)
delta<- (xgrid1[2] - xgrid1[1])
xraw<- exp(Sfit2$predicted[,1])
yden2<- splint( xraw, Sfit2$predicted[,2], xgrid)
yden2<- exp( yden2)/xgrid1
yden2<- yden2/ sum( yden2*delta)

##################################
brk<- seq( thresh2*(1- 1e-3), 25,,50)
hist( y[ind3], breaks=brk, xlim=c(0,15), prob=TRUE, col="green4", border="green2", xlab="Daily Precip (cm)", ylab="Density", main="", ylim=c(0,.8) )
lines( xgrid1, yden, lwd=3)
#dev.copy2pdf(file="/Users/nychka/Home/Current_talks/AGUExtremes/pix/Boulder3.pdf")

#################################
par( mar=c( 5,5,1,5))
fields.style()
matplot( xgrid1, cbind( yden ,yden2), lwd=2, col="grey40", type="l", ylab="density", xlab="(cm)")
lines(xgrid+thresh, (ydenGP), col="black", lwd=2)

usrPar<- par()$usr
usrPar[3:4] <- c( -9, 2)
par( usr=usrPar)


matlines( xgrid1, cbind(log( yden), log(yden2)), lwd=2, col="green2", type="l")
lines(xgrid+thresh, log(ydenGP), col="green4", lwd=2)
axis( side=4)
mtext( side=4, text="log Density", line = 3, cex=1.5, col="green4")

#dev.copy2pdf(file="/Users/nychka/Home/Current_talks/AGUExtremes/pix/Boulder4a.pdf")

