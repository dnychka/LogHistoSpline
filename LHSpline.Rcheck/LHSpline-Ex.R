pkgname <- "LHSpline"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('LHSpline')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("poissonSpline")
### * poissonSpline

flush(stderr()); flush(stdout())

### Name: poissonSpline
### Title: Estimates a Poisson cubic smoothing spline
### Aliases: poissonSpline
### Keywords: spatial

### ** Examples

set.seed(12)
x<- seq( 0,1,length.out= 250)
alphaTrue<-  exp( 25*x^3*(1-x) )
#alphaTrue<-  alphaTrue/ max( alphaTrue)
gTrue<- log(alphaTrue )
y<- rpois(length( alphaTrue), lambda=alphaTrue)


lambdaGrid<- 10**( seq( -3, -1,length.out=150))
fit1<- poissonSpline( x,y, lambdaGrid= lambdaGrid)

set.panel(3,1)
plot( x,y, pch=16, col="grey", cex=.5)
lines( x, alphaTrue, col="red", lty=1)
lines( x, fit1$fitted.values, col="blue", lty=1)
title("Counts and fitted curve")

plot( x,fit1$gHat, type="l" , col="blue")
lines(x , gTrue, col="red")
title("Estimated (blue) and true (red) log intensity")

plot( fit1$info[,1], fit1$info[,4] ,
type="b",
pch=16,cex=.5,
xlab="log lambda", ylab="CV function", log="x")
xline( (fit1$fit$lambda), col="blue")
title("CV criterion with minimum")





cleanEx()
nameEx("splineDensity")
### * splineDensity

flush(stderr()); flush(stdout())

### Name: splineDensity
### Title: Log density estimate as a cubic spline
### Aliases: splineDensity

### ** Examples

# an easy one
set.seed( 122)
Y<- rnorm( 2000)
fitY<- splineDensity( Y)
set.panel(2,2 )
plot( fitY)
set.panel()

# compact support
set.seed( 322)
Y <- rbeta( 2000, shape1= 5, shape2=.5)
fitY<- splineDensity( Y, xlim=c(0,1), threshold=c(TRUE,TRUE) )
set.panel(2,2 )
plot( fitY)
set.panel()
xgrid<- seq( 0,1,,200)
hist( Y, prob=TRUE)
lines( xgrid, dbeta( xgrid, shape1= 5, shape2=.5), col="grey")
lines( xgrid, dsplineDensity( fitY, xgrid) , col="red")

# return levels assuming exponential tails

set.seed( 122)
Y<- rgamma( 5000, shape=2.5)
fitY<- splineDensity( Y,  threshold=c(TRUE,FALSE))
returnLevels<- seq( 5,12,,200)
Probs<- psplineDensity( fitY, returnLevels)
plot( 1/( 1-Probs), returnLevels, log="x")



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
