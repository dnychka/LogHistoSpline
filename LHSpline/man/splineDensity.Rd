\name{splineDensity}
\alias{splineDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Log density estimate as a cubic spline
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splineDensity(Y, xlim = NULL,
 threshold = c(FALSE, FALSE), N = 150, lambda = NA, nGrid = 10, verbose = FALSE, debug = FALSE, epsilon = 0, gStart = NULL, expandRange = NULL, whichCV = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Data to fit a density function
%%     ~~Describe \code{Y} here~~
}
  \item{xlim}{ Limits for the support of the density. 
  If not specified these are set from the data range and expandRange.
%%     ~~Describe \code{xlim} here~~
}
  \item{threshold}{ 
   
%%     ~~Describe \code{threshold} here~~
}
  \item{N}{ Number of histogram bins. 
%%     ~~Describe \code{N} here~~
}
  \item{lambda}{ Cubic spline smoothing parameter.
%%     ~~Describe \code{lambda} here~~
}
  \item{nGrid}{Number of lambda values used in initial grid search. 
%%     ~~Describe \code{nGrid} here~~
}
  \item{verbose}{ If TRUE then print out intermediate information.
%%     ~~Describe \code{verbose} here~~
}
  \item{debug}{
  If TRUE print out much intermediate objects for debugging.
%%     ~~Describe \code{debug} here~~
}
  \item{epsilon}{
  
%%     ~~Describe \code{epsilon} here~~
}
  \item{gStart}{
  If supplied the starting values for the log density at the hisotgram midpoints. This is useful if one has the results from an estimate that differs but just the value of lamdba used. 
%%     ~~Describe \code{gStart} here~~
}
  \item{expandRange}{
  Expands the range of the historgram beyond the range of the data.(And fills these bins with zero counts.)
%%     ~~Describe \code{expandRange} here~~
}
  \item{whichCV}{
  Determines which approximate cross-validation criterion to use for minimization. 
%%     ~~Describe \code{whichCV} here~~
}
}
\details{
This function fits a density to histogram counts. It assumes that the counts follow a Poisson distribution and the log intensity is represented as a cubic smoothing spline with knots at the histogram mid points. lambda is the usual smoothing parameter. The final estimate is normalized to integrate to one and also setup so that the density outisde the range of the data follows an exponential distribution. This extrapolation is consistent with the spline extrapolating linearly as the log density. 


%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class splineDensity.
Similar to the families of parametric density there are functions for generating a random sample, evaluating the density, cumulative distribution function and the quantile function. 

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Huang, Whitney K., Douglas W. Nychka, and Hao Zhang. "Estimating precipitation extremes using the log‐histospline." Environmetrics 30.4 (2019): e2543.
}
\author{
Douglas Nychka
%%  ~~who you are~~
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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

xlim=c(0,1) )
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
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
