\name{poissonSpline}
\alias{poissonSpline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimates a Poisson cubic smoothing spline
%%  ~~function to do ... ~~
}
\description{This is a lean function to estimate the log intensity using a cubic smoothing spline. The model is assumed to be Poission distributed counts at locations and the smoothing parameter (lambda) is estimated by several different approximate versions of cross-validation. This function can be adapted as to also estimate a continuous probability density from a histogram.


}
\usage{
poissonSpline(x, y, lambdaGrid, gStart = NULL, nIteration = 50, tol = 1e-06,
epsilon = 0, verbose = FALSE, whichCV = 4)

poissonSplineLambda(x, y, gStart, lambda, SigmaInv, epsilon,
 Q, nIteration = 10, tol = 1e-06, verbose = FALSE, gFloor = -15, gCeiling = 10)
}
\arguments{
  \item{x}{
  Locations of the counts.

}
  \item{y}{
  The Poisson counts. 

}
  \item{lambdaGrid}{
  A grid of smoothing parameters (typically on a log spacing) to evaluate the cross-validation functions. 

}
  \item{gStart}{
  Initial guess at the log intensity estimates at locations \code{x}. If missing then starting density is an exponential distribution. 

}
  \item{nIteration}{
  Number of iterations to use to find the log spline for a given lambda value. 
}
  \item{tol}{
  Tolerance for converenge for a given lambda. 

}
  \item{epsilon}{
  Precision parameter for prior on the constant and linear terms for the cubic spline. These should be small reflecting a diffuse prior on these parameters. 
}
  \item{verbose}{
  IF TRUE prints out intermediate information. 

}
  \item{whichCV}{
  Selects version of approximate cross-validation to minimize.
}

}
\details{
For fixed lambda find the cubic smoothing spline that solves:
\deqn{  min_{g} - \sum_{i=1}^n l( y_i, \exp (g(x_i)) ) +
\lambda \int (g^{\prime \prime} (x))^2 dx 
}

Here \eqn{ l(y,\alpha)} is the log Poission likelihood for the observation 
\eqn{y} and with intensity \eqn{\alpha}. 

The form of solution is cubic smoothing spline with knots at \eqn{ \{ x_i \} }
and natural boundary conditions. The standard way to find this solution is
by iterative Fisher method of scoring where each iteration is a weighted cubic
smoothing spline fit to psuedo observations. Details can be found in the
reference given below. 

This computation is implemented directly using the \code{splint} cubic spline
function from fields. 
However, for comparison it is also available through the function 
\code{mcgcv:gam } and using the objects in this function would follow the call:

\preformatted{
out<- gam(  y~ s(x, k= N-3) + x ,family=poisson, sp=lambda) 
}

\strong{About the approximate CV criteria} In this function 4 choices are found
and the R code is abstracted below. Here \code{gHatI} is a linear approximation to estimating
g at each observation having omitted it. \code{bPrime = exp(gHat)} from the full estimate
\preformatted{
 # minus log Poission likelihood of CV estimates.     
  CVlnLike<- -1*sum( y*gHatI - exp(gHatI) ) 
  # CV suggested by Green and Silverman 
  CVap<- mean( bPrime*(z- gHat)^2 / (1-diagAMatrix)^2 )
  # Approx log Likelihood from linear approximation
  CVap2<- mean( bPrime*( y - exp(gHatI) )^2 )
  # leaving out the weights 
  CVap3<- mean(( y - exp(gHatI) )^2 )
}
}

\value{
\item{x}{Data locations}
\item{y}{Data counts}
\item{gHat}{Estimated log intensities at \code{x} at cross-validation estimate of \code{lambda} }
\item{lambdaGrid}{Grid of lambdas searched}
\item{lambdaEst}{Lambda that minimizes CV criterion}

\item{info}{A table summarizing CV criteria for each value of lambda in the grid.
info[,5] is the number of iterations for convergence, info[,6] approximate degrees of freedom in the curve estimate based on linearing the problem at convergence. 
}
\item{z}{Pseudo data at convergence. \code{z<- (y - exp(gHat))/ exp(gHat) + gHat } }
\item{fit}{Details from fitting at the CV lambda}
\item{fitted.values}{\code{exp(gHat)}}
\item{gStart}{Starting values used.}

}
\references{
See
Green, Peter J., and Bernard W. Silverman. Nonparametric regression and generalized linear models: a roughness penalty approach. CRC Press, 1993  and specially  page 100.
}
\author{
Douglas Nychka
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{splineDensity}}, \code{\link{splint}}
}
\examples{
set.seed(123)
x<- seq( 0,1,length.out= 150)
alphaTrue<-  exp( 6*x^2*(1-x))
alphaTrue<-  25*alphaTrue/ max( alphaTrue)
gTrue<- log(alphaTrue )
y<- rpois(length( alphaTrue), lambda=alphaTrue)
lambdaGrid<- 10**( seq( -3, 0,length.out=50))
fit1<- poissonSpline( x,y, lambdaGrid= lambdaGrid)
set.panel(3,1)
plot( x,y)
lines( x, fit1$fitted.values)
title("Counts and fitted curve")

plot( x,fit1$gHat, type="l" )
lines(x , gTrue, col="red")
title("Estimated and true (red) log intensity functions")

plot( fit1$info[,6], fit1$info[,fit1$whichCV] ,
type="b",
pch=16,cex=.5,
xlab="effective df", ylab="CV function")
xline( (fit1$fit$trA))
title("CV criterion with minimum")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}% use one of  RShowDoc("KEYWORDS")

