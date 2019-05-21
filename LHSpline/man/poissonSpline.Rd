\name{poissonSpline}
\alias{poissonSpline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Esimates Poisson spline
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
poissonSpline(x, y, lambdaGrid, gStart = NULL, nIteration = 50, tol = 1e-06,
epsilon = 0, verbose = FALSE, whichCV = 4)
poissonSplineLambda(x, y, gStart, lambda, SigmaInv, epsilon,
 Q, nIteration = 10, tol = 1e-06, verbose = FALSE, gFloor = -15, gCeiling = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{lambdaGrid}{
%%     ~~Describe \code{lambdaGrid} here~~
}
  \item{gStart}{
%%     ~~Describe \code{gStart} here~~
}
  \item{nIteration}{
%%     ~~Describe \code{nIteration} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{whichCV}{
%%     ~~Describe \code{whichCV} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}% use one of  RShowDoc("KEYWORDS")

