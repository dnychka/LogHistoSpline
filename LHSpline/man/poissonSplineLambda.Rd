\name{poissonSplineLambda}
\alias{poissonSplineLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
poissonSplineLambda(x, y, gStart, lambda, SigmaInv, epsilon, Q, nIteration = 10, tol = 1e-06, verbose = FALSE, gFloor = -15, gCeiling = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{gStart}{
%%     ~~Describe \code{gStart} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{SigmaInv}{
%%     ~~Describe \code{SigmaInv} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{Q}{
%%     ~~Describe \code{Q} here~~
}
  \item{nIteration}{
%%     ~~Describe \code{nIteration} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{gFloor}{
%%     ~~Describe \code{gFloor} here~~
}
  \item{gCeiling}{
%%     ~~Describe \code{gCeiling} here~~
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

## The function is currently defined as
function (x, y, gStart, lambda, SigmaInv, epsilon, Q, nIteration = 10, 
    tol = 1e-06, verbose = FALSE, gFloor = -15, gCeiling = 10) 
{
    gOLD <- gStart
    N <- length(y)
    Omega <- matrix(0, N, N)
    Omega[3:N, 3:N] <- lambda * SigmaInv
    Omega[1, 1] <- epsilon
    Omega[2, 2] <- epsilon
    Gamma <- Q \%*\% Omega \%*\% t(Q)
    for (iteration in 1:nIteration) {
        gOLD <- ifelse(gOLD > gCeiling, gCeiling, gOLD)
        gOLD <- ifelse(gOLD < gFloor, gFloor, gOLD)
        bPrime <- exp(gOLD)
        z <- (y - bPrime)/bPrime + gOLD
        weights <- bPrime
        if (verbose) {
            print(lambda)
            print(stats(cbind(gOLD, weights, z)))
        }
        gNEW <- splint(x, z, xgrid = x, wt = weights, lambda = lambda/N)
        testTol <- sqrt(mean((gNEW - gOLD)^2)/mean(gOLD^2))
        if (verbose) {
            cat("iteration: ", iteration, "  ", sqrt(mean((gNEW - 
                gOLD)^2)), testTol, fill = TRUE)
        }
        if (testTol < tol) 
            break
        gOLD <- gNEW
    }
    gHat <- gNEW
    hHat <- exp(gNEW)
    beta <- Q \%*\% gHat
    bPrime <- exp(gHat)
    weights <- bPrime
    W <- diag(weights)
    z <- (y - bPrime)/bPrime + gHat
    diagAMatrix <- diag(solve(W + Gamma) \%*\% W)
    CVres <- (z - gHat)/(1 - diagAMatrix)
    gHatI <- z - CVres
    CVlnLike <- -1 * sum(y * gHatI - exp(gHatI))
    CVap <- mean(bPrime * (z - gHat)^2/(1 - diagAMatrix)^2)
    CVap2 <- mean(bPrime * (y - exp(gHatI))^2)
    CVap3 <- mean((y - exp(gHatI))^2)
    lnLikeLinear <- NA
    lnLike <- NA
    if (any(gHat == gFloor)) {
        warning("gHat set to the lower bound at some values")
    }
    if (any(gHat == gCeiling)) {
        warning("gHat set to the upper bound at some values")
    }
    return(list(nit = iteration, gHat = gHat, fitted.values = exp(gHat), 
        lambda = lambda, beta = beta, lnLike = lnLike, CVlnLike = CVlnLike, 
        CVap = CVap, CVap2 = CVap2, CVap3 = CVap3, epsilon = epsilon, 
        lnLikeLinear = lnLikeLinear, z = z, gStart = gStart))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
