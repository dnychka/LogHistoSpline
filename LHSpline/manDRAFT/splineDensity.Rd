\name{splineDensity}
\alias{splineDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splineDensity(Y, xlim = NULL, threshold = c(FALSE, FALSE), N = 150, lambda = NA, nGrid = 10, verbose = FALSE, debug = FALSE, epsilon = 0, gStart = NULL, expandRange = NULL, whichCV = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{nGrid}{
%%     ~~Describe \code{nGrid} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{gStart}{
%%     ~~Describe \code{gStart} here~~
}
  \item{expandRange}{
%%     ~~Describe \code{expandRange} here~~
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

## The function is currently defined as
function (Y, xlim = NULL, threshold = c(FALSE, FALSE), N = 150, 
    lambda = NA, nGrid = 10, verbose = FALSE, debug = FALSE, 
    epsilon = 0, gStart = NULL, expandRange = NULL, whichCV = 4) 
{
    if (is.na(lambda[1])) {
        lambda <- 10^seq(-2, 2, , nGrid)
    }
    require(fields)
    if (is.null(xlim)) {
        xlim <- range(Y, na.rm = TRUE)
    }
    if (!is.null(expandRange)) {
        xlim[2] <- xlim[2] * expandRange
        N <- round(N * expandRange)
    }
    ind <- Y < xlim[1] | Y > xlim[2] | is.na(Y)
    Y <- Y[!ind]
    binWidth <- (xlim[2] - xlim[1])/(N - 1)
    brk <- seq(xlim[1] - binWidth/2, xlim[2] + binWidth/2, length.out = (N + 
        1))
    binMidpoints <- brk[1:N] + (brk[2] - brk[1])/2
    histObject <- hist(Y, breaks = brk, plot = FALSE)
    fitMgcv <- NULL
    if (debug) {
        require(mgcv)
        dataLocal <- data.frame(x = binMidpoints, counts = histObject$counts)
        fitMgcv <- gam(counts ~ s(x, k = N - 3) + binMidpoints, 
            family = poisson, data = dataLocal, ...)
    }
    fit <- poissonSpline(x = binMidpoints, y = histObject$counts, 
        lambdaGrid = lambda, verbose = verbose, gStart = gStart, 
        whichCV = whichCV)
    info <- fit$info
    lambdaEst <- fit$lambdaEst
    gHat <- fit$gHat
    xr <- range(binMidpoints)
    dx <- xr[2] - xr[1]
    yEnd <- splint(binMidpoints, gHat, c(xr[1], xr[2]))
    yD <- splint(binMidpoints, gHat, c(xr[1], xr[2]), derivative = 1)
    I1 <- I3 <- 0
    if (!threshold[1]) {
        if (yD[1] < 0) {
            warning("left derviative at end point is negative")
        }
        I1 <- exp(yEnd[1])/yD[1]
    }
    if (!threshold[2]) {
        I3 <- -exp(yEnd[2])/yD[2]
        if (yD[2] > 0) {
            warning("right derviative at end point is positive")
        }
    }
    NG <- 40 * N
    xGrid <- seq(xr[1], xr[2], , NG)
    dx0 <- binMidpoints[2] - binMidpoints[1]
    dx <- xGrid[2] - xGrid[1]
    yGrid <- splint(binMidpoints, gHat, xGrid)
    cdfTemp <- c(0, cumsum(0.5 * (exp(yGrid[-1]) + exp(yGrid[-NG]))) * 
        dx)
    I2 <- cdfTemp[NG]
    cdfTemp <- cdfTemp/I2
    nConstant <- I1 + I2 + I3
    probLeft <- I1/nConstant
    probMiddle <- I2/nConstant
    probRight <- I3/nConstant
    cdfTemp <- (probLeft + probMiddle * cdfTemp)
    gHat <- gHat - log(nConstant)
    bTail <- yEnd - log(nConstant)
    aTail <- yD
    cdf <- splint(xGrid, cdfTemp, binMidpoints)
    obj <- list(Y = Y, fitPoisson = fit, hist = histObject, logPdf = list(x = binMidpoints, 
        y = gHat), pdf = list(x = binMidpoints, y = exp(gHat)), 
        cdf = list(x = binMidpoints, y = cdf), cdfRoot = list(x = xGrid, 
            y = cdfTemp), xr = xr, probLeft = probLeft, probRight = probRight, 
        probMiddle = probMiddle, I = c(I1, I2, I3), aTail = aTail, 
        bTail = bTail, threshold = threshold, lambda = lambda, 
        lambdaEst = fit$lambdaEst, fitMgcv = fitMgcv, info = info, 
        xlim = xlim, whichCV = whichCV)
    class(obj) <- "splineDensity"
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
