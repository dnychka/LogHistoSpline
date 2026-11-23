
splintLinear<- function(x,y, xGrid ){
  N<- length( x)
yGrid<- approx( x, y, xout=xGrid)$y
ind<- xGrid < x[1]
bslope <- (y[2]- y[1])/ (x[2]- x[1])
yGrid[ind]<- bslope*(xGrid[ind]- x[1]) + y[1]

ind<- xGrid >  x[N]
bslope <- (y[N]- y[N-1])/ (x[N]- x[N-1])
yGrid[ind]<- bslope*(xGrid[ind]- x[N]) + y[N]
return(yGrid )
}