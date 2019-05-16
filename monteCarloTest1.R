library( extRemes)

pR<- 1/(c(5,10,25,50,100)*365)
shapePar<- 3
scalePar<- 1/2
N<- 10*365
returnTrue<- qgamma( pR, shape=shapePar, scale=scalePar, lower.tail=FALSE)
tailProb<- .01
nReps<- 200
out<- matrix( NA, ncol=length(pR), nrow=nReps)
out2<- matrix( NA, ncol=length(pR), nrow=nReps)
set.seed( 124)
for( k in 1:nReps){
  cat(k," ")
Y<- rgamma(N, shape=shapePar, scale=scalePar)
threshold<- quantile( Y, 1- tailProb)
GPFit<- fevd( Y, threshold=threshold, type="GP",
              method="MLE")
out[k,]<- c(return.level(GPFit, c(5,10,25,50,100), do.ci=FALSE ))
}

cat(" ", fill=TRUE)
error<-  (t(t( out)/ returnTrue) -1)*100
stats(error )
bplot( error)

