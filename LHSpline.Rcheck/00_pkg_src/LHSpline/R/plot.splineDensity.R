
plot.splineDensity<- function( obj, 
                               which=c( TRUE, TRUE, TRUE, TRUE),
                               truePdf=NULL,
                               xlim =NULL, ylim=NULL
){
  #  set.panel( 3,1)
  if( which[1]){
    plot( obj$hist, freq=FALSE )
    lines(obj$pdf, co="blue", lwd=2)
    if( !is.null( truePdf)){
      lines( truePdf, col="grey20", lty=2, lwd=2)
    }
    
  }
  if( which[2]){
    plot(obj$logPdf, col="blue", type="l", lwd=2,
         xlim=xlim, ylim =ylim,
         ylab="log density function", xlab=" x ")
    xline( range(obj$Y), col="red")
    if( !is.null( truePdf)){
      lines( truePdf$x,log(truePdf$y), col="grey20", lty=2, lwd=2 )
    }
  }
  ## QQ plot 
  if( which[3]){
    N<- length( obj$Y)
    sortY<- sort(obj$Y)
    pY<- psplineDensity( obj, sortY)
    pGrid<-  ( (1:N) -.5)/N
    QY<- qsplineDensity(obj, pGrid)
    plot( QY, sortY, pch=16,
          xlab= "Fitted Quantiles",
          ylab="Sorted Sample", col="grey30")
    abline( 0,1,col="red", lwd=2) 
  }
  
  if( which[4]& (nrow(obj$info) > 1 ) ){
    whichCV<- obj$whichCV
    CVlabels<- c("-lnLike", "CVGreenSilverman",
                 "CVIntensity", "CVUnweighted")
    plot( obj$info[,1], obj$info[,whichCV], type="b",
          lwd=1, col="grey", pch=16, log="x",
          main= CVlabels[whichCV-2], 
          ylab="CV score",
          xlab="log lambda"
          )
    xline( obj$lambdaEst)
  }
}
