print.splineDensity<- function( obj, ...){
  
    cat("Data summary:", fill=TRUE)
     temp<- as.data.frame( t(stats( obj$Y)) )
     row.names( temp)<- " "
     print( temp)
      cat(" ", fill=TRUE)
      N<- length( obj$hist$mids)
      cat("Histogram range:", obj$xr,fill=TRUE)
      cat("Number of bins:", N, fill=TRUE)
      cat(" ", fill=TRUE)
      whichCV<- obj$whichCV
      CVlabels<- c("-lnLike", "CVGreenSilverman",
                   "CVIntensity", "CVUnweighted")
      cat("Smoothing parameter (lambda):", obj$lambdaEst,
          fill=TRUE)
      cat("Found by:", CVlabels[whichCV], fill=TRUE)
      
      cat( "thresholds (Right, Left):",
           obj$threshold, fill=TRUE)
      
     }
