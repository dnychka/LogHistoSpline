makeWendlandKernel<- function(z_center,
                              X, 
                              fracClip=.75,
                              NThin=5, 
                              overlap= 3*NThin){
  p<- ncol(X)
  QR <- qr(X)
  bigQ <- qr.Q(QR, complete = TRUE)
  
  Q1 <- bigQ[, 1:p , drop = FALSE]
  Q2 <- bigQ[, (p + 1):nrow(X), drop = FALSE]
  
  N<- length( z_center)
  thinIndex<- seq( 1, round(N*(1-fracClip)), NThin)
  binWidth<- z_center[2] - z_center[1]
  GPKernel<- Wendland( rdist(z_center, z_center[thinIndex]),
                       dimension=2, k=2, aRange=overlap*binWidth)
  vals_k<- rep( 1, ncol( GPKernel))
  vecs_k =  t(Q2)%*%GPKernel
 out<- 
   list( 
     Q1 = Q1,
     Q2 = Q2,
     vals_k = vals_k,
     vecs_k = vecs_k,
     transform_z= z_center,
     GPKernel= GPKernel
   )
  return(out)
}