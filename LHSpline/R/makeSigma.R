

makeSigma<- function(z_center,
                     X){
p<- ncol(X)
QR <- qr(X) 
bigQ <- qr.Q(QR, complete = TRUE)

Q1 <- bigQ[, (1:p)]
Q2 <- bigQ[, (p + 1):nrow(X)]

GPKernel<- Matern( rdist( z_center, z_center), aRange= .5, 
                   smoothness =1.5)

 K<- GPKernel%*%t( GPKernel)
 Sigma <- t(Q2) %*% K %*% Q2
 Sigma_scale <- mean(diag(Sigma))
 Sigma <- (Sigma + t(Sigma)) / 2 / Sigma_scale
 Sigma <- Sigma + 1e-6 * diag(nrow(Sigma))
#
 eig <- eigen(Sigma, symmetric = TRUE)
 vals <- eig$values
 vecs <- eig$vectors

# keep_pos <- vals > 1e-10
# vals <- vals[keep_pos]
# vecs <- vecs[, keep_pos, drop = FALSE]
# 
cumprop <- cumsum(vals) / sum(vals)
r99 <- which(cumprop >= 0.999)[1]
r_k <- min(length(vals), max(10, r99))

vals_k <- vals[(1:r_k)]
vecs_k <- vecs[, (1:r_k)]


return(
  list( 
    Q1= Q1,
    Q2=Q2,
    vals_k = vals_k,
    vecs_k = vecs_k,
    transform_z= z_center
   )
  )
}