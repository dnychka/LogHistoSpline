
fit_bayes_GPKernel<- function(x,
                                   N = 150,
                                   nMCMC = 5000,
                                expandRange = 1.0,
                                tau_beta = 40,
                                #tau_vec = c(25, 25, 4, 25),
                                c_val = 0.85,
                                a0 ,
                                burn.in= .1,
                                verbose = FALSE,
                                b0=NULL,
                                threshold = c( TRUE, FALSE),
                                zlim=NULL, 
                                a= 1.5, # range for GPkernel on [0,1]
                                b= 1.0 # endpoint GPKernel
                                ) {
  z <- x
  n_eff<- length( z)

  countObj<- createCounts( z, N, expandRange= expandRange)
  y<- countObj$y
  
  X<- cbind( 1,countObj$z_center,
             countObj$z_center^2,
             log10(countObj$z_grid)
            )
  
  p <- ncol(X)
  binWidth<- countObj$binWidth
  
  objSigma<- makeWendlandKernel(countObj$z_center, X)
  #objSigma<- makeSigma( countObj$z_center, X)
  vals_k <-objSigma$vals_k
  vecs_k <-objSigma$GPKernel # use raw basis functions. 
  #XVersion<- objSigma$Q1
  
  Q1<- objSigma$Q1
  Q2<- objSigma$Q2
  XVersion<- X
  #XVersion<- Q1
  
  
  glmOffset<- rep(log(n_eff * countObj$binWidth), N)
  glm_init <- glm(
    y ~ XVersion - 1 + offset(glmOffset),
    family = poisson()
     )
  beta <- coef(glm_init)
  betaInit<- beta
  mu_glm <- fitted(glm_init)
  
  raw_gamma <- (y - mu_glm) / sqrt(pmax(mu_glm, 1e-6))
  raw_gamma <- as.vector(raw_gamma - Q1 %*% (t(Q1) %*% raw_gamma))
  
  sg <- sd(raw_gamma)
  raw_gamma <- if (is.finite(sg) && sg > 0){ (raw_gamma / sg)}
  else {
    rep(0, length(raw_gamma))
  }
  
  u_init <- as.vector(t(Q2) %*% (0.02 * raw_gamma))
  a <- as.vector(t(vecs_k) %*% u_init / sqrt(vals_k))
  
  tmp <- build_gamma_from_a(a, Q2, vecs_k, vals_k)
  u <- tmp$u
  gamma <- tmp$gamma
  
  #
  # inverse gamma  mean= b0/( a0-1)
  #                sd = sqrt( b0^2/( (a0-1)^2*(a0-2)))
  
  
  meanPriorSigma2<- b0/( a0-1)
  sigma2 <- b0/( a0-1)
  sdPriorSigma2<- sqrt( b0^2/( (a0-1)^2*(a0-2)))
  
  betaMCMC  <- matrix(NA_real_, nMCMC, p)    # keep for interpretability
  gammaMCMC <- matrix(NA_real_, nMCMC, length(y))
  aMCMC <- matrix(NA_real_, nMCMC, length(a))
  sigma2MCMC <- numeric(nMCMC)
  ess_steps_a <- integer(nMCMC)
  ess_steps_beta <- integer(nMCMC)
  
  ############################################################
  ###### begin MCMC loop
  ############################################################
  for (k in (1:nMCMC) ) {
    #if (verbose && k %% 250 == 0) cat(k, "")
    
    # Gamma/a block by ESS.
    # ESS uses the Gaussian prior through the ellipse, so the slice target
    # must be the log-likelihood only.
    loglik_current <- loglik_gamma_lowrank(
      a, beta, XVersion, y, Q2, 
      vecs_k, 
      vals_k, 
      n_eff,
      countObj$binWidth
    )
    
    nu <- sqrt(sigma2) * rnorm(length(a))
    logy <- loglik_current + log(runif(1))
    theta <- runif(1, 0, 2 * pi)
    theta_min <- theta - 2 * pi
    theta_max <- theta
    n_try <- 0L
    
    repeat {
      #cat(n_try,  a, fill=TRUE)
      n_try <- n_try + 1
      a_prop <- a * cos(theta) + nu * sin(theta)
      
      loglik_prop <- loglik_gamma_lowrank(
        a_prop, beta, XVersion, y, Q2, vecs_k, vals_k, 
        n_eff, 
        countObj$binWidth
      )
      
      if (is.finite(loglik_prop) && loglik_prop >= logy) {
        a <- a_prop
        tmp <- build_gamma_from_a(a, Q2, vecs_k, vals_k)
        u <- tmp$u
        gamma <- tmp$gamma
        break
      }
      
      if (theta < 0) theta_min <- theta else theta_max <- theta
      theta <- runif(1, theta_min, theta_max)
    }
    ess_steps_a[k] <- n_try
    
    # Gibbs update for sigma2 under:
    #   a_l | sigma2 ~ N(0, sigma2)
    #   sigma2 ~ IG(a0, b0)
    sigma2 <- draw_sigma2(a, a0, b0)
    #cat( k, sigma2, fill=TRUE)
    
    # Beta block by ESS.
    # Again, use log-likelihood only because beta prior is the ellipse.
    loglik_current <- loglik_counts(beta, gamma, XVersion, y, n_eff,
                                    countObj$binWidth)
    
    # with coefficient-specific scales:
  
   # nu <- sqrt(tau_vec) * rnorm(p)
    nu <- sqrt(tau_beta) * rnorm(p)
    logy <- loglik_current + log(runif(1))
    theta <- runif(1, 0, 2 * pi)
    theta_min <- theta - 2 * pi
    theta_max <- theta
    n_try <- 0L
    
    repeat {
      n_try <- n_try + 1L
      beta_prop <- beta * cos(theta) + nu * sin(theta)
      
      loglik_prop <- loglik_counts(beta_prop, gamma, 
                                   XVersion, y, n_eff, binWidth)
      
      if (is.finite(loglik_prop) && loglik_prop >= logy) {
        beta <- beta_prop
        break
      }
      
      if (theta < 0) theta_min <- theta else theta_max <- theta
      theta <- runif(1, theta_min, theta_max)
    }
    ess_steps_beta[k] <- n_try
    
    betaMCMC[k, ] <- beta
    gammaMCMC[k, ] <- gamma
    aMCMC[k, ] <- a
    sigma2MCMC[k] <- sigma2
  }
  

  # zlim_obs: where data are and where MCMC should fit
  # zlim: final support for extended density output
  ######################################################
  # DWN  all output w/r to zlim 
  ######################################################
  
  eta_draws_z <- XVersion %*%t(betaMCMC) +  t( gammaMCMC)
  gammaDraws<- t( gammaMCMC)
  # trim output  for burn in
  
  if( !is.null( burn.in)){
    idx<- ( round(nMCMC*burn.in) : nMCMC)
    eta_draws_z<-  eta_draws_z[,idx]
    betaMCMC <- betaMCMC[idx,]
    gammaMCMC <- gammaMCMC[idx,]
    aMCMC <- aMCMC[idx,]
    sigma2MCMC <- sigma2MCMC[idx]
  }
  
  etaHat<- rowMeans( eta_draws_z) 
  gammaHat<- rowMeans( gammaDraws) 
  predictList<- setupDensityEvaluation( etaHat, 
                                        countObj$z_grid,
                                    threshold=threshold
                                    )
  
  out<- c( 
    list(
    z=z, # original data 
    y=y, # histogram possibly with extended bins
    z_grid = countObj$z_grid,
    brk = countObj$brk,
    binWidth = countObj$binWidth,
    X = X,
    XVersion=XVersion,
    Q1=Q1, 
    Q2 = Q2,
    vals_k = vals_k,
    vecs_k = vecs_k,
    zlim_obs = NA,
    nMCMC= nMCMC, 
    eta_draws_z= eta_draws_z,
    # MCMCs
    betaMCMC = betaMCMC,
    gammaMCMC = gammaMCMC,
    aMCMC = aMCMC,
    sigma2MCMC = sigma2MCMC,
    
    # ESS steps
    ess_steps_a = ess_steps_a,
    ess_steps_beta = ess_steps_beta,
  # 
    betaInit= betaInit,
    gammaDraws= gammaDraws, 
    gammaHat= gammaHat, 
    idx= idx,
    gHat = etaHat,
    xr = range(countObj$brk)
  ), 
  predictList

  )
  
}

