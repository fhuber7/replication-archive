bvar.mix <- function(Xraw, 
                     nsave = 500,     # Nsave is the number of saved posterior draws.
                     nburn = 500,     # Nburn is the number of burn-in iterations.
                     G.max = 15,      # G.max is the maximum number of groups.
                     p = 2,           # P is the lag order of the VAR model.
                     SV.idi = FALSE,  # SV.idi is the indicator for stochastic volatility.
                     fhorz = 12,      # Fhorz is the forecast horizon.
                     irf.hor = 24,    # Irf.hor is the impulse response function horizon.
                     NG.idi = FALSE,  # NG.idi is the indicator for Normal-Gamma prior. Put to FALSE (not thoroughly tested)
                     prior.cov = "AR",# Prior.cov is the prior covariance structure.
                     decomp = "chol", # Decomp is the decomposition method. Chol for Cholesky or Eigen 
                     sample.w = "conditional", # Sample.w is the sampling method for the random effects. 
                     flat = FALSE     # Flat is the indicator for flat prior on the VAr coefficients
) 
  {
  require(mvtnorm)
  require(mvnfast)
  ntot <- nsave+nburn

  #Creates design matrices for the VAR
  M <- ncol(Xraw)
  Xraw <- embed(Xraw,dimension = p+1)
  Y <- Xraw[, 1:M]
  X <- cbind(Xraw[, (M+1):((p+1)*M)], 1)
  
  # Some dimensions such as the length of the sample in the estimation part and the number of covariates in X
  T <- nrow(X)
  K <- ncol(X)
  
  # Computes quantities that are useful later on or initialize the sampler
  A <- solve(crossprod(X) + 1e-4 * diag(K))%*%crossprod(X, Y)
  Sigma <- crossprod(Y - X%*%A)/(T-ncol(X))
  ystar <- Y - X%*%A
  
  # Here we set a few prior hyperparameters
  A.prior <- matrix(0, K, M)
  A.prior[1:M, 1:M] <- diag(M)*0.8 #Prior mean for the VAR coefficients. Forces the series towards persistence
  
  v_1 <- v_2 <- 0.6 #Shrinkage hyperparameters for the Normal-Gamma prior on the common mean term
  
  #Prior on alpha (standard choice stipulated in Escobar & West)
  a_alpha <- 2
  b_alpha <- 4
  
  #Prior on the the component means and variances
  b0 <- rep(0, M) #This is a draw for the common mean of the Gaussian distribution of the mixtures
  B0 <- diag(M)*1e-1 #Prior variance on the common mean
  B0inv <- solve(B0)
  
  delta <- 10^2
  D0 <- delta*diag(M) #Prior variance on the prior on the common mean
  
  #Prior on the Wishart components
  c0 <- 2*(2.5 + (M-1)/2)
  g0 <- 2*(0.5 + (M-1)/2)
  C_0 <- diag(M)
  
  #Construct scaling matrix for common distribution of the Sigmas using an AR(p)  regression
  if (prior.cov == "SFS" || prior.cov == "AR"){
    G0 <- diag(M)
    for (jj in 1:M) G0[jj,jj] <- summary(lm(Y[(p+1):T,jj]~embed(Y[,jj], dimension = p+1)[, 2:(p+1)]-1))$sigma^2
    Omega <- G0
    C_0 <- solve(G0) 
  }else{
    Omega <- diag(M)*1e-4
    C_0 <- Omega
    G0 <- Omega
  }
  
  
  #Starting values for the prior on A
  Theta <- matrix(1e-1, M, K) 
  B.j <- rep(1, M)
  
  #Initialize stuff related to the mixture model
  G <- G.max 
  kappa <- 0.8 #truncation parameter for the slice sampler 
  xi <- (1-kappa)*kappa^(seq(1, G)-1)
  #--------------------------Prior hyperparameters --------------------------
  alpha <- 0.5 # initialize at expected value
  a_k <- 1
  b_k <- alpha
  
  # Initial values for the prior
  a_tau <- matrix(1, M, 1)
  count_a <- matrix(0, M, 1)
  scale_a <-  matrix(0.1, M, 1)
  #-----------------------------Starting values for the infinite mixture -----------------------------
  if (G.max > 1){
    
    # Classification
    si <- matrix(0, T, 3)
    for (t in 1:T) si[t, sample(1:3,1)] <- 1 
    S <- apply(si, 1, function(x) which(x==max(x)))
    
    mu.G <- matrix(0, M, G)
    Sigma.G <- array(0, c(M,M, G))
    Sigmainv.G <- array(0, c(M, M, G))
    Q.G <- array(0, c(M,M, G))
    for (jj in seq_len(ncol(si))){
      mu.G[ , jj] <- apply(ystar[S==jj,],2,mean)
      Sigma.G[ , , jj] <- cov(ystar[S==jj,]) + diag(M)
      if (decomp == "chol"){
        Q.G[ , , jj] <- t(chol(Sigma.G[,,jj]+diag(M), pivot = TRUE))
      }else{
        eigs.Sig <- eigen(Sigma.G[,,jj])
        Q.G[ , , jj] <- eigs.Sig$vectors%*%diag(sqrt(eigs.Sig$values))
      }
      Sigmainv.G[ , , jj] <- solve(Sigma.G[,,jj]+diag(M))
    }
    
    counts <- matrix(0, G, 1)
    S.tab <- table(S)
    counts[as.numeric(names(S.tab))] <- S.tab
    Nclust <- length(counts)	  # initial number of clusters
    
    mu.t <- matrix(0, T, M)
    Q.t <- array(0, c(T, M, M))
    
    for (t in seq_len(T)){
      s.i <- S[[t]]
      mu.t[t, ] <- mu.G[ , s.i]
      if (decomp == "chol"){
      Q.t[t, , ] <- t(chol(Sigma.G[ , , s.i], pivot=TRUE))
      }else{
        eigs.Sig <- eigen(Sigma.G[,,s.i])
        Q.t[ t, , ] <- eigs.Sig$vectors%*%diag(sqrt(eigs.Sig$values))
      }
    }
    
    nu <- rbeta(G,a_k, b_k)
    nu[[G]] <- 1
    
    eta <- rep(1, G)
    for (j in 1:G){
      if (j==1){
        eta[[j]] <- nu[[j]]
      }else if (j==2){
        eta[[j]] <- (1-nu[[j-1]])*nu[[j]]
      }else if (j>2){
        eta[[j]] <-  nu[[j]]*prod(1-nu[1:(j-1)])
      }
    }
  }else{
    counts <- matrix(T, 1, 1)
    S <- matrix(1, T, 1)
    
    mu.G <- matrix(0, M, G)
    Sigma.G <- array(0, c(M,M, G))
    Sigmainv.G <- array(0, c(M, M, G))
    Q.G <- array(0, c(M,M, G))
    for (jj in 1:G.max){
      mu.G[ , jj] <- apply(ystar[S==jj,],2,mean)
      Sigma.G[ , , jj] <- cov(ystar[S==jj,])
      
      if (decomp == "chol"){
        
        Q.G[ , , jj] <- t(chol(Sigma, pivot=TRUE))
        
      }else{
        
        eigs.Sig <- eigen(Sigma)
        Q.G[ , , jj] <- eigs.Sig$vectors%*%diag(sqrt(eigs.Sig$values))
        
      }
      Sigmainv.G[ , , jj] <- solve(Sigma)
    }
  }
  
  v <- matrix(rnorm(T*M, 0.01), T, M)
  w <- matrix(0, T, M)
  
  # Create storage matrices used to store MCMC output
  A.store <- array(NA, c(nsave, K, M))
  mu.store <- array(NA, c(nsave, M, G.max))
  Sigma.store <- array(NA, c(nsave, M, M, G.max))
  Omega.store <- array(NA, c(nsave, M))
  fit.store <- array(NA, c(nsave, T, M))
  G.store <- array(NA, c(nsave, 2),dimnames=list(NULL, c("G", "G+")))
  S.store <- array(NA, c(nsave, T))
  eta.store <- array(NA, c(nsave, G.max))
  vola.store <- array(NA, c(nsave, T, M))
  count.store <- array(NA, c(nsave, G.max))
  OT.store <- array(NA, c(nsave, M))
  vola.vola.store <- array(NA, c(nsave, M))
  
  if (fhorz > 0){
    fcst.store <- array(NA, c(nsave, M, fhorz))
  }else{
    fcst.store <- NULL
  }
  
  if (irf.hor > 0){
    IRF.store <- array(NA, c(nsave, M, M,irf.hor, G.max))
  }else{
    IRF.store <- NULL
  }
  
  acc <- 0
  
  if (SV.idi){
    para.list <- list()
    prior.list <- list()
    C.mat <- matrix(10, M, 1)
    for (jj in seq_len(M)){
      startpara <- list(mu = 0, phi = 0.9, sigma = 0.1,
                        nu = Inf, rho = 0, beta = NA,
                        latent0 = 0)
      
      prior.list[[jj]] <- specify_priors(
        mu = sv_normal(mean = 0, sd = 100),
        phi = sv_beta(shape1 = 25, shape2 = 1e-5),
        sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*C.mat[jj])),
        nu = sv_infinity(),
        rho = sv_constant(0),
        latent0_variance = "stationary",
        beta = sv_multinormal(mean = 0, sd = 10000, dim = 1)
      )
      para.list[[jj]] <- startpara
    }
    H <- log(t(matrix(rep(diag(Omega), T), M, T)))
    H.store <- array(NA, c(nsave, T, M))
    
    # Some starting values used for the shrinkage factors in the log-vola equation
    a_vola <- 0.1
    count_vola <- 0
    scale_vola <- 0.25
    
  }else{
    H <- log(t(matrix(rep(diag(Omega), T), M, T)))
  }
  
  for (irep in seq_len(ntot)){
    #Start sampling the mixture model
    #Step 1: Sample mu.k and sd.K
    ystar <- Y - X%*%A - v #Conditional on subtracting the measurement errors here
  
    #Step 1a:) Sample the component moments of the infinite mixture
    for (s in seq_len(G)){
      if (counts[s,1]>0){
        #Samples the component means from multivariate Gaussian posteriors
        T_k <- sum(S==s)
        
        if (sum(S==s) > 1) mu.OLS <-  apply(ystar[S==s,],2,mean)  else mu.OLS <- ystar[S==s,]
        V.mu <- solve(Sigmainv.G[,,s] * T_k + B0inv)
        mu.mean <- V.mu %*% (B0inv %*% b0 + (Sigmainv.G[,,s] * T_k) %*%mu.OLS)
        mu.mean.draw <- try(mu.mean + t(chol(V.mu)) %*% rnorm(M), silent=TRUE)
        if (is(mu.mean.draw, "try-error")) mu.mean.draw <- mu.mean + t(chol(V.mu, pivot=TRUE)) %*% rnorm(M)
        
        mu.G[, s] <- mu.mean.draw
        
        #Samples the component variances
        yss <- (ystar[S==s,] - mu.G[,s])
        if (sum(S==s) >1) SSEwish <- crossprod(yss) else SSEwish <- tcrossprod(yss)
        
        C_k <- 1/2* SSEwish + C_0
        c_k <- c0 + T_k/2
        
        Sigmainv_k <- rWishart(1, c_k, solve(C_k))[,,1] #Draw from the Wishart posterior for Sigma^-1
        Sigma_k <- solve(Sigmainv_k) #Recover Sigma
        
        if (decomp == "chol"){
          Q_k <- try(t(chol(Sigma_k)))
        
          if (is(Q_k, "try-error")) Q_k <- t(chol(Sigma_k, pivot=TRUE))
        }else{
          eigen.k <- eigen(Sigma_k)
          Q_k <- eigen.k$vectors%*%diag(sqrt(eigen.k$values))
        }
        Q.G[,,s] <- Q_k
        
        Sigma.G[, , s] <- Sigma_k
        Sigmainv.G[, , s] <- Sigmainv_k
      }else{
        #Sample from the prior if no observation is allocated
        mu.G[, s] <- b0 + rnorm(M, 0, diag(sqrt(B0)))
        
        #Samples the component variances
        C_k <- C_0
        c_k <- c0 
        
        Sigmainv_k <- rWishart(1, c_k, solve(C_k))[,,1]# Draw from the posterior (which equals the prior)
        Sigma_k <- solve(Sigmainv_k)
        
        if (decomp == "eigen"){
          Q.G[,,s] <- t(chol(Sigma_k))
        }else{
          eigen.k <- eigen(Sigma_k)
          Q.G[,,s] <- eigen.k$vectors%*%diag(sqrt(eigen.k$values))
        }
        Sigma.G[, , s] <- Sigma_k
        Sigmainv.G[, , s] <- Sigmainv_k
      }
    }
    
    if (prior.cov == "SFS"){
      #Step 1a:) Estimate the common mean matrix for the Wishart prior
      g1 <- g0 + G*c0
      G1 <- G0 + apply(Sigmainv.G, c(1,2), sum)
      C_0 <- bayesm::rwishart(2*g1, 0.5*G1)$W#rWishart(1, g1, solve(G1))[,,1]
    }
    #Step 1b:) Sample the common means for the priors
    #Common prior mean
    V.b0 <- 1/G *B0
    b0 <- apply(mu.G, 1, mean) + rnorm(M, 0, sqrt(diag(V.b0)))
    
    #Step 1c:) Sample the prior variances for the common Gaussian prior
    for (jj in seq_len(M)){
      p_0 <- 2*v_2
      p_K <- v_1 - G/2
      
      z_j <- sum((mu.G[jj,] - b0[jj])^2)
      
      B0[jj,jj] <- GIGrvg::rgig(1, p_K, z_j+1e-8, p_0)
    }
    B0.diag <- diag(B0)
    
    # This is just for stability reasons
    B0.diag[B0.diag < 1e-7] <- 1e-7
    B0.diag[B0.diag > 3] <- 3
    # ---
    
    B0inv <- diag(1/B0.diag)
    
    
    #Step 2:  Sample the sticks from beta distributions
    XA <- (X%*%A)
    if (G.max > 1){
      for (k in seq_len(G-1)){
        #This block samples the sticks
        a_nuk <- a_k +   counts[k,1]
        b_nuk <- b_k + sum(counts[(k+1):G, 1])
        nu.k <- rbeta(1, a_nuk, b_nuk) 
        
        nu[[k]] <- nu.k
        
        #This block exploits the stick-breaking representation (see Eq. 5 in FS & MW)
        if (k==1){
          eta[[k]] <- nu[[k]]
        }else if (k==2){
          eta[[k]] <- (1-nu[[k-1]])*nu[[k]]
        }else if (k>2){
          eta[[k]] <-  nu[[k]]*prod(1-nu[1:(k-1)])
        }
      }
      eta[[G]] <- prod(1-nu[1:(G-1)])
      
      #Step 3: Sample the classifications S
      u <- rep(0, T)
      fit <- matrix(0, T, M)
      epsilon <- matrix(0, T, M)
      shock.epsilon <- matrix(0, T, M)
      for (t in seq_len(T)){
        s.i <- S[[t]]
        xi.t <- xi[s.i]
        u.i <- runif(1, 0, xi.t)
        
        lik <- matrix(0, G, 1)
        for (jj in 1:G){
	        lik.jj <- mvnfast::dmvn(ystar[t,], mu.G[,jj], Sigma.G[,,jj]+ 1e-6*diag(M), log=TRUE)#try(mvnfast::dmvn(ystar[t,], mu.G[,jj], Sigma.G[,,jj], log=TRUE), silent=TRUE)
          lik[jj,] <- lik.jj
        } 
      
        prob.non <- log((u.i < xi[1:G])/xi[1:G] * eta[1:G]) + lik #+ log(1e-32)
        prob.non <- exp(prob.non - max(prob.non))/sum(exp(prob.non - max(prob.non)))
        
        s.i <- sample(1:G,1, prob=prob.non/sum(prob.non), replace=TRUE)
        
        scale <- exp(-.5*H[t,]) # Scaling by the idiosyncratic variances
        
        #Samples the latent factors w(t) either through a regression representation or through the full cond. distr
        if (sample.w == "regression"){
          yss <- (Y[t, ] - XA[t,] - mu.G[,s.i]) * scale 
          xss <- t(t(Q.G[,,s.i]) * scale) 
          
          V_u <- try(solve(crossprod(xss) + diag(M)), silent=TRUE)
          if (is(V_u, "try-error")) V_u <- MASS::ginv(crossprod(xss) + diag(M))
          mu_u <- V_u %*% crossprod(xss, yss)
          
          w.t <- try(t(mu_u + t(chol(V_u))%*%rnorm(M)), silent=TRUE)
          if (is(w.t, "try-error")) w.t <- t(mu_u + t(chol(V_u, pivot=TRUE))%*%rnorm(M))
          w[t,] <- w.t
          epsilon[t,] <- Q.G[,,s.i]%*%w[t,] + mu.G[,s.i]
        }else{
          #  Sample from the full conditionals
          # epsilon_t | y_t ~ N( mu_t, V_t)
          # mu_t = Sigma (Sigma + H)^-1 (yt - B X(t))
          # V_t = Sigma - Sigma (Sigma + H)^-1 Sigma
          
          yss <- (Y[t, ] - XA[t,] - mu.G[,s.i])
          V_prec <- try(Sigma.G[,, s.i] %*% solve(Sigma.G[,, s.i] + diag(exp(H[t,]))), silent=TRUE)
          if (is(V_prec, "try-error")) V_prec <- Sigma.G[,, s.i] %*% MASS::ginv(Sigma.G[,, s.i] + diag(exp(H[t,])))
          V_u <- Sigma.G[,,s.i] -  V_prec %*% Sigma.G[,, s.i]
          mu_u <-   mu.G[,s.i] + V_prec %*% yss
	        e.draw <- try(as.numeric(mu_u + t(chol(V_u + diag(M)*1e-5))%*%rnorm(M)), silent=TRUE)
	        if (is(e.draw, "try-error")) e.draw <- as.numeric(MASS::mvrnorm(1, mu_u, V_u, tol=1e-5))

          epsilon[t, ] <- e.draw
          w[t, ] <- solve(Q.G[,,s.i])%*%epsilon[t, ] #maps back into the standard normal terms
        }
        
        mu.t[t,] <- mu.G[,s.i]
        Q.t[t,,] <- Q.G[,,s.i]
        shock.epsilon[t, ] <- mu.t[t,] + Q.t[t,,] %*%rnorm(M) #Only for computing in-sample predictive densities
        
        
        S[[t]] <- s.i
        u[[t]] <- u.i
      }
      
      S.tab <- table(S)
      counts <- matrix(0,G.max,1)
      counts[as.numeric(names(S.tab))] <- S.tab
      G.plus <- sum(counts>0)
      
      #Does the truncation of G
      G <- sum((1-cumsum(eta)) < min(u))
      if (G.max == 1) G <- 1
      if (G==0) G <- 1
      
      #Step 4: Sample the precision parameter through an MH step
      l_alpha_star <- log(alpha) + rnorm(1,0,0.25)
      
      post.old <- G.plus*log(alpha)+(lgamma(alpha+1) - log(alpha))-lgamma(T+alpha)+sum(lgamma(counts[counts>0]))  + dgamma(alpha, a_alpha, b_alpha, log=TRUE)
      post.new <- G.plus*l_alpha_star+(lgamma(exp(l_alpha_star)+1) - l_alpha_star)-lgamma(T+exp(l_alpha_star))+sum(lgamma(counts[counts>0]))  + dgamma(exp(l_alpha_star), a_alpha, b_alpha, log=TRUE)
      
      acc <- post.new - post.old
      if (is.nan(acc)) acc <- 0
      
      if (acc > log(runif(1))){
        alpha <- exp(l_alpha_star)
      }
      b_k <- alpha
    }else{
      epsilon <- matrix(0, T, M)
      shock.epsilon <- matrix(0, T, M)
      G.plus <- 1
      if (sample.w == "regression"){ 
        #Corresponds to sampling the random effects using a regression interpretation
        for (t in seq_len(T)){
          scale <- exp(-.5*H[t,])
          yss <-  (Y[t, ] - XA[t,] - mu.G) * scale
          xss <- t(t(Q.G[,,1]) * scale)
          
          V_u <- try(solve(crossprod(xss) + diag(M)), silent=TRUE)
          if (is(V_u, "try-error")) V_u <- MASS::ginv(crossprod(xss) + diag(M))
          
          mu_u <- V_u %*% crossprod(xss, yss)
          w.draw <- try(t(mu_u + t(chol(V_u))%*%rnorm(M)), silent=TRUE)
          if (is(w.draw, "try-error")) w.draw <- t(mu_u + t(chol(V_u, pivot=TRUE))%*%rnorm(M))
          w[t,] <- w.draw
          epsilon[t,] <- Q.G[,,1]%*%w[t,] + mu.G
          shock.epsilon[t, ] <-mu.G + Q.G[,,1] %*%rnorm(M)
        }
      }else{ #Samples the random effects by exploiting properties of the conditional distributions
        # Sample from the full conditionals
        # epsilon_t | y_t ~ N( mu_t, V_t)
        # mu_t = Sigma (Sigma + H)^-1 (yt - B X(t))
        # V_t = Sigma - Sigma (Sigma + H)^-1 Sigma
        
        for (t in seq_len(T)){
          yss <- (Y[t, ] - XA[t,] - mu.G[,1])
          V_prec <- try(Sigma.G[,, 1] %*% solve(Sigma.G[,, 1] + diag(exp(H[t,]))), silent=TRUE)
          if (is(V_prec, "try-error")) V_prec <-Sigma.G[,, 1] %*% MASS::ginv(Sigma.G[,, 1] + diag(exp(H[t,])))
          V_u <- Sigma.G[,,1] -  V_prec %*% Sigma.G[,, 1]
          mu_u <- mu.G[,1] + V_prec %*% yss
	        e.draw <- as.numeric(mu_u + t(chol(V_u))%*%rnorm(M))#try(as.numeric(mu_u + t(chol(V_u + diag(M)*1e-6))%*%rnorm(M)), silent=TRUE)
	        if (is(e.draw, "try-error")) e.draw <- as.numeric(MASS::mvrnorm(1, mu_u, V_u + diag(M)*1e-4, tol=1e-5))

          epsilon[t, ] <- e.draw
          w[t, ] <- solve(Q.G[,,1])%*%epsilon[t, ] #maps back into the standard normal terms
        }
      }
    }
    
    #Sample the VAR coefficients conditional on w 
    yhat <- (Y - epsilon)

    for (i in seq_len(M)){
      normalizer <- 1/exp(.5*H[,i])
      y.i <- yhat[,i] * normalizer
      X.i <- X * normalizer
      
      V.a <- solve(crossprod(X.i) + diag(1/Theta[i,]))
      m.a <- V.a %*% (crossprod(X.i, y.i) + diag(1/Theta[i,])%*%A.prior[,i])
      a.draw <- try(m.a + t(chol(V.a))%*%rnorm(K), silent=TRUE)
      
      if (is(a.draw, "try-error")) a.draw <- m.a + t(chol(V.a, pivot=TRUE))%*%rnorm(K)
      
      A[ ,i] <- a.draw
      
      #Sample prior variances on the VAR coefficients
      if (!flat){
        prior.A <- get.ng(A[,i]-A.prior[,i], Theta[,i],a_tau=a_tau[i,1],acc=count_a[i,1], pscale=scale_a[i,1], sample_a=TRUE,irep=irep,nburn=nburn)
        a_tau[i,1] <- prior.A$a_tau
        count_a[i,1] <- prior.A$acc
        scale_a[i,1] <- prior.A$pscale
        
        Theta[i , ] <- t(prior.A$psi)
      }
    }

    Theta[Theta < 1e-5] <- 1e-5
    Theta[Theta > 3] <- 3
    
    #Sample the idiosyncratic error variances 
    v <- yhat - X%*%A
    
    if (!SV.idi){
      if (!NG.idi){
        for (i in seq_len(M)){
          a.i <- T/2 + 0.01
          b.i <- sum(v[,i]^2)/2 + 0.01
          Omega[i,i] <- 1/rgamma(1, a.i, b.i)
        }
      }else{
        for (i in seq_len(M)){
          q.i <- 0 #implies a G(1/2, 1/2) prior
          nu.i <- Omega[i,i]
          kn.i <- 2*1/2
          B.j[[i]] <- GIGrvg::rgig(1, lambda=q.i, chi=nu.i, psi=kn.i)
          
          
          a.i <- T/2 + 1/2
          b.i <- sum(v[,i]^2)
          c.i <- 1/(2*B.j[[i]])
          
          Omega[i,i] <- GIGrvg::rgig(1, lambda = a.i, chi = b.i, psi = c.i)
        }
      }
      H <- log(t(matrix(rep(diag(Omega), T), M, T)))
    } else {
      
      sv.vola <- matrix(NA, M, 1)
      for (i in seq_len(M)){
        res <- svsample_fast_cpp(v[,i], startpara = para.list[[i]], startlatent = H[,i], priorspec = prior.list[[i]])
        para.list[[i]][c("mu","phi","sigma")] <- as.list(res$para[, c("mu", "phi", "sigma")])
        latent <- drop(res$latent)
        latent[exp(latent) < 1e-6] <- log(1e-6)
        H[,i] <- latent
        
        sv.vola[i, 1] <- para.list[[i]]["sigma"]$sigma
      }

      # Draw the shrinkage factors for the log-volatility equation
      switch <- sample(c(-1,1), M, replace=TRUE)
      ng.vola <- get.ng(switch*sv.vola, C.mat,a_tau=a_vola,acc=count_vola, pscale=scale_vola, sample_a=TRUE,irep=irep,nburn=nburn)
      C.mat <- ng.vola$psi
      a_vola <- ng.vola$a_tau
      count_vola <- ng.vola$acc
      scale.vola <- ng.vola$pscale
      
      for (i in seq_len(M)) prior.list[[i]]$sigma2$rate <- 1/(2*C.mat[[i]]) 
    }

    #Start storing stuff after the burn-in stage
    if (irep > nburn){
      A.store[irep-nburn,,] <- A
      Sigma.store[irep-nburn,,,] <- Sigma.G 
      mu.store[irep-nburn,,] <- mu.G
      if (!SV.idi) Omega.store[irep-nburn,] <- diag(Omega) else H.store[irep-nburn,,] <- exp(.5*H)
      G.store[irep-nburn,] <- c(G, G.plus)
      S.store[irep-nburn,] <- t(S)
      count.store[irep-nburn,] <- t(counts)
      if (G.max > 1) eta.store[irep-nburn,] <- as.numeric(eta)
      OT.store[irep-nburn,] <- exp(H[T,])
      if (SV.idi) vola.vola.store[irep-nburn, ] <- switch*sv.vola else vola.vola.store <- NULL
      
      fit.store[irep-nburn,,] <- X%*%A + shock.epsilon + exp(.5*H)*matrix(rnorm(T * M), T, M)# %*% diag(sqrt(diag(Omega)))
      if (G.max > 1) for (j in seq_len(M)) vola.store[irep-nburn,,j] <- Q.t[,j,j]^2/(Q.t[,j,j]^2 + Omega[j,j])
      
      #Forecasts go here
      if (fhorz > 0){
        #Create draw from the predictive density
        y.pred <- matrix(NA, M, fhorz)
        XTpi <- c(Y[T,], X[T, 1:(M*(p-1))],1)
        for (ihorz in seq_len(fhorz)){
          
          if (G.max > 1) Sth <- sample(1:G.max, 1, prob=eta) else Sth <- 1
          
          YTpi <- as.vector(XTpi%*%A) + mu.G[, Sth]+Q.G[,,Sth]%*%rnorm(M) + diag(exp(.5*H[T,]))%*%rnorm(M) # epsilon(T+1) ~ N(0, Sigma), Sigma = Nprime Nprime'
          XTpi <- c(YTpi, XTpi[1:(M*(p-1))], 1)
          y.pred[, ihorz] <- YTpi
        }
        fcst.store[irep-nburn,,] <- y.pred
      }
      
      #IRFs go here
      if (irf.hor > 0){
        #Compute IRFs
        PHI_array <- array(0,c(M,M,p))
        for (ss in 1:p){
          PHI_array[,,ss] <- t(A[((ss-1)*M+1):(ss*M),])
        }
        
        for (ss in seq_len(G.max)){
          if (counts[ss,1]>0){
            irf.s <- impulsdtrf(PHI_array, Q.G[,,ss], irf.hor)
            IRF.store[irep-nburn,,,,ss] <- irf.s
          }
        }
      }
    }
    
    if (irep %% 50 == 0){
      print(paste0("MCMC sampling started. Currently, ", round(irep/ntot*100,digits=2), "% finished"))
    }
  }
  A.median <- apply(A.store, c(2,3), median)
  
  #Sort regimes by number of observations to solve label switching ex post
  S.re.ordered <- S.store
  eta.re.ordered <- eta.store
  
  for (irep in seq_len(nsave)){
    ind.sort <- sort(count.store[irep, ], decreasing = TRUE, index.return=T)$ix
    if (irf.hor > 1) IRF.store[irep,,,, ] <- IRF.store[irep,,,, ind.sort]
    Sigma.store[irep,,,] <- Sigma.store[irep,,,ind.sort]
    
    S.tt <- matrix(0, T, length(ind.sort))
    for (jj in 1:length(ind.sort)){
      S.tt[,jj] <- S.store[irep,]==jj*1
    }
    
    S.tt <- S.tt[,ind.sort]
    S.tt.re <- apply(t(t(S.tt)*1:length(ind.sort)),1,sum)
    S.re.ordered[irep, ] <- S.tt.re
    eta.re.ordered[irep,] <- eta.re.ordered[irep, ind.sort]
  }
  
  ret.list <- list("A_coefs"=A.store, "mu_G"=mu.store, "Sigma_G"=Sigma.store, "G"=G.store, "S"=S.re.ordered, "IRF"=IRF.store, "p(y(t+h)| Rest)"= fcst.store, "fit"=fit.store, "vola"=vola.store, "eta"=eta.re.ordered, "Omega.T"=OT.store, "C0"=C_0, "vola.vola.store"=vola.vola.store)
  return(ret.list)
}

# ------------------------------------------------------------------------------------
# MCMC setup
flexBART <- function(Yraw, nburn=5000, nsave = 5000, thinfac = 1, prior="HS", prior.sig, model ="mixBART", sv = "SV", fc.approx="exact", restr.var = NULL, fhorz = 3, quiet = FALSE, p=2){
  require(MASS)
  require(dbarts)
  require(stochvol)
  #Yraw <- Xraw; nburn=50; nsave = 50; thinfac = 1; prior="HS"; prior.sig=NULL; model ="BVAR"; sv = "SV"; fc.approx="exact"; restr.var = NULL; fhorz = 12; quiet = FALSE; p=2
  
  # auxiliary functions
  #source("aux_func.R")
  Z.sv <- NULL
  
  # standardize data
  Ymu <- apply(Yraw, 2, mean,na.rm=T)
  Ysd <- apply(Yraw, 2, sd,na.rm=T)
  Yraw <- apply(Yraw, 2, function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})
  
  iter.update <- 250
  
  nthin <- round(thinfac * nsave)
  ntot <- nburn + nsave
  thin.set <- floor(seq(nburn+1,ntot,length.out=nthin))
  in.thin <- 0
  
  cons <- FALSE # whether constant should be included
  id.L <- FALSE # whether factor model for the errors should be identified
  
  # create design matrices Y/X
  if(cons){
    X <- cbind(mlag(Yraw,p),1)[(p+1):nrow(Yraw),]
  }else{
    X <- cbind(mlag(Yraw,p))[(p+1):nrow(Yraw),]
  }
  Y <- Yraw[(p+1):nrow(Yraw), ] # leave original input matrix unchanged
  
  T <- nrow(Y)
  M <- ncol(Y)
  K <- ncol(X)
  
  # initialize draws
  PHI_draw  <- A_draw <- A_approx <- A_OLS <- ginv(crossprod(X))%*%crossprod(X,Y)
  A_prior <- matrix(0, K, M)
  Sig_OLS <- crossprod(Y-X%*%A_OLS)/T
  
  # factor structure in the errors
  eta <- E.hat <- Y - X%*%PHI_draw
  Q <- factorstochvol::ledermann(M)
  Lambda <- matrix(0,M,Q)
  Ft <- prcomp(eta)$x[,1:Q]
  Omega <- matrix(1,T,Q)
  theta.Lambda <- Lambda^0
  if(id.L){
    id.Lambda <- rbind(lower.tri(diag(Q)),matrix(TRUE,M-Q,Q))
    Lambda[1:Q,1:Q] <- diag(Q)
  }else{
    id.Lambda <- matrix(TRUE,M,Q)
  }
  ASV_approx <- matrix(0,K+1,M)
  AFSV_approx <- matrix(0,K+1,Q)
  
  # HS prior for factor loadings (by column)
  lambda.Lambda <- matrix(1,M,Q)
  nu.Lambda <- matrix(1,M,Q)
  zeta.Lambda <- tau.Lambda <- rep(1,Q)
  
  # covariance related objects
  eta <- matrix(NA,T,M)
  H <- matrix(-3,T,M)
  
  # objects for SV
  sv_draw <- list()
  sv_latent <- list()
  for (mm in seq_len(M)){
    sv_draw[[mm]] <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    sv_latent[[mm]] <- rep(0,T)
  }
  
  svfac_draw <- list()
  svfac_latent <- list()
  for (qq in seq_len(Q)){
    svfac_draw[[qq]] <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    svfac_latent[[qq]] <- rep(0,T)
  }
  
  sv_priors <- list()
  if(sv == "SV"){
    for(mm in 1:M){
      sv_priors[[mm]] <- specify_priors(
        mu = sv_normal(mean = 0, sd = 1),
        phi = sv_beta(shape1 = 5, shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*0.01)),
        nu = sv_infinity(),
        rho = sv_constant(0)
      ) 
    }
  }
  
  svfac_priors <- list()
  if(sv == "SV"){
    for(qq in 1:Q){
      svfac_priors[[qq]] <- specify_priors(
        mu = sv_normal(mean = 0, sd = 1),
        phi = sv_beta(shape1 = 5, shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*0.01)),
        nu = sv_infinity(),
        rho = sv_constant(0)
      ) 
    }
  }
  
  # priors
  theta_A <- matrix(10,K,M)
  a0sig <- 3
  b0sig <- 0.3
  
  if(prior=="Minn"){
    # MH preliminaries
    accept <- 0
    accept2 <- 0
    scale1 <- .43
    scale2 <- .45
    
    shrink.1 <- 0.1
    shrink.2 <- 0.2
    
    ind <- matrix(0,M,p)
    for (i in 1:M){
      ind[i,] <- seq(i,K-cons, by=M)
    }
    
    sigma_sq  <- matrix(0,M,1) #vector which stores the residual variance
    for (i in 1:M){
      Ylag_i <- mlag(Yraw[,i],p)
      Ylag_i <- Ylag_i[(p+1):nrow(Ylag_i),]
      Y_i <- Yraw[(p+1):nrow(Yraw),i]
      alpha_i <- solve(t(Ylag_i)%*%Ylag_i)%*%t(Ylag_i)%*%Y_i
      sigma_sq[i,1] <- (1/(nrow(Yraw)-p))*t(Y_i-Ylag_i%*%alpha_i)%*%(Y_i-Ylag_i%*%alpha_i) 
    }
    
    V_i <- get_V(shrink.1,shrink.2,ind1=ind,sigma_sq1=sigma_sq,p1=p,set.cons=cons,M1=M,K1=K)
    theta_PHI <- theta_A <- V_i
  }else if (prior=="NG"){
    #Prior for the VAR coefficients /initial values
    a_tau <- matrix(1, M, 1)
    count_a <- matrix(0, M, 1)
    scale_a <-  matrix(0.1, M, 1)
    
    theta_PHI <- theta_A <- matrix(10,K,M)
  }
  
  # BART initialization
  cgm.level <- 0.95
  cgm.exp <- 2
  sd.mu <- 2
  num.trees <- 250
  
  control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                           keepTrees = FALSE, n.samples = ntot,
                           n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
                           n.threads = 1, n.thin = 1L, printEvery = 1,
                           printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                           updateState = FALSE)
  
  sampler.list <- list()
  sampler.hetero.list <- list()
  sampler.fac.list <- list()
  svdraw.list <- list()
  
  for (jj in seq_len(M)){
    if (sv == "heteroBART" || sv == "SV"){ 
      prior.sig = c(10000^50, 0.5)
      sigma.init <- 1
    }else{
      sigma.init <- sqrt(Sig_OLS[jj,jj]) 
    }
    
    sampler.list[[jj]] <- dbarts(Y[,jj]~X, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=sigma.init, resid.prior = chisq(prior.sig[[1]], prior.sig[[2]]))
    if (sv == "heteroBART"){
      sampler.hetero.list[[jj]] <-  dbarts(E.hat[, jj] ~ cbind(seq(1, T), X), control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=3)
    }
  }
  
  for (qq in seq_len(Q)){
    if (sv == "heteroBART" || sv == "SV"){ 
      prior.sig = c(10000^50, 0.5)
      sigma.init <- 1
    }else{
      sigma.init <- 1
    }
    
    if (sv == "heteroBART"){
      sampler.fac.list[[qq]] <-  dbarts(Ft[,qq] ~ cbind(seq(1, T), X), control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=3)
    }
  }
  
  sampler.run <- list()
  sampler.sv.run <- list()
  sampler.fsv.run <- list()
  if (sv == "heteroBART"){
    # Stuff necessary for the 10-component mixture)
    m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
    v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
    q_st  = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)
  }
  
  sigma.fac.mat <- matrix(1,Q,1)
  sigma.mat <- matrix(1, M, 1)
  count.mat <- matrix(0, M*p, M)
  
  # initialize HS prior on linear VAR coefs
  lambda.VAR <- 1
  nu.VAR <- 1
  tau.VAR <- 1
  zeta.VAR <- 1
  sv_parms_mat <- matrix(NA, M, 4)
  svfac_parms_mat <- matrix(NA, Q, 4)
  Sig_t <- array(0,dim=c(T,M,M))
  
  Y.fit.BART <- Y.fit.VAR <- Y.fit.full <- Y*0
  
  # storage objects
  if (!is.null(restr.var)){
    sl.quants <- seq(0, 1, by=0.05)
    range.conditional <- quantile(Y[,restr.var], sl.quants)
    R <- length(range.conditional)
    H_store <- array(NA, dim=c(nthin, T, M)) #Store stochastic volatilities
    fcst_store <- Hfcst_store <- array(NA,dim=c(nthin,fhorz,M,R))
  }else{
    H_store <- array(NA, dim=c(nthin, T, M)) #Store stochastic volatilities
    fcst_store <- Hfcst_store <- array(NA,dim=c(nthin,fhorz,M))
    PHI_store <- array(NA, dim=c(nthin, K, M))
  }
  
  # -----------------------------------------------------------------------------
  # start Gibbs sampler
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = ntot, style = 3)
    start  <- Sys.time()
  }
  
  error.full <- eta - X%*%PHI_draw
  X.ginv <- MASS::ginv(X)
  F.ginv <- MASS::ginv(cbind(seq(1, T), X))
  count.dirt <- 0
  for (irep in 1:ntot){
    check <- TRUE # algorithm stability
    
    # 1) sample model coefficients (either linear VAR or BART)
    Y_ <- Y - tcrossprod(Ft,Lambda)
    while (check){
      for (mm in seq_len(M)){
        if (model == "mixBART" || model == "BART"){
          sampler.list[[mm]]$setResponse(Y_[,mm])
          
          # This part estimates the eq-specific BART model
          rep_mm <- sampler.list[[mm]]$run(0L, 1L)
          sampler.run[[mm]] <- rep_mm
          
          sig.draw <- rep_mm$sigma
          sigma.mat[mm,] <- sig.draw
          A_approx[,mm] <- X.ginv%*%rep_mm$train
          
          if (any(is.na(rep_mm$train))) break
          Y.fit.BART[,mm] <- rep_mm$train
          eta[,mm] <- Y_[,mm] - rep_mm$train
          A_draw[,mm] <- X.ginv%*%rep_mm$train
          count.mat[,mm] <- rep_mm$varcount
        }else{
          A_draw <- matrix(0, K, M)
          eta <- Y_
        }
        
        # Here we simulate the regression part of the model
        norm_mm <- as.numeric(exp(-.5*H[, mm]))
        eta_mm <- eta[,mm] * norm_mm
        
        if (model=="mixBART" || model=="BVAR"){
          #This part allows for estimating a VAR model with BART shocks that depend on the lagged values of Y as well
          u_mm <- X*norm_mm
          theta_PHI <- theta_A[,mm]
          v0.inv <- diag(1/c(theta_PHI)) 
          
          V.cov <- solve(crossprod(u_mm) + v0.inv)
          mu.cov <- V.cov %*% (crossprod(u_mm, eta_mm) + v0.inv%*%c(A_prior[, mm]))
          
          mu.phi.a0.draw <- try(mu.cov + t(chol(V.cov)) %*% rnorm(ncol(V.cov)), silent=TRUE)
          if (is(mu.phi.a0.draw, "try-error")){
            mu.phi.a0.draw <-  as.numeric(mvtnorm::rmvnorm(1, mu.cov, as.matrix(Matrix::forceSymmetric((V.cov)))))
          } 
          
          error.full[,mm] <- eta[,mm]
          PHI_draw[,mm] <- mu.phi.a0.draw[1:K]
        }else{
          PHI_draw <- matrix(0, K, M)
        }
      }
      
      max.eigen <- (max(abs(eigen(get_companion(PHI_draw, c(M, 0, p))$MM)$values)))
      if (irep > 0.05*nburn){
        if (max.eigen < 10^4) check <- FALSE
      }else{
        check <- FALSE
      }
    }
    Y.fit.VAR <- X%*%PHI_draw
    Y.fit.full <- Y.fit.BART + Y.fit.VAR
    
    shock_norm <- Y - Y.fit.full - tcrossprod(Ft,Lambda)
    for (mm in seq_len(M)){
      if (sv == "SV"){
        svdraw_mm <- svsample_general_cpp(shock_norm[,mm]/sigma.mat[mm], startpara = sv_draw[[mm]], startlatent = sv_latent[[mm]], priorspec = sv_priors[[mm]])
        sv_draw[[mm]][c("mu", "phi", "sigma")] <- as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
        sv_latent[[mm]] <- svdraw_mm$latent
        sv_parms_mat[mm, ] <- c(svdraw_mm$para[, c("mu", "phi", "sigma")], svdraw_mm$latent[T])
        
        normalizer <- as.numeric(exp(-.5*svdraw_mm$latent))
        weights.new <- as.numeric(exp(-svdraw_mm$latent))
        sampler.list[[mm]]$setWeights(weights.new)
        H[,mm] <- log(sigma.mat[mm]^2) + svdraw_mm$latent
      }else if (sv =="heteroBART"){
        eta.star <- shock_norm[ , mm]/sigma.mat[mm]
        eta.off <- 1e-15 
        eta.star <- log((eta.star - eta.off)^2)
        
        # Sample the mixture components (previous draw of ht)
        z <- sapply(eta.star-H[,mm], ncind, m_st, sqrt(v_st2), q_st)
        
        # Subset mean and variances to the sampled mixture components; (n x p) matrices
        m_st_all = matrix(m_st[z],T, 1)
        v_st2_all = matrix(v_st2[z], T, 1)
        eta.star <- matrix((eta.star - m_st_all))
        
        # Do BART for conditional variances
        sampler.hetero.list[[mm]]$setResponse(eta.star)
        sampler.hetero.list[[mm]]$setWeights(1/v_st2_all)
        
        rep.hetero.i <- sampler.hetero.list[[mm]]$run(0L, 1L)
        sampler.sv.run[[mm]] <- rep.hetero.i
        
        ht.i <- rep.hetero.i$train
        ht.i[ht.i < -10] <- -10
        ASV_approx[,mm] <- F.ginv%*%ht.i
        
        weights.new <- as.numeric(exp(-ht.i))
        sampler.list[[mm]]$setWeights(weights.new)
        H[ , mm] <- ht.i + log(as.numeric(sigma.mat[mm])^2)
      }else{
        H[,mm] <- log(sigma.mat[mm]^2)
        if(model=="BVAR"){
          H[,mm] <- log(1/rgamma(1,0.01+T/2,0.01+sum(shock_norm[,mm]^2)))
        }
      }
    }
    H[H<log(1e-6)] <- log(1e-6)
    
    # draw factors and loadings
    eps <- Y - Y.fit.full
    Ft <- get.factors(eps,S=exp(H),H=exp(Omega),L=Lambda,q=Q,t=T)
    if(!id.L) Ft <- apply(Ft,2,function(x) (x-mean(x))/sd(x)) # normalize factor draw
    Lambda <- get.Lambda(eps,fac=Ft,S=exp(H),pr=theta.Lambda,m=M,q=Q,id.fac=!id.Lambda[1,Q])
    
    # shrinkage on loadings by columns
    for(qq in 1:Q){
      hs_draw <- get.hs(bdraw=as.numeric(Lambda[id.Lambda[,qq],qq]),
                        lambda.hs=as.numeric(lambda.Lambda[id.Lambda[,qq],qq]),
                        nu.hs=nu.Lambda[id.Lambda[,qq],qq],tau.hs=tau.Lambda[qq],zeta.hs=zeta.Lambda[qq])
      theta.Lambda[id.Lambda[,qq],qq] <- hs_draw$psi
      lambda.Lambda[id.Lambda[,qq],qq] <- hs_draw$lambda
      nu.Lambda[id.Lambda[,qq],qq] <- hs_draw$nu
      tau.Lambda[qq] <- hs_draw$tau
      zeta.Lambda[qq] <- hs_draw$zeta
    }
    theta.Lambda[theta.Lambda<1e-5] <- 1e-5
    
    # draw variances of factors
    if(sv=="SV"){
      for(qq in 1:Q){
        svdraw_qq <- svsample_general_cpp(Ft[,qq]/sigma.fac.mat[qq], startpara = svfac_draw[[qq]], 
                                          startlatent = svfac_latent[[qq]], priorspec = svfac_priors[[qq]])
        svfac_draw[[qq]][c("mu", "phi", "sigma")] <- as.list(svdraw_qq$para[, c("mu", "phi", "sigma")])
        svfac_latent[[qq]] <- svdraw_qq$latent
        svfac_parms_mat[qq, ] <- c(svdraw_qq$para[, c("mu", "phi", "sigma")], svdraw_qq$latent[T])
        Omega[,qq] <- log(sigma.fac.mat[qq]^2) + svdraw_qq$latent
      }
    }else if(sv=="heteroBART"){
      for(qq in 1:Q){
        eta.star <- Ft[,qq]/sigma.fac.mat[qq]
        eta.off <- 1e-15 
        eta.star <- log((eta.star - eta.off)^2)
        
        # Sample the mixture components (previous draw of ht)
        z <- sapply(eta.star-Omega[,qq],ncind,m_st,sqrt(v_st2),q_st)
        
        # Subset mean and variances to the sampled mixture components; (n x p) matrices
        m_st_all = matrix(m_st[z],T,1)
        v_st2_all = matrix(v_st2[z],T,1)
        eta.star <- matrix((eta.star - m_st_all))
        
        # Do BART for conditional variances
        sampler.fac.list[[qq]]$setResponse(eta.star)
        sampler.fac.list[[qq]]$setWeights(1/v_st2_all)
        
        rep.fac.q <- sampler.fac.list[[qq]]$run(0L, 1L)
        sampler.fsv.run[[qq]] <- rep.fac.q
        
        omega.q <- rep.fac.q$train
        omega.q[omega.q < -8] <- -8
        AFSV_approx[,qq] <- F.ginv%*%omega.q
        Omega[,qq] <- omega.q + log(as.numeric(sigma.fac.mat[qq])^2)
      }
    }else{
      for(qq in 1:Q){
        Omega[,qq] <- 0 # log(1/rgamma(1,3+T/2,0.3+sum(Ft[,qq]^2)))
      }
    }
    
    for(tt in seq_len(T)){
      if (Q  > 1 ){
        S.t <- Lambda%*%tcrossprod(diag(exp(Omega[tt,])),Lambda) + diag(exp(H[tt,]))
      }else{
        S.t <- exp(Omega[tt,]) * (Lambda %*% t(Lambda)) + diag(exp(H[tt,]))
      }
      Sig_t[tt,,] <- S.t
    }
    
    # 2) updating shrinkage priors
    if (model=="mixBART" || model=="BVAR"){
      #Update HS for the VAR coefficients
      if(prior=="HS"){
        hs_draw <- get.hs(bdraw=as.numeric(PHI_draw - A_prior),lambda.hs=lambda.VAR,nu.hs=nu.VAR,tau.hs=tau.VAR,zeta.hs=zeta.VAR)
        lambda.VAR <- hs_draw$lambda
        nu.VAR <- hs_draw$nu
        tau.VAR <- hs_draw$tau
        zeta.VAR <- hs_draw$zeta
        theta_A <- matrix(hs_draw$psi,K,M)
      }else if(prior=="Minn"){
        shrink.1.prop <- exp(rnorm(1,0,scale1))*shrink.1
        theta.prop <- get_V(shrink.1.prop,shrink.2,ind1=ind,sigma_sq1=sigma_sq,p1=p,set.cons=cons,M1=M,K1=K)
        
        post1.prop <- sum(dnorm(as.vector(PHI_draw),as.vector(A_prior),sqrt(as.vector(theta.prop)),log=TRUE))+dgamma(shrink.1.prop,0.01,0.01,log=TRUE)
        post1 <- sum(dnorm(as.vector(PHI_draw),as.vector(A_prior),sqrt(as.vector(theta_A)),log=TRUE))+dgamma(shrink.1,.01,0.01,log=TRUE)
        if ((post1.prop-post1)>log(runif(1,0,1))){
          shrink.1 <- shrink.1.prop
          theta_A <- theta.prop
          accept <- accept+1
        }
        
        shrink.2.prop <- exp(rnorm(1,0,scale2))*shrink.2
        theta.prop <- get_V(shrink.1,shrink.2.prop,ind1=ind,sigma_sq1=sigma_sq,p1=p,set.cons=cons,M1=M,K1=K)
        
        post2.prop <- sum(dnorm(as.vector(PHI_draw),as.vector(A_prior),sqrt(as.vector(theta.prop)),log=TRUE))+dgamma(shrink.2.prop,0.01,0.01,log=TRUE)
        post2 <- sum(dnorm(as.vector(PHI_draw),as.vector(A_prior),sqrt(as.vector(theta_A)),log=TRUE))+dgamma(shrink.2,0.01,0.01,log=TRUE)
        if ((post2.prop-post2)>log(runif(1,0,1))){
          shrink.2 <- shrink.2.prop
          theta_A <- theta.prop
          accept2 <- accept2+1
        }
        
        if (irep<(0.5*nburn)){
          if ((accept/irep)>0.3) scale1 <- 1.01*scale1
          if ((accept/irep)<0.15) scale1 <- 0.99*scale1
          if ((accept2/irep)>0.3) scale2 <- 1.01*scale2
          if ((accept2/irep)<0.15) scale2 <- 0.99*scale2
        }
      }else if (prior =="NG"){
        for (i in 1:M){
          prior.A <- get.ng(PHI_draw[,i], theta_A[,i],a_tau=a_tau[i,1],acc=count_a[i,1], pscale=scale_a[i,1], sample_a=TRUE,irep=irep,nburn=nburn)
          a_tau[i,1] <- prior.A$a_tau
          count_a[i,1] <- prior.A$acc
          scale_a[i,1] <- prior.A$pscale
          theta_A[,i] <- t(prior.A$psi)
        }
      }
      theta_A[theta_A < 1e-12] <- 1e-12
    }
    
    if(irep %in% thin.set){
      in.thin <- in.thin+1
      H_store[in.thin,,] <- H
      PHI_store[in.thin,,] <- PHI_draw
      
      # include forecast loop here
      if(fhorz>0 & is.null(restr.var)){
        Yfc <- matrix(NA,fhorz,M)
        Hfc <- matrix(NA,fhorz,M)
        
        if (cons){
          X.hat <- c(Y[T,],X[T,1:(M*(p-1))],1)
        }else{
          X.hat <- c(Y[T,],X[T,1:(M*(p-1))])  
        }
        
        if (sv == "SV" || sv == "heteroBART"){
          HT <-  H[T, ] - log(as.numeric(sigma.mat)^2)
          OT <- Omega[T, ] - log(as.numeric(sigma.fac.mat)^2)
        }else{
          HT <- H[T,]
          OT <- Omega[T, ]
        }
        
        Sig_T <- Sig_t[T,,] # use final observation for Sigma
        tree.pred <- matrix(0, M)
        for (hh in seq_len(fhorz)){
          if (sv == "SV"){
            HT  <- log(as.numeric(sigma.mat)^2) + (sv_parms_mat[, 1] + sv_parms_mat[ , 2] * (HT - sv_parms_mat[,1]) + sv_parms_mat[ , 3]*rnorm(M))
            OT <- log(as.numeric(sigma.fac.mat)^2) + (svfac_parms_mat[,1] + svfac_parms_mat[,2] * (OT - svfac_parms_mat[,1]) + svfac_parms_mat[,3]*rnorm(Q))
          }else if (sv == "heteroBART"){
            vola.predict.tree <- matrix(0, M, 1)
            for (nn in seq_len(M)){
              vola.predict.tree[nn,] <- sampler.hetero.list[[nn]]$predict(c(T+hh, X.hat))#
            }
            HT <- as.numeric(log(as.numeric(sigma.mat)^2) + vola.predict.tree)
            
            vola.predict.fac <- matrix(0, Q, 1)
            for (qq in seq_len(Q)){
              vola.predict.fac[qq,] <- sampler.fac.list[[qq]]$predict(c(T+hh,X.hat))
            }
            OT <- as.numeric(log(as.numeric(sigma.fac.mat)^2) + vola.predict.fac)
          }
          Hfc[hh,] <- exp(HT)
          
          if (model == "mixBART" || model == "BART"){
            for (j in seq_len(M)) tree.pred[j] <- sampler.list[[j]]$predict(X.hat)
          }else{
            tree.pred <- rep(0, M)
          }
          if (Q>1) Sig_T <- Lambda%*%tcrossprod(diag(exp(OT)),Lambda) + diag(exp(HT)) else Sig_T <- exp(OT) * (Lambda%*%t(Lambda)) + diag(exp(HT))
          
          if(fc.approx=="exact"){
            Y.tp1 <- try(as.numeric(X.hat%*%PHI_draw)+ as.numeric(tree.pred) + t(chol(Sig_T))%*%rnorm(M), silent=TRUE)
            if (is(Y.tp1, "try-error")) Y.tp1 <- as.numeric(mvtnorm::rmvnorm(1, tree.pred + as.numeric(X.hat%*%PHI_draw), Sig_T))
          }else if(fc.approx=="approx"){
            Y.tp1 <- try(as.numeric(X.hat%*%(A_approx+PHI_draw)) + t(chol(Sig_T))%*%rnorm(M), silent=TRUE)
            if (is(Y.tp1, "try-error")) Y.tp1 <- as.numeric(mvtnorm::rmvnorm(1, as.numeric(X.hat%*%(A_approx+PHI_draw)), Sig_T))
          }
          
          if (cons){
            X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))],1)
          }else{
            X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))])
          }
          Yfc[hh,] <- Y.tp1
        }
        fcst_store[in.thin,,] <- (Yfc*t(matrix(Ysd,M,fhorz)))+t(matrix(Ymu,M,fhorz))
        Hfcst_store[in.thin,,] <- Hfc*t(matrix(Ysd,M,fhorz))
      }else if(fhorz>0 & !is.null(restr.var)){
        Yfc <- array(NA,c(fhorz,M, R))
        Hfc <- array(NA,c(fhorz,M, R))
        sl.cond <- which(colnames(Y)==restr.var)
        
        for (r in seq_len(R)){
          cond.vals <- range.conditional[[r]]
          sl.cond.lags <- seq(sl.cond, K, by=M)
          if (cons){
            X.hat <- c(Y[T,],X[T,1:(M*(p-1))],1)
          }else{
            X.hat <- c(Y[T,],X[T,1:(M*(p-1))])  
          }
          X.hat[sl.cond.lags] <- cond.vals
          
          if (sv == "SV" || sv == "heteroBART"){
            HT <-  H[T, ] - log(as.numeric(sigma.mat)^2)
            OT <- Omega[T, ] - log(as.numeric(sigma.fac.mat)^2)
          }else{
            HT <- H[T,]
            OT <- Omega[T, ]
          }
          
          Sig_T <- Sig_t[T,,] # use final observation for Sigma
          tree.pred <- matrix(0, M)
          for (hh in seq_len(fhorz)){
            if (sv == "SV"){
              HT  <- log(as.numeric(sigma.mat)^2) + (sv_parms_mat[, 1] + sv_parms_mat[ , 2] * (HT - sv_parms_mat[,1]) + sv_parms_mat[ , 3]*rnorm(M))
              OT <- log(as.numeric(sigma.fac.mat)^2) + (svfac_parms_mat[,1] + svfac_parms_mat[,2] * (OT - svfac_parms_mat[,1]) + svfac_parms_mat[,3]*rnorm(Q))
            }else if (sv == "heteroBART"){
              vola.predict.tree <- matrix(0, M, 1)
              for (nn in seq_len(M)){
                vola.predict.tree[nn,] <- sampler.hetero.list[[nn]]$predict(c(T+hh, X.hat))
              }
              HT <- as.numeric(log(as.numeric(sigma.mat)^2) + vola.predict.tree)
              
              vola.predict.fac <- matrix(0, Q, 1)
              for (qq in seq_len(Q)){
                vola.predict.fac[qq,] <- sampler.fac.list[[qq]]$predict(c(T+hh,X.hat))
              }
              OT <- as.numeric(log(as.numeric(sigma.fac.mat)^2) + vola.predict.fac)
            }
            Hfc[hh,,r] <- exp(HT)*Ysd
            
            if (model == "mixBART" || model == "BART"){
              for (j in seq_len(M)) tree.pred[j] <- sampler.list[[j]]$predict(X.hat)
              tree.pred[sl.cond] <- cond.vals
            }else{
              tree.pred <- rep(0, M)
            }
            if (Q > 1 ) Sig_T <- Lambda%*%tcrossprod(diag(exp(OT)),Lambda) + diag(exp(HT)) else Sig_T <- exp(OT) * (Lambda %*% t(Lambda)) + diag(exp(HT))
            Y.tp1 <- try(as.numeric(X.hat%*%PHI_draw)+ as.numeric(tree.pred) + t(chol(Sig_T))%*%rnorm(M), silent=TRUE)
            if (is(Y.tp1, "try-error")) Y.tp1 <- as.numeric(mvtnorm::rmvnorm(1, tree.pred + as.numeric(X.hat%*%PHI_draw), Sig_T))
            Y.tp1[sl.cond] <- cond.vals
            if (cons){
              X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))],1)
            }else{
              X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))])
            }
            Yfc[hh,,r] <- Y.tp1*Ysd + Ymu
          }
        }
        fcst_store[in.thin,,,] <- Yfc
        Hfcst_store[in.thin,,,] <- Hfc
      }
    }
    
    if(!quiet){
      setTxtProgressBar(pb, irep)
      if (irep %% iter.update==0){
        end <- Sys.time()
        message(paste0("\n Average time for single draw over last ",iter.update," draws ", round(as.numeric(end-start)/iter.update, digits=4), " seconds, currently at draw ", irep))
        start <- Sys.time()
        ts.plot(cbind(Y.fit.full, Y), col=c(rep(1, M), rep(2, M)))
      }
    }
  }
  dimnames(fcst_store) <- dimnames(Hfcst_store) <- list(paste0("mcmc",1:nthin),paste0("fhorz", 1:fhorz),colnames(Y))
  return_obj <- list("fcst"=fcst_store,"Hfcst"=Hfcst_store, "PHI"=PHI_store)
  return(return_obj)
}


impulsdtrf <- function(B,smat,nstep)
  ### By:             As emerges from rfvar, neqn x nvar x lags array of rf VAR coefficients.
  ### smat:           nshock x nvar matrix of initial shock vectors.  To produce "orthogonalized
  ###                 impulse responses" it should have the property that crossprod(t(smat))=sigma,
  ###                 where sigma is the Var(u(t)) matrix and u(t) is the rf residual vector.  One
  ###                 way to get such a smat is to set smat=t(chol(sigma)).  To get the smat
  ###                 corresponding to a different ordering, use
  ###                 smat = t(chol(P %*% Sigma %*% t(P)) %*% P), where P is a permutation matrix.
  ###                 To get impulse responses for a structural VAR in the form A(L)y=eps, with
  ###                 Var(eps)=I, use B(L)=-A_0^(-1)A_+(L) (where A_+ is the coefficients on strictly
  ###                 positive powers of L in A), smat=A_0^(-1).
  ###                 In general, though, it is not required that smat be invertible.
### response:       nvar x nshocks x nstep array of impulse responses.
###
### Code written by Christopher Sims, based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.
{
  ##-----debug--------
  ##browser()
  ##------------------
  neq <- dim(B)[1]
  nvar <- dim(B)[2]
  lags <- dim(B)[3]
  dimnB <- dimnames(B)
  if(dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
  response <- array(0,dim=c(neq,nvar,nstep+lags-1));
  response[ , , lags] <- smat
  response <- aperm(response, c(1,3,2))
  irhs <- 1:(lags*nvar)
  ilhs <- lags * nvar + (1:nvar)
  response <- matrix(response, ncol=neq)
  B <- B[, , seq(from=lags, to=1, by=-1)]  #reverse time index to allow matrix mult instead of loop
  B <- matrix(B,nrow=nvar)
  for (it in 1:(nstep-1)) {
    #browser()
    response[ilhs, ] <- B %*% response[irhs, ]
    irhs <- irhs + nvar
    ilhs <- ilhs + nvar
  }
  ## for (it in 2:nstep)
  ##       {
  ##         for (ilag in 1:min(lags,it-1))
  ##           response[,,it] <- response[,,it]+B[,,ilag] %*% response[,,it-ilag]
  ##       }
  dim(response) <- c(nvar, nstep + lags - 1, nvar)
  response <- aperm(response[ , -(1:(lags-1)), ], c(1, 3, 2)) #drop the zero initial conditions; array in usual format
  dimnames(response) <- list(dimnB[[1]], dimnames(smat)[[2]], NULL)
  ## dimnames(response)[2] <- dimnames(smat)[1]
  ## dimnames(response)[1] <- dimnames(B)[2]
  return(response)
}

# get a draw for the normal-gamma prior
get.ng <- function(bdraw,psi,pr=0,e0=0.01,e1=0.01,sample_a,a_tau,acc=0,pscale=0.1,irep=1,nburn=1){
  k <- NROW(bdraw)
  if(pr==0){pr <- rep(0,k)}
  
  # global shrinkage parameter
  e0_po <- e0 + a_tau*k
  e1_po <- e1 + (a_tau * sum(psi))/2
  lambda <- rgamma(1,e0_po,e1_po)
  
  # local shrinkage scalings
  psi <- matrix(0,k,1)
  for(kk in 1:k){
    scale <- ifelse((bdraw[kk]-pr[kk])^2<=1e-12, 1e-12, (bdraw[kk]-pr[kk])^2) # for stability
    psi[kk] <- GIGrvg::rgig(n=1, lambda=a_tau-0.5, chi=scale, psi=lambda*a_tau)
  }
  
  if(sample_a){
    atau_c <- a_tau # current value
    latau_p <- log(atau_c) + rnorm(1,0,pscale) # proposal
    
    atau_p <- exp(latau_p)
    lpost_c <- sum(dgamma(as.vector(psi),atau_c,(atau_c*lambda/2),log=TRUE))+dexp(atau_c,rate=1,log=TRUE)  
    lpost_p <- sum(dgamma(as.vector(psi),atau_p,(atau_p*lambda/2),log=TRUE))+dexp(atau_p,rate=1,log=TRUE)  
    alpha1 <- min(1,exp(lpost_p-lpost_c+log(atau_p)-log(atau_c)))
    
    if(alpha1>runif(1)){
      a_tau <- atau_p
      acc <- acc+1
    }
    
    if (irep<(0.5*nburn)){
      if((acc/irep)>0.35) pscale <- 1.01*pscale
      if((acc/irep)<0.15) pscale <- 0.99*pscale
    }
  }
  
  return(list("psi"=psi,"a_tau"=a_tau,"acc"=acc,"pscale"=pscale))
}


#------ These helper functions come from some other project on BART
# ledermann bound (maximum number of factors)
ledermann <- function(m) {
  as.integer(floor((2*m+1)/2 - sqrt((2*m+1)^2/4 - m^2 + m)))
}

# functions for calculating KL divergence
my.ecdf  <-  function(x)   {
  x   <-   sort(x)
  x.u <-   unique(x)
  n  <-  length(x) 
  x.rle  <-  rle(x)$lengths
  y  <-  (cumsum(x.rle)-0.5) / n
  FUN  <-  approxfun(x.u, y, method="linear", yleft=0, yright=1,
                     rule=2)
  return(FUN)
}   

KL_est  <-  function(x, y)   {
  dx  <-  diff(sort(unique(x)))
  dy  <-  diff(sort(unique(y)))
  ex  <-  min(dx) ; ey  <-  min(dy)
  e   <-  min(ex, ey)/2
  n   <-  length(x)    
  P  <-   my.ecdf(x) ; Q  <-  my.ecdf(y)
  KL  <-  sum( log( (P(x)-P(x-e))/(Q(x)-Q(x-e)))) / n
  return(KL)              
}

# function for quantile weighted CRPS
qwCRPS_sample <-function(true,mcmc,tau,weighting="none"){
  require(pracma)
  
  tau_len <- length(tau)
  Q.tau <- quantile(mcmc,probs=tau)
  true_rep <- rep(true,tau_len)
  QS.vec <- (true_rep-Q.tau)*(tau-((true_rep<=Q.tau)*1))
  
  weights <- switch(tolower(weighting),
                    "none" = 1,
                    "tails" = (2*tau-1)^2,
                    "right" = tau^2,
                    "left" = (1-tau)^2,
                    "center" = tau*(1-tau))
  wghs <- QS.vec*weights
  return(pracma::trapz(tau,wghs))
}

# This file contains several auxiliary functions.
ncind <- function(y,mu,sig,q){
  sample(1:length(q),
         size = 1,
         prob = q*dnorm(y,mu,sig))
}

# get posteriors for the horseshoe prior (see Makalic & Schmidt, 2015)
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

# lag variables
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(NA,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

# companion matrix
get_companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0, (nl-1)*nn, nn)))
  
  return(list(MM=MM,Jm=Jm))
}
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 15 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# prior moments for Minnesota prior
get_V <- function(a_bar_1,a_bar_2,ind1,sigma_sq1,p1,set.cons,M1,K1){
  V_i <- matrix(0,K1,M1)
  #this double loop fills the prior covariance matrix
  for (i in 1:M1){ #for each i equation
    for (j in 1:K1){ #for each variable on the rhs
      if (set.cons==1) {
        if (j==K1) {
          V_i[j,i] <- 10^3*sigma_sq1[i,1] #variance on constant, trend, dummies and on ex variables         
        }    
        else if (any(j==ind1[i,])) {
          ll <- which(ind1[i,]==j)
          V_i[j,i] <- a_bar_1/(ll^2) #variance on own lags
        }else{
          ll <- which(ind1==j,arr.ind=TRUE)[2]
          kj <- ind1[which(ind1==j,arr.ind=TRUE)[1],1]
          V_i[j,i] <- (a_bar_2*sigma_sq1[i,1])/((ll^2)*sigma_sq1[kj,1])
        }
      }else{
        if (any(j==ind1[i,])) {
          ll <- which(ind1[i,]==j)
          V_i[j,i] <- a_bar_1/(ll^2) #variance on own lags
        } 
        else{
          ll <- which(ind1==j,arr.ind=TRUE)[2]
          kj <- ind1[which(ind1==j,arr.ind=TRUE)[1],1]
          V_i[j,i] <- (a_bar_2*sigma_sq1[i,1])/((ll^2)*sigma_sq1[kj,1])
        }
      }
    }
  }  
  return(V_i)
}

# -----------------------------------------------------------------------------------------------
# function to draw the factor loadings (basic linear regression)
get.facload <- function(yy,xx,l_sd){
  V_prinv <- diag(NCOL(xx))/l_sd
  V_lambda <- solve(crossprod(xx) + V_prinv)
  lambda_mean <- V_lambda %*% (crossprod(xx,yy))
  
  lambda_draw <- lambda_mean + t(chol(V_lambda)) %*% rnorm(NCOL(xx))
  return(lambda_draw)
}

# factor loadings draw
get.Lambda <- function(eps,fac,S,pr,m,q,id.fac){
  L <- matrix(0,m,q)
  if(id.fac){
    for(jj in 1:m){
      if (jj<=q){
        normalizer <- exp(0.5*S[,jj])
        yy0 <- (eps[,jj]-fac[,jj])/normalizer
        xx0 <- fac[,1:(jj-1),drop=FALSE]/normalizer
        if (jj>1){
          l_sd <- pr[jj,1:(jj-1)]
          lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
        }else{
          lambda0 <- 1
        }
        
        if (jj>1){
          L[jj,1:(jj-1)] <- lambda0
          L[jj,jj] <- 1
        }else if (jj==1){
          L[jj,jj] <- 1
        }
      }else{
        normalizer <- exp(0.5*S[,jj])
        yy0 <- (eps[,jj])/normalizer
        xx0 <- fac[,,drop=FALSE]/normalizer
        l_sd <- pr[jj,]
        lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
        L[jj,] <- lambda0
      }
    }
  }else{
    for(jj in 1:m){
      normalizer <- exp(0.5*S[,jj])
      yy0 <- (eps[,jj])/normalizer
      xx0 <- fac[,,drop=FALSE]/normalizer
      l_sd <- pr[jj,]
      lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
      L[jj,] <- lambda0
    }
  }
  return(L)
}

# sample the latent factors
get.factors <- function(e,S,H,L,q,t){
  F_raw <- matrix(0,t,q)
  for (tt in 1:t){
    normalizer <- exp(-S[tt,]/2)
    Lt <- L*normalizer
    yt <- e[tt,]*normalizer
    
    if (q==1) fac.varinv <-  1/H[tt] else fac.varinv <- diag(q)/H[tt]
    fac.Sigma <-  solve(crossprod(Lt)+fac.varinv)
    fac.mean <- fac.Sigma%*%crossprod(Lt,yt)
    
    F_temp <- try(fac.mean + t(chol(fac.Sigma)) %*% rnorm(q),silent=TRUE)
    if (is(F_temp,"try-error")) F_temp <- fac.mean + t(chol(fac.Sigma+diag(q)*1e-6)) %*% rnorm(q)
    F_raw[tt,] <- F_temp
  }
  return(F_raw)
}
