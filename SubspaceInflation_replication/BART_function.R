BART_bnp <- function(y, X, yho, Xho,nburn=5000, ntot=10000, thinning=1, PDP=FALSE, sample.omega=TRUE, tau2.start=0.5, G.max=25, fcst=TRUE, a0=1, a1=1, PCA=TRUE, mix.sv =FALSE){
#  y<- y.tilde; X<-X.new; yho <- NULL; Xho<- NULL; nburn=5000; ntot=10000; thinning=1; PDP=FALSE; sample.omega=TRUE; tau2.start=0.5; G.max=25; fcst=TRUE; a0=1; a1=1; PCA=TRUE
#  y=y; X=X; yho=yho; Xho=Xho;nburn=nburn; ntot=ntot; thinning=2; PDP=FALSE; sample.omega=TRUE; tau2.start=0.5; G.max=G.slct; fcst=TRUE; a0=1; a1=1; PCA=PCA
   nsave <- floor(ntot / thinning)

  
  #G.max <- "SV"
  require(Matrix)
  require(dbarts)
  n <- dim(X)[1] #Sample size
  p <- dim(X)[2] #Number of covariates
  
  a <- 1/2
  b <- 1/2
 

  if (PCA){
    #The Optimal Hard Threshold for Singular Values is 4/sqrt(3)
   qstar <- min(6, ncol(X)) #in case of very few columns, do not use dim reduction
   
   svd.x <- svd(X)
   Xhat <- svd.x$u
   Xhat <- Xhat[, 1:qstar, drop=FALSE]
  
   K <- ncol(Xhat)
   
  }else{
    Xhat <- X

    K <- ncol(Xhat)
  }
  Q0 <- Xhat%*%solve(crossprod(Xhat)+diag(ncol(Xhat))*1e-4)%*%t(Xhat)
  Xhat <- X
  sigma2.draw <- rep(1e-3, n)

  if (G.max =="SV"){
    require(stochvol)
    G.max <- 1
    SV.ident <- TRUE
    
    sv.draw <- list(mu = 0, phi = 0.6, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    sv.latent <- rep(0, n)
    
    sv_priors <- specify_priors(
      mu = sv_normal(mean = 0, sd = 10), # -4, 1e-3
      phi = sv_beta(shape1 = 25, shape2 = 5),
      sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*1)),
      nu = sv_infinity(),
      rho = sv_constant(0))
    
    startpara <- list(mu = 0, phi = 0.9, sigma = 0.1,
                      nu = Inf, rho = 0, beta = NA,
                      latent0 = 0)
  }else{
    SV.ident <- FALSE
    G.max <- G.max
  }
  G <- G.max 
  kappa <- 0.8 #truncation parameter for the slice sampler 
  xi <- (1-kappa)*kappa^(seq(1, G)-1)
  #--------------------------Prior hyperparameters --------------------------
  alpha <- 0.5 # initialize at expected value
  a_k <- 1
  b_k <- alpha
  
  #Prior on alpha (standard choice stipulated in Escobar & West)
  a_alpha <- 2
  b_alpha <- 4
  
  #Prior on the the component means and variances
  mu.0 <- 0
  
  if (G==1) zeta.0 <- 10^5 else zeta.0 <- 1
  
  
  g.0 <- 10
  g.1 <- 5
  
  #-----------------------------Starting values for the infinite mixture -----------------------------
  if (G.max > 1){
    #Classification
    start.flexmix <- flexmix::flexmix((y-Q0%*%y) ~ 1, k=3)
    si <- start.flexmix@posterior$scaled
    S <- apply(si, 1, function(x) which(x==max(x)))#apply(t(t(si) * seq(1, ncol(si))), 1, sum)#rep(1, n)
    
    mu.G <- rep(0, G)
    sd.G <- rep(1, G)
    
    for (jj in seq_len(ncol(si))){
      mu.G[[jj]] <- parameters(start.flexmix)[1,jj]
      sd.G[[jj]] <- parameters(start.flexmix)[2,jj]
    }
    
    counts <- matrix(0, G, 1)
    S.tab <- table(S)
    counts[as.numeric(names(S.tab))] <- S.tab
    Nclust <- length(counts)	  # initial number of clusters
    
    mu.t <- rep(0,n)
    sd.t <- rep(0,n)
    
    for (t in seq_len(n)){
      s.i <- S[[t]]
      mu.t[[t]] <- mu.G[[s.i]]
      sd.t[[t]] <- sd.G[[s.i]]
    }
  }else{
    mu.t <- rep(0,n)
    sd.t <- rep(1e-3,n)
    mu.G <- rep(0, G)
    sd.G <- rep(1, G)
    
    counts <- matrix(0, 1, 1)
    counts[1,1] <- n
    S <- matrix(1, n, 1)
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
  
  #Conditional mean modeling through BART
  cgm.level <- 0.95
  cgm.exp   <- 2
  sd.mu     <- 2
  num.trees <- 250
  
  prior.sig <- c(10000^50, 0.5)
  sigma.init <- 1
  
  control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                           keepTrees = FALSE, n.samples = ntot,
                           n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
                           n.threads = 1, n.thin = 1L, printEvery = 1,
                           printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                           updateState = FALSE)
  
  sampler.list <- dbarts(y~X, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,n), sigma=sigma.init, resid.prior = chisq(prior.sig[[1]], prior.sig[[2]]))
  

  
  #Storage matrices
  S.store <- matrix(NA, nsave, n)
  mu.store <- matrix(NA, nsave, n)
  sd.store <- matrix(NA, nsave, n)
  yhat.store <- shock.store <- fit.store <- matrix(NA, nsave, n)
  alpha.store <- matrix(NA, nsave, 1)
  G.store <- matrix(NA, nsave, 1)
  omega.store <- matrix(NA, nsave, 1)
  theta.store <- matrix(NA, nsave, 2)
  if (fcst){
    pred.store <- array(NA, c(nsave, 1))
    fevd.store <- array(NA, c(nsave, 1))
    lpl.store <- array(NA, c(nsave, 1))
  }
  
  
  count <- 0
  sd.prop <- 0.2
  max.count <- 500
  for (irep in -(nburn-1):ntot){
    #Step I: Sample draws from the posterior p(F| bullet) where F is BART
    y.hat <- y-mu.t
    
    weights.draw <- 1/sd.t^2
    
    sampler.list$setResponse(y.hat)
    sampler.list$setWeights(weights.draw)
    bart.draw <- sampler.list$run(0L, 1L) #single draw of conditional mean (BART)
    # Full fit
    f.draw <- bart.draw$train

    #Step III: Sample the infinite dim mixture model
    y.hat <- y - f.draw
    
    #Step 1: Sample mu.k and sd.K
    for (s in seq_len(G)){
      if (counts[s,1]>0){
        #Samples the component means
        v.1 <- 1/(counts[[s]] + zeta.0)
        mu.1 <- v.1*sum(y.hat[S==s])
        
        mu.draw <- mu.1 + rnorm(1, 0, sqrt(v.1)*sd.G[[s]])
        mu.G[[s]] <- mu.draw
        
        #Samples the component variances
        e.1 <- counts[[s]]/2 + g.0
        d.1 <- sum((y.hat[S==s]-mu.draw)^2)/2 + g.1
        sd.G[[s]] <- sqrt(1/rgamma(1, e.1, d.1))
      }else{
        mu.G[[s]] <- rnorm(1, mu.0, sqrt(1/zeta.0))
        sd.G[[s]] <- sqrt(1/rgamma(1, g.0, g.1))
      }
    }
    
 
    #Step 2:  Sample the sticks from beta distributions
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
    }
    #Step 3: Sample the classifications S
    u <- rep(0, n)
    fit <- rep(0,n)
    mu.t <- rep(0,n)
    sd.t <- rep(0,n)
    for (t in seq_len(n)){
      s.i <- S[[t]]
      xi.t <- xi[s.i]
      u.i <- runif(1, 0, xi.t)
      
      prob.non <- (u.i < xi[1:G])/xi[1:G] * eta[1:G] * dnorm(y.hat[[t]], mu.G[1:G], sd.G[1:G]) + 1e-32
      
      s.i <- sample(1:G,1, prob=prob.non/sum(prob.non), replace=TRUE)
      
      fit[[t]] <- mu.G[[s.i]] + rnorm(1, 0, sd.G[[s.i]])
      
      mu.t[[t]] <- mu.G[[s.i]]
      sd.t[[t]] <- sd.G[[s.i]]
      
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
    
    post.old <- G.plus*log(alpha)+(lgamma(alpha+1) - log(alpha))-lgamma(n+alpha)+sum(lgamma(counts[1:G.plus]))  + dgamma(alpha, a_alpha, b_alpha, log=TRUE)
    post.new <- G.plus*l_alpha_star+(lgamma(exp(l_alpha_star)+1) - l_alpha_star)-lgamma(n+exp(l_alpha_star))+sum(lgamma(counts[1:G.plus]))  + dgamma(exp(l_alpha_star), a_alpha, b_alpha, log=TRUE)
    
    acc <- post.new - post.old
    if (is.nan(acc)) acc <- 0
    
    if (acc > log(runif(1))){
      alpha <- exp(l_alpha_star)
    }
    b_k <- alpha
    
    #This part does SV estimation (if turned on) 
    if (SV.ident){
      svdraw <- svsample_general_cpp(y-f.draw, startpara = startpara, startlatent = sv.latent, priorspec = sv_priors)
  
      startpara[c("mu", "phi", "sigma")] <- as.list(svdraw$para[, c("mu", "phi", "sigma")])
    
      sv.latent <- svdraw$latent
      sd.t <- exp(as.numeric(sv.latent)/2)
      sd.G[[1]] <- sd.t[[n]]
    }
    
    shock <- mu.t + rnorm(n, 0, sd.t)
    
    #Store posterior draws
    if (irep > 0 && irep %% thinning == 0){
      alpha.store[irep/thinning,] <- alpha
      fit.store[irep/thinning,] <- f.draw
      yhat.store[irep/thinning,] <- f.draw + shock#mu.t + rnorm(1, 0, sd.t)#+ 
      shock.store[irep/thinning,] <- shock
      S.store[irep/thinning,] <- t(S)
      G.store[irep/thinning,] <- G.plus
      mu.store[irep/thinning,] <- mu.t
      sd.store[irep/thinning,] <- sd.t

      if (fcst){
        f.star.draw <- sampler.list$predict(as.numeric(Xho[1:(length(Xho)-1)])) #Put in stuff from BART here
        
        #Draw from the infinite mixture distribution based on the eta's
        s.star <- sample(1:G, 1, prob=eta[1:G])
        
        y.pred <- f.star.draw + rnorm(1, mu.G[[s.star]], sd.G[[s.star]])
        
        
        pred.store[irep/thinning, ] <- y.pred
        lpl.store[irep/thinning, ] <- dnorm(yho[[2]], f.star.draw + mu.G[[s.star]], sd.G[[s.star]], log=TRUE) #second element --> final vintage
      }
    }
    
    #print(irep)
    # print(count/(irep+nburn))
  }
  
  #Computes some posterior quantities
  LR <- ES <- matrix(NA, n, 1)
  for (i in 1:n){
    ES[i,] <- mean(yhat.store[yhat.store[,i] < quantile(yhat.store[,i],0.05),i])
    LR[i,] <- mean(yhat.store[yhat.store[,i] > quantile(yhat.store[,i],0.95),i])
  }
  
  S.median <- apply(S.store, 2, mean)#ts(apply(S.store, 2, mean), start=c(1973,2), frequency = 4)
  mu.mean <- apply(mu.store, 2, mean)#ts(apply(mu.store, 2, mean), start=c(1973,2), frequency = 4)
  sd.mean <- apply(sd.store, 2, mean)#ts(apply(sd.store, 2, mean), start=c(1973,2), frequency = 4)
  shock.mean <- t(apply((shock.store), 2, function(x) quantile(x, c(0.05, 0.5, 0.95))))#ts(t(apply((shock.store)*sd.y, 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))), start=c(1973,2), frequency = 4)
  
  if (fcst==FALSE){
    pred.store <- NULL
    fevd.store <- NULL
    lpl.store <- NULL
  } 
  return.list <- list("predictions"=pred.store, "fevds"=fevd.store, "lpl"=lpl.store, "ES_LR"=cbind(ES, LR), "S"=S.median, "mu.mix"=mu.mean, "sd.mix"=sd.mean, "shock"=shock.mean, "fit"=fit.store, "yho"=yho[[2]])
  
  return(return.list)
}


BART_mixSV <- function(y, X, yho, Xho,nburn=5000, ntot=10000, thinning=1, PDP=FALSE, sample.omega=TRUE, tau2.start=0.5, G.max=25, fcst=TRUE, a0=1, a1=1, PCA=TRUE, mix.sv=FALSE){
  #  y<- y.tilde; X<-X.new; yho <- NULL; Xho<- NULL; nburn=5000; ntot=10000; thinning=1; PDP=FALSE; sample.omega=TRUE; tau2.start=0.5; G.max=25; fcst=TRUE; a0=1; a1=1; PCA=TRUE
  #y=y; X=X; yho=yho; Xho=Xho;nburn=nburn; ntot=ntot; thinning=2; PDP=FALSE; sample.omega=TRUE; tau2.start=0.5; G.max=25; fcst=TRUE; a0=1; a1=1; PCA=PCA; mix.sv <- TRUE
  nsave <- floor(ntot / thinning)

  require(Matrix)
  require(dbarts)
  n <- dim(X)[1] #Sample size
  p <- dim(X)[2] #Number of covariates
  
  a <- 1/2
  b <- 1/2
  c.0 <- c.1 <- 0.6
  
  if (PCA){
    #The Optimal Hard Threshold for Singular Values is 4/sqrt(3)
    qstar <- min(6, ncol(X)) #in case of very few columns, do not use dim reduction
    
    svd.x <- svd(X)
    Xhat <- svd.x$u
    Xhat <- Xhat[, 1:qstar, drop=FALSE]
  
    K <- ncol(Xhat)
  }else{
    Xhat <- X
    K <- ncol(Xhat)
  }
  Q0 <- Xhat%*%solve(crossprod(Xhat)+diag(ncol(Xhat))*1e-4)%*%t(Xhat)
  Xhat <- X
  
  sigma2.draw <- rep(1e-3, n)
  Q0.ident <- (diag(n)-Q0)
  
  if (G.max =="SV"){
    require(stochvol)
    G.max <- 1
    SV.ident <- TRUE
    
    sv.draw <- list(mu = 0, phi = 0.6, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    sv.latent <- rep(0, n)
    
    sv_priors <- specify_priors(
      mu = sv_normal(mean = 0, sd = 10), # -4, 1e-3
      phi = sv_beta(shape1 = 25, shape2 = 5),
      sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*1)),
      nu = sv_infinity(),
      rho = sv_constant(0))
    
    startpara <- list(mu = 0, phi = 0.9, sigma = 0.1,
                      nu = Inf, rho = 0, beta = NA,
                      latent0 = 0)
  }else{
    SV.ident <- FALSE
    G.max <- G.max
  }
  
  #r.t <- matrix(0, n, 1)
  r.t <- rep(0, n)
  if (mix.sv){
    rt.draw <- list(mu = 0, phi = 0.6, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    
    
    sv_priors <- specify_priors(
      mu = sv_normal(mean = 0, sd = 10), # -4, 1e-3
      phi = sv_beta(shape1 = 25, shape2 = 5),
      sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*1)),
      nu = sv_infinity(),
      rho = sv_constant(0))
    
    startpara.rt <- list(mu = 0, phi = 0.9, sigma = 0.01,
                         nu = Inf, rho = 0, beta = NA,
                         latent0 = 0)
    
  }
  
  G <- G.max 
  kappa <- 0.8 #truncation parameter for the slice sampler 
  xi <- (1-kappa)*kappa^(seq(1, G)-1)
  #--------------------------Prior hyperparameters --------------------------
  alpha <- 0.5 # initialize at expected value
  a_k <- 1
  b_k <- alpha
  
  #Prior on alpha (standard choice stipulated in Escobar & West)
  a_alpha <- 2
  b_alpha <- 4
  
  #Prior on the the component means and variances
  mu.0 <- 0
  
  if (G==1) zeta.0 <- 10^5 else zeta.0 <- 0.01
  
  
  v <- 100
  
  mu.o <- 100
  phi <- 1
  o.1 <- 10#mu.o * phi
  o.2 <- 5#phi
  
  tau0 <- 2
  #-----------------------------Starting values for the infinite mixture -----------------------------
  if (G.max > 1){
    #Classification
    start.flexmix <- flexmix::flexmix((y-Q0%*%y) ~ 1, k=3)
    si <- start.flexmix@posterior$scaled
    S <- apply(si, 1, function(x) which(x==max(x)))#apply(t(t(si) * seq(1, ncol(si))), 1, sum)#rep(1, n)
    
    mu.G <- rep(0, G)
    sd.G <- rep(1, G)
    
    for (jj in seq_len(ncol(si))){
      mu.G[[jj]] <- parameters(start.flexmix)[1,jj]
      sd.G[[jj]] <- parameters(start.flexmix)[2,jj]
    }
    
    counts <- matrix(0, G, 1)
    S.tab <- table(S)
    counts[as.numeric(names(S.tab))] <- S.tab
    Nclust <- length(counts)	  # initial number of clusters
    
    mu.t <- rep(0,n)
    sd.t <- rep(0,n)
    
    for (t in seq_len(n)){
      s.i <- S[[t]]
      mu.t[[t]] <- mu.G[[s.i]]
      sd.t[[t]] <- sd.G[[s.i]]
    }
  }else{
    mu.t <- rep(0,n)
    sd.t <- rep(1e-4,n)
    mu.G <- rep(0, G)
    sd.G <- rep(1, G)
    
    counts <- matrix(0, 1, 1)
    counts[1,1] <- n
    S <- matrix(1, n, 1)
  }
  
  sd.raw.t <- sd.t
  
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
  
  #Conditional mean modeling through BART
  cgm.level <- 0.95
  cgm.exp   <- 2
  sd.mu     <- 2
  num.trees <- 250
  
  prior.sig <- c(10000^50, 0.5)
  sigma.init <- 1
  
  control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                           keepTrees = FALSE, n.samples = ntot,
                           n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
                           n.threads = 1, n.thin = 1L, printEvery = 1,
                           printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                           updateState = FALSE)
  
  sampler.list <- dbarts(y~X, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,n), sigma=sigma.init, resid.prior = chisq(prior.sig[[1]], prior.sig[[2]]))
  
  
  
  #Storage matrices
  S.store <- matrix(NA, nsave, n)
  mu.store <- matrix(NA, nsave, n)
  sd.store <- matrix(NA, nsave, n)
  yhat.store <- shock.store <- fit.store <- matrix(NA, nsave, n)
  alpha.store <- matrix(NA, nsave, 1)
  G.store <- matrix(NA, nsave, 1)
  omega.store <- matrix(NA, nsave, 1)
  theta.store <- matrix(NA, nsave, 2)

  if (fcst){
    pred.store <- array(NA, c(nsave, 1))
    fevd.store <- array(NA, c(nsave, 1))
    lpl.store <- array(NA, c(nsave, 1))
  }
  
  
  count <- 0
  sd.prop <- 0.2
  max.count <- 500
  for (irep in -(nburn-1):ntot){
    #Step I: Sample draws from the posterior p(F| bullet)
    y.hat <- y-mu.t
    
    weights.draw <- 1/sd.t^2
    
    sampler.list$setResponse(y.hat)
    sampler.list$setWeights(weights.draw)
    bart.draw <- sampler.list$run(0L, 1L) #single draw of conditional mean (BART)
  
    f.draw <-  bart.draw$train
    #mean((fhat - Q0%*%y)^2)
  
    #Step III: Sample the infinite dim mixture model
    y.hat <- y - f.draw
    
    #Step 1: Sample mu.k and sd.K
    for (s in seq_len(G)){
      if (counts[s,1]>0){
        #Samples the component means
        v.1 <- 1/(crossprod((S==s)*1/exp(r.t/2)) + zeta.0)
        mu.1 <- v.1*sum(y.hat[S==s]/exp(r.t[S==s]/2))
        
        mu.draw <- as.numeric(mu.1 + rnorm(1, 0, sqrt(v.1)*sd.G[[s]]))
        mu.G[[s]] <- mu.draw
        
        #Samples the component variances
        e.1 <- counts[[s]]/2 + o.1#v/2
        d.1 <- sum(((y.hat[S==s]-mu.draw)*1/exp(r.t[S==s]/2))^2)/2 + o.2#v/(2*tau0)
        
        if (!mix.sv) sd.i <- sqrt(1/rgamma(1, e.1, d.1)) else sd.i <- 1
        
        sd.G[[s]] <- sd.i#sqrt(1/rgamma(1, e.1, d.1))
      }else{
        
        if (!mix.sv) sd.i <- sqrt(1/rgamma(1, o.1,  o.2)) else sd.i <- 1
        
        mu.G[[s]] <- rnorm(1, mu.0, sqrt(1/zeta.0))
        sd.G[[s]] <- sd.i #1/rgamma(1, v/2,  v/(2*tau0))
      }
    }
    
    #Step 2:  Sample the sticks from beta distributions
    if (G.max > 1 && G > 1){
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
    }
    
    #Step 3: If we use a mixture SV model sample the latent components and associated parameters
    if (mix.sv){
      rt.draw <- svsample_general_cpp((y-f.draw-mu.t)/sd.raw.t, startpara = startpara.rt, startlatent = r.t, priorspec = sv_priors)
      
      startpara.rt[c("mu", "phi", "sigma")] <- as.list(rt.draw$para[, c("mu", "phi", "sigma")])
      #  print(rt.draw$para)
      
      r.t <- as.numeric(rt.draw$latent)
      r.t[r.t < -6] <- -6
      r.t[r.t > 1] <- 1
    }
    
    #Step 4: Sample the classifications S
    u <- rep(0, n)
    fit <- rep(0,n)
    mu.t <- rep(0,n)
    sd.t <- rep(0,n)
    sd.raw.t <- rep(0,n)
    for (t in seq_len(n)){
      s.i <- S[[t]]
      xi.t <- xi[s.i]
      u.i <- runif(1, 0, xi.t)
      
      prob.non <- (u.i < xi[1:G])/xi[1:G] * eta[1:G] * dnorm(y.hat[[t]]/exp(r.t[[t]]/2), mu.G[1:G], sd.G[1:G]) + 1e-32
      
      s.i <- sample(1:G,1, prob=prob.non/sum(prob.non), replace=TRUE)
      
      fit[[t]] <- mu.G[[s.i]] + rnorm(1, 0, sd.G[[s.i]])
      
      mu.t[[t]] <- mu.G[[s.i]]
      sd.t[[t]] <- exp(r.t[[t]]/2)*sd.G[[s.i]]
      sd.raw.t[[t]] <- sd.G[[s.i]]
      S[[t]] <- s.i
      u[[t]] <- u.i
    }
    
    #sd.t[sd.t<1e-4] <- 1e-4
    
    
    S.tab <- table(S)
    counts <- matrix(0,G.max,1)
    counts[as.numeric(names(S.tab))] <- S.tab
    G.plus <- sum(counts>0)
    #Does the truncation of G
    G <- sum((1-cumsum(eta)) < min(u))
    if (G.max == 1) G <- 1
    if (G==0) G <- 1
    
    
    #Step 5: Sample the precision parameter through an MH step
    l_alpha_star <- log(alpha) + rnorm(1,0,0.25)
    
    post.old <- G.plus*log(alpha)+(lgamma(alpha+1) - log(alpha))-lgamma(n+alpha)+sum(lgamma(counts[1:G.plus]))  + dgamma(alpha, a_alpha, b_alpha, log=TRUE)
    post.new <- G.plus*l_alpha_star+(lgamma(exp(l_alpha_star)+1) - l_alpha_star)-lgamma(n+exp(l_alpha_star))+sum(lgamma(counts[1:G.plus]))  + dgamma(exp(l_alpha_star), a_alpha, b_alpha, log=TRUE)
    
    acc <- post.new - post.old
    if (is.nan(acc)) acc <- 0
    
    if (acc > log(runif(1))){
      alpha <- exp(l_alpha_star)
    }
    b_k <- alpha
    
    #This part does SV estimation (if turned on) 
    if (SV.ident){
      svdraw <- svsample_general_cpp((y-f.draw), startpara = startpara, startlatent = sv.latent, priorspec = sv_priors)
      
      startpara[c("mu", "phi", "sigma")] <- as.list(svdraw$para[, c("mu", "phi", "sigma")])
      
      sv.latent <- svdraw$latent
      sd.t <- exp(as.numeric(sv.latent)/2)
      sd.G[[1]] <- sd.t[[n]]
    }
    
    shock <- mu.t + rnorm(n, 0, sd.t)
    
   
    
    #Store posterior draws
    if (irep > 0 && irep %% thinning == 0){
      alpha.store[irep/thinning,] <- alpha
      fit.store[irep/thinning,] <- f.draw
      yhat.store[irep/thinning,] <- f.draw + shock#mu.t + rnorm(1, 0, sd.t)#+ 
      shock.store[irep/thinning,] <- shock
      S.store[irep/thinning,] <- t(S)
      G.store[irep/thinning,] <- G.plus
      mu.store[irep/thinning,] <- mu.t
      sd.store[irep/thinning,] <- sd.t
      #Compute partial effects of the NFCI on GDP growth
      if (fcst){
        f.star.draw <- sampler.list$predict(as.numeric(Xho[1:(length(Xho)-1)])) #Put in stuff from BART here
        
        #Draw from the infinite mixture distribution based on the eta's
        s.star <- sample(1:G, 1, prob=eta[1:G])
        
        y.pred <- f.star.draw + rnorm(1, mu.G[[s.star]], sd.G[[s.star]])
        
        
        pred.store[irep/thinning, ] <- y.pred
        lpl.store[irep/thinning, ] <- dnorm(yho[[2]], f.star.draw + mu.G[[s.star]], sd.G[[s.star]], log=TRUE) #second element --> final vintage
      }
    }
    
    #print(irep)
    # print(count/(irep+nburn))
  }
  
  #Computes some posterior quantities
  LR <- ES <- matrix(NA, n, 1)
  for (i in 1:n){
    ES[i,] <- mean(yhat.store[yhat.store[,i] < quantile(yhat.store[,i],0.05),i])
    LR[i,] <- mean(yhat.store[yhat.store[,i] > quantile(yhat.store[,i],0.95),i])
  }
  
  S.median <- apply(S.store, 2, mean)#ts(apply(S.store, 2, mean), start=c(1973,2), frequency = 4)
  mu.mean <- apply(mu.store, 2, mean)#ts(apply(mu.store, 2, mean), start=c(1973,2), frequency = 4)
  sd.mean <- apply(sd.store, 2, mean)#ts(apply(sd.store, 2, mean), start=c(1973,2), frequency = 4)
  shock.mean <- t(apply((shock.store), 2, function(x) quantile(x, c(0.05, 0.5, 0.95))))#ts(t(apply((shock.store)*sd.y, 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))), start=c(1973,2), frequency = 4)
  
  if (fcst==FALSE){
    pred.store <- NULL
    fevd.store <- NULL
    lpl.store <- NULL
  } 
  return.list <- list("predictions"=pred.store, "fevds"=fevd.store, "lpl"=lpl.store, "ES_LR"=cbind(ES, LR), "S"=S.median, "mu.mix"=mu.mean, "sd.mix"=sd.mean, "shock"=shock.mean, "fit"=fit.store, "yho"=yho[[2]])
  
  return(return.list)
}