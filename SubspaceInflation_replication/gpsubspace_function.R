gp_bnp <- function(y, X, yho, Xho,nburn=5000, ntot=10000, thinning=1, PDP=FALSE, sample.omega=TRUE, tau2.start=0.5, G.max=25, fcst=TRUE, a0=1, a1=1, PCA=TRUE, mix.sv =FALSE){
#  y<- y.tilde; X<-X.new; yho <- NULL; Xho<- NULL; nburn=5000; ntot=10000; thinning=1; PDP=FALSE; sample.omega=TRUE; tau2.start=0.5; G.max=25; fcst=TRUE; a0=1; a1=1; PCA=TRUE
#  y=y; X=X; yho=yho; Xho=Xho;nburn=nburn; ntot=ntot; thinning=2; PDP=FALSE; sample.omega=TRUE; tau2.start=0.5; G.max=G.slct; fcst=TRUE; a0=1; a1=1; PCA=PCA
   nsave <- floor(ntot / thinning)

  #G.max <- "SV"
  #---- Stuff for partial dependence plots
  if (PDP){
  sl.p <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9,0.95, 0.99)#seq(0.05, 0.95, length.out=10)
  nfci.grid <- quantile(X[,3], sl.p)
  y.grid <- quantile(X[,2], sl.p)#seq(min(X[,2]), max(X[,2]), length.out=10)
  
  N.0  <- length(nfci.grid)
  N.1 <- length(y.grid)
  }
  require(Matrix)
  n <- dim(X)[1] #Sample size
  p <- dim(X)[2] #Number of covariates
  
  a <- 1/2
  b <- 1/2
  
  med.mat <- matrix(NA, n, n)
  for (j in seq_len(n)){
    for (i in seq_len(n)){
      med.mat[j,i] <- sqrt(sum((X[j,] - X[i,])^2))
    }
  }
  scale <- 1/median(med.mat[lower.tri(med.mat)])
  scale.linear <- 0.999
  
  Kn <- GaussKernel_star(t(X), t(X), scale, scale.linear)
  diag(Kn) <- 1  
  Kn.inv <- solve(Kn)
  
  Khat <- t(chol(Kn, pivot=TRUE))
  
  tau2 <- tau2.start
  if (PCA){
    #The Optimal Hard Threshold for Singular Values is 4/sqrt(3)
   qstar <- min(6, ncol(X)) #in case of very few columns, do not use dim reduction
   
   svd.x <- svd(X)
   Xhat <- svd.x$u
   Xhat <- Xhat[, 1:qstar, drop=FALSE]
   
   Xhat.ho <- svd(rbind(Xho, X))$u
   Xhat.ho <- Xhat.ho[, 1:qstar, drop=FALSE] 
   K <- ncol(Xhat)
  }else{
    Xhat <- X
    Xhat.ho <- rbind(X, Xho)
    K <- ncol(Xhat)
  }
  Q0 <- Xhat%*%solve(crossprod(Xhat)+diag(ncol(Xhat))*1e-4)%*%t(Xhat)
  X.big <- Xhat.ho%*%solve(crossprod(Xhat.ho)+diag(ncol(Xhat))*1e-4)%*%t(Xhat.ho)
  
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
    sd.t <- rep(0,n)
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

  
  #Storage matrices
  S.store <- matrix(NA, nsave, n)
  mu.store <- matrix(NA, nsave, n)
  sd.store <- matrix(NA, nsave, n)
  yhat.store <- shock.store <- fit.store <- matrix(NA, nsave, n)
  alpha.store <- matrix(NA, nsave, 1)
  G.store <- matrix(NA, nsave, 1)
  omega.store <- matrix(NA, nsave, 1)
  theta.store <- matrix(NA, nsave, 2)
  if (PDP) pd.store <- array(NA, c(nsave, N.0, N.1))
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
    K1 <- try(solve(Kn.inv + Q0.ident/tau2), silent=TRUE)
    if (is(K1, "try-error")) K1 <- MASS::ginv(Kn.inv + Q0.ident/tau2)
    Dinv <- solve(K1 + diag(sd.t^2)) #
    fhat <- K1 %*% Dinv %*% (y-mu.t)#solve(Kn + diag(sigma2,n), y)
    Vhat <- K1 - K1 %*% Dinv %*% t(K1)#solve(Kn+diag(sigma2,n),Kn)
    f.draw <-  fhat + t(chol(Vhat, pivot=TRUE))%*%rnorm(n)
    mean((fhat - Q0%*%y)^2)
    
    #Step II: Sample omega using inverse transform sampling (notice that prior does not depend on sigma2)
    if (sample.omega){
      check <- TRUE
      omega.count <- 0
      while (check){ #  && omega.count < max.count
        omega.count <- omega.count+1
        eta.j <- 1/tau2
        t.j <- (eta.j+1)^(-a-b)
        
        u.j <- runif(1, 0, t.j)
        t.j.star <- u.j^(1/(-(a+b)))-1
        
        a.0 <- a+(n-K)/2
        b.0 <- as.numeric(t(f.draw)%*%(Q0.ident)%*%f.draw)/2
        
        eta.j.star <- rgamma(1, a.0, b.0)
        if (eta.j.star < t.j.star) check <- FALSE 
      }
      tau2 <- 1/sqrt(eta.j.star)
    }else{
      tau2 <- tau2.start
    }
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
    if (G == 0) G <- 1
    
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
    
    #Step 5: Sample the hyperparameters of the Kernel matrix
    y.hat <- y-mu.t
    ml.old <-  - 1/2 * t(y.hat) %*% Dinv %*% y.hat - 1/2 * log(1/det(Dinv)) - n/2 * log(2*pi)+ dbeta(scale, a0,a1,log=TRUE) + dbeta(scale.linear, a0, a1, log=TRUE)
    
    #Construct proposal distribution using a two-dimensional random walk MH update
    theta.star <- rnorm(2, c(scale, scale.linear), sd.prop)
    scale.prop <- theta.star[[1]]; scale.linear.prop <- theta.star[[2]]
    Kn.prop <- GaussKernel_star(t(X), t(X), scale.prop, scale.linear.prop)
    diag(Kn.prop) <- 1 
    
    Kn.inv.prop <- try(solve(Kn.prop), silent=TRUE)
    if (is(Kn.inv.prop, "try-error")) ml.new <- -Inf else{
      K1.prop <- solve(Kn.inv.prop + Q0.ident/tau2)
      Dinv.prop <- solve(K1.prop + diag(sd.t^2))
      ml.new <-   -1/2 * t(y.hat) %*% Dinv.prop %*% y.hat - 1/2 * log(1/det(Dinv.prop)) - n/2 * log(2*pi) + dbeta(scale.prop, a0,a1,log=TRUE) + dbeta(scale.linear.prop, a0, a1, log=TRUE)
      if (is.na(ml.new)) ml.new <- -Inf
    }
    
    log.post.ratio <- (ml.new - ml.old)
    
    if (is.na(log.post.ratio)) log.post.ratio <- -Inf
    
    if (log.post.ratio > log(runif(1))){
      scale <- scale.prop
      scale.linear <- scale.linear.prop
      Kn <- Kn.prop
      Kn.inv <- Kn.inv.prop
      K1 <- K1.prop
      count <- count+1
      #Dinv <- solve(K1 + diag(sd.t^2))
    }
    
    if (irep < 0){
      if (count/(irep+nburn) < 0.2) sd.prop <- 0.99 * sd.prop
      if (count/(irep+nburn) > 0.4) sd.prop <- 1.01 * sd.prop
    }
    
    #Store posterior draws
    if (irep > 0 && irep %% thinning == 0){
      alpha.store[irep/thinning,] <- alpha
      fit.store[irep/thinning,] <- f.draw
      yhat.store[irep/thinning,] <- f.draw + shock#mu.t + rnorm(1, 0, sd.t)#+ 
      shock.store[irep/thinning,] <- shock
      S.store[irep/thinning,] <- t(S)
      G.store[irep/thinning,] <- G.plus
      omega.store[irep/thinning,] <- 1/(1+tau2)
      mu.store[irep/thinning,] <- mu.t
      sd.store[irep/thinning,] <- sd.t
      theta.store[irep/thinning,] <- c(scale, scale.linear)
      #Compute partial effects of the NFCI on GDP growth
      if (fcst){
        K.star <- GaussKernel_star(t(rbind(Xho,X)), t(rbind(Xho,X)), scale, scale.linear)
        diag(K.star) <- 1
        
        K.big <- solve(solve(K.star) + (diag(n+1)-X.big)/tau2)
        Kn.star <- K.big[1, 2:(n+1), drop=F]
        KK.star <- K.big[1, 1]
        
        f.star.mean <- Kn.star %*% Dinv %*% (y - mu.t)
        V.f.mean <- KK.star - Kn.star %*% Dinv %*% t(Kn.star)
        
        #if (V.f.mean < 0) V.f.mean <- 1e-6
        
        f.star.draw <- f.star.mean + rnorm(1, 0, sqrt(V.f.mean))
        
        #Draw from the infinite mixture distribution based on the eta's
        s.star <- sample(1:G, 1, prob=eta[1:G])
        
        y.pred <- f.star.draw + rnorm(1, mu.G[[s.star]], sd.G[[s.star]])
        
        
        pred.store[irep/thinning, ] <- y.pred
        fevd.store[irep/thinning, ] <- sqrt(V.f.mean)/(sqrt(V.f.mean) + sd.G[[s.star]])
        lpl.store[irep/thinning, ] <- dnorm(yho[[2]], f.star.draw + mu.G[[s.star]], sd.G[[s.star]], log=TRUE) #second element --> final vintage
      }
     # browser()
      
      if (PDP){
        pd.mat <- matrix(NA, N.0, N.1)
        for (ff in seq_len(N.1)){  
          X.star <- cbind(1,y.grid[[ff]], nfci.grid)
          X.big <- rbind(X.star, X)%*%solve(crossprod(rbind(X.star, X)))%*%t(rbind(X.star, X))
          # Kn.star <- GaussKernel_star(t(X.star), t(X), 0.25, 0.5)
          # KK.star <- GaussKernel_star(t(X.star), t(X.star), 0.25, 0.5)
          # 
          K.star <- GaussKernel_star(t(rbind(X.star,X)), t(rbind(X.star,X)), scale, scale.linear)
          #diag(K.star) <- 1
          
          K.big <- solve(solve(K.star) + (diag(n+N.0)-X.big)/tau2)
          Kn.star <- K.big[1:N.0, (N.0+1):(n+N.0)]
          KK.star <- K.big[1:N.0, 1:N.0]
          
          f.star.mean <- Kn.star %*% Dinv %*% (y - mu.t)
          V.f.mean <- KK.star - Kn.star %*% Dinv %*% t(Kn.star)
          f.star.draw <- f.star.mean + t(chol(V.f.mean, pivot=TRUE)) %*%rnorm(N.0)
          
          pd.mat[,ff] <- f.star.draw
        }
        pd.store[irep/thinning,,] <- pd.mat
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