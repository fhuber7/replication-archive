ucsv <- function(y,  yho, fhorz, G.max ="SV", mix.SV=FALSE){
# y <- rnorm(150)
# yho <- c(1,2)
# fhorz <- 4
# G.max= "SV" 
# mix.SV = FALSE

  ck1994 <- function(y,Z,Ht,Q,m,p,t,B0,V0){
    # Carter and Kohn (1994), On Gibbs sampling for state space models.
    # Kalman Filter
    bp <- B0 #the prediction at time t=0 is the initial state
    Vp <- V0 #Same for the variance
    bt <- matrix(0,t,m) #Create vector that stores bt conditional on information up to time t
    Vt <- matrix(0,m^2,t) #Same for variances
    for (i in 1:t){
      R <- Ht[i] #CHK LATER NOISE INNOVATION VARIANCE TO MEASUREMENT ERRORS
      Qt <- Q[i]
      H <- Z[i]
      
      #Prediction steps
      cfe <- y[i] - H%*%bp   # conditional forecast error
      f <- H%*%Vp%*%t(H) + R    # variance of the conditional forecast error
      
      #Updating steps
      inv_f <- t(H)%*%solve(f)
      
      btt <- bp + Vp%*%inv_f%*%cfe  #updated mean estimate for btt Vp * inv_F is the Kalman gain
      Vtt <- Vp - Vp%*%inv_f%*%H%*%Vp #updated variance estimate for btt
      if (i < t){
        bp <- btt # prediction of tau(t+1)|info up to time t = tau(t)|info up to time t
        Vp <- Vtt + Qt
      }
      bt[i,] <- t(btt)
      Vt[,i] <- matrix(Vtt,m^2,1)
    }
    
    # draw the final value of the latent states using the moments obtained from the KF filters' terminal state
    bdraw <- matrix(0,t,m)
    bdraw[t,] <- btt+t(chol(Vtt))%*%rnorm(nrow(Vtt))
    
    #Now do the backward recurssions
    for (i in 1:(t-1)){
      Qt <- Q[t-i]
      bf <- t(bdraw[t-i+1,])
      btt <- t(bt[t-i,])
      Vtt <- matrix(Vt[,t-i,drop=FALSE],m,m)
      f <- Vtt + Qt
      
      inv_f <- Vtt%*%solve(f)
      cfe <- bf - btt
      
      bmean <- t(btt) + inv_f%*%t(cfe)
      bvar <- Vtt - inv_f%*%Vtt
      
      bdraw[t-i,] <- bmean+t(chol(bvar))%*%rnorm(nrow(bvar))
    }
    
    return(bdraw)
  }
  #Loads some packages needed
  require(stochvol) #to sample the log-volatilities and the coefficients of the state equation
  require(bvarsv)
  
  #The model we are going to estimate is the unobserved components stochastic volatility model of stock and watson 2007
  #                   pi(t) = tau(t) + epsilon(t)
  #                   tau(t)=tau(t-1)+v(t)
  #                   epsilon(t) ~ N(0, e^h(t))
  #                   v(t) ~ N(0, e^s(t))
  
  #Some MCMC preliminaries used for estimation
  nsave <- 5000
  nburn <- 5000
  ntot <- nsave+nburn
  
  #Specifies some priors
  b0 <- 0 #initial state equals zero
  Vb <- 10 #variance on initial state equals ten, quite uninformative
  
  
  #Get some key dimensions of the data
  n <- length(y) #we need to substract 1 below because we lose one observation when we calculate f(t) -f(t-1)
  
  #Before we start the Gibbs loop, we need to create storage matrices to store the MCMC output
  f_store <- matrix(NA,nsave,n) #stores the latent factors
  h_store <- matrix(NA,nsave,n) #stores the log-volatilities of the measurement errors
  s_store <- matrix(NA,nsave,n) #stores the log-volatilities of the errors to the states
  y_store <- matrix(NA, nsave, n) #stores the predicted value of inflation
  pred.store <- array(NA, c(nsave, 1))
  lpl.store <- array(NA, c(nsave, 1))
  
  #Initialize  inflation equation SV part
  sv.pi.draw <- list(mu = 0, phi = 0.9, sigma = 0.1, nu = Inf, rho = 0, beta = NA, latent0 = 0)
  sv.pi.latent <- rep(0, n)
  
  #Sets priors for measurement error variance
  if (G.max == "SV") B <- 0.01 else B <- 1e-10
  prior.pi <- specify_priors(
    mu = sv_normal(mean = 0, sd = 10),
    phi = sv_beta(shape1 = 25, shape2 = 1.5),
    sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*B)), #G(1/2, 1/(2*B)) \equiv  +/- sigma ~ N(0, B)
    nu = sv_infinity(),
    rho = sv_constant(0),
    latent0_variance = "stationary",
    beta = sv_multinormal(mean = 0, sd = 10000, dim = 1)
  )
  
  
  #Sets priors for state innovation error variances
  sv.tau.draw <- list(mu = 0, phi = 0.9, sigma = 0.1, nu = Inf, rho = 0, beta = NA, latent0 = 0)
  sv.tau.latent <- rep(0, n)
  C <- 0.01
  prior.tau <- specify_priors(
    mu = sv_normal(mean = 0, sd = 10),
    phi = sv_beta(shape1 = 25, shape2 =  1.5),
    sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*C)),
    nu = sv_infinity(),
    rho = sv_constant(0),
    latent0_variance = "stationary",
    beta = sv_multinormal(mean = 0, sd = 10000, dim = 1)
  )
  
  mu.t <- rep(0, n)
  
  if (G.max != "SV" || !is.na(as.numeric(G.max))){
  G <- G.max
  #G <- G.max 
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
  
  
  #v <- 100
  
  #mu.o <- 100
  #phi <- 1
  o.1 <- 10#mu.o * phi
  o.2 <- 5#phi
  
 # tau0 <- 2
  
  
  #This code is necessary if we estimate a mixture model for the shocks of the obs eq
  #Classification
  start.flexmix <- flexmix::flexmix((y) ~ 1, k=3)
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
  counts <- matrix(counts)
  Nclust <- length(counts)	  # initial number of clusters
  
  mu.t <- rep(0,n)
  sd.t <- rep(0,n)
  
  for (t in seq_len(n)){
    s.i <- S[[t]]
    mu.t[[t]] <- mu.G[[s.i]]
    sd.t[[t]] <- sd.G[[s.i]]
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
  
  vola.draw <- list(mu = 0, phi = 0.6, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
  
  
  sv_priors <- specify_priors(
    mu = sv_normal(mean = 0, sd = 10), # -4, 1e-3
    phi = sv_beta(shape1 = 25, shape2 = 5),
    sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*1)),
    nu = sv_infinity(),
    rho = sv_constant(0))
  
  startpara.vola <- list(mu = 0, phi = 0.9, sigma = 0.01,
                       nu = Inf, rho = 0, beta = NA,
                       latent0 = 0)
  
  r.t <- rep(0, n)
  sv.pi.latent <- log(sd.t^2)
  
  }
  
  #Now we start the big MCMC loop
  for (irep in 1:ntot){
    #Step I: Draw the latent non-stationary component using the CK/FS 1994 algorithm
    f.draw <- ck1994(y = (y-mu.t), Z = matrix(1,n,1), Ht = exp(sv.pi.latent), Q = exp(sv.tau.latent),m=1,p = 1,t=n,B0 = b0,V0 = Vb)
    
    #Step II: Draw log-volatilities and coefficients associated with the state EQ of the volas in the measurement equation
    if (G.max == "SV" || G.max == 1){
      u.draw <- y- f.draw - mu.t
      h_parm <- svsample_fast_cpp(u.draw,startpara=sv.pi.draw, startlatent=sv.pi.latent, priorspec = prior.pi)
    
      sv.pi.draw[c("mu", "phi", "sigma")] <- as.list(h_parm$para[, c("mu", "phi", "sigma")])
      sv.pi.latent <- h_parm$latent
      sv.pi.latent[exp(sv.pi.latent) < 1e-5] <- log(1e-5)
    }else{
      y.hat <-  y- f.draw - mu.t
      
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
          
          if (!mix.SV) sd.i <- sqrt(1/rgamma(1, e.1, d.1)) else sd.i <- 1
          
          sd.G[[s]] <- sd.i#sqrt(1/rgamma(1, e.1, d.1))
        }else{
          
          if (!mix.SV) sd.i <- sqrt(1/rgamma(1, o.1,  o.2)) else sd.i <- 1
          
          mu.G[[s]] <- rnorm(1, mu.0, sqrt(1/zeta.0))
          sd.G[[s]] <- sd.i #1/rgamma(1, v/2,  v/(2*tau0))
        }
      }
      
      if (mix.SV){
          vola.draw <- svsample_general_cpp(y.hat, startpara = startpara.vola, startlatent = r.t, priorspec = sv_priors)
          
          startpara.vola[c("mu", "phi", "sigma")] <- as.list(vola.draw$para[, c("mu", "phi", "sigma")])
          #  print(rt.draw$para)
          
          r.t <- as.numeric(vola.draw$latent)
          r.t[r.t < -6] <- -6
          r.t[r.t > 1] <- 1
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
          }else if (j==2){
            eta[[k]] <- (1-nu[[k-1]])*nu[[k]]
          }else if (j>2){
            eta[[k]] <-  nu[[k]]*prod(1-nu[1:(k-1)])
          }
        }
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
      sv.pi.latent <- log(sd.t^2)
      
      S.tab <- table(S)
      counts <- matrix(0,G.max,1)
      counts[as.numeric(names(S.tab))] <- S.tab
      G.plus <- sum(counts>0)
      #Does the truncation of G
      G <- sum((1-cumsum(eta)) < min(u))
      if (G.max == 1) G <- 1
      if (G == 0) G <- 1
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
      
    }
    
    #Step III: Draw log-volatilities and coeffs in the state equation
    e.draw <- f.draw[2:n,1]-f.draw[1:(n-1),1]
    e.draw <- c(e.draw[1],e.draw) #DIRTY HACK, assumes that in the initial state the error equals the error in the second point in time
    
    s_parm <-   svsample_fast_cpp(e.draw,startpara=sv.tau.draw, startlatent=sv.tau.latent, priorspec = prior.tau)
    sv.tau.draw[c("mu", "phi", "sigma")] <- as.list(s_parm$para[, c("mu", "phi", "sigma")])
    sv.tau.latent <- s_parm$latent
    sv.tau.latent[exp(sv.tau.latent) < 1e-5] <- log(1e-5)
    
    if (irep>nburn){ #If irep exceeds the amounts of burn-ins needed, we start to store our draws
      f_store[irep-nburn,] <- f.draw
      h_store[irep-nburn,] <- t(sv.pi.latent)
      s_store[irep-nburn,] <- t(sv.tau.latent)
      y_store[irep-nburn,] <- f.draw + rnorm(n, 0, exp(.5*t(sv.pi.latent[n])))
      
      #compute h-step-ahead predictions for the volatilities
      parms.tau <- sv.tau.draw[c("mu","phi","sigma")]
      if (G.max == "SV"){ 
        parms.sig <- sv.pi.draw[c("mu","phi","sigma")] 
        sv.sig <- sv.pi.latent[n]
      }else{
        parms.sig <- startpara.vola[c("mu","phi","sigma")] 
        sv.sig <- r.t[n]
      } 
      sv.h <- sv.tau.latent[n]
      sv.mat.h <- matrix(NA, fhorz,1)
      sv.mat.sig <- matrix(NA, fhorz, 1)
      for (h in 1:fhorz){
        sv.h <- parms.tau$mu + parms.tau$phi*(sv.h - parms.tau$mu) + rnorm(1,0, parms.tau$sigma)
        if (G.max == "SV" || mix.SV){
          if (G.max == "SV"){
            sv.sig <- parms.sig$mu + parms.sig$phi*(sv.sig - parms.sig$mu) + rnorm(1,0, parms.sig$sigma)
          }else{
            sv.sig <- parms.sig$mu + parms.sig$phi*(sv.sig - parms.sig$mu) + rnorm(1,0, parms.sig$sigma)
          } 
        } else {
          s.star <- sample(1:G, 1, prob=eta[1:G])
          sv.sig <- log(sd.G[[s.star]]^2)
        }
        
        sv.mat.h[h,] <- sv.h
        sv.mat.sig[h,] <- sv.sig
      }
      
      pred.store[irep-nburn,] <- f.draw[[n]]+rnorm(1, 0, sqrt(sum(exp(sv.mat.h)))) + rnorm(1, 0, sqrt(sum(exp(sv.mat.sig))))
      lpl.store[irep-nburn,] <- dnorm(yho[[2]], f.draw[[n]], sqrt(sum(exp(sv.mat.h))+sum(exp(sv.mat.sig))), log=TRUE) 
    }
    # print(irep) 
  }
  
  return.list <- list("predictions"=pred.store,  "lpl"=lpl.store, "yho"=yho[[2]])
  
  return(return.list)
}
