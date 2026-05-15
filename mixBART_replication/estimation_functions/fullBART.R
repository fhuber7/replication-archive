# ------------------------------------------------------------------------------------
# MCMC setup
fullBART <- function(Yraw, nburn=5000, nsave = 5000, thinfac = 1, prior.sig, sv = TRUE, fc.approx="exact", restr.var = NULL, fhorz = 3, quiet = TRUE){
  Z.sv <- NULL
  source("aux_func.R")
  # standardize data
  Ymu <- apply(Yraw, 2, mean,na.rm=T)
  Ysd <- apply(Yraw, 2, sd,na.rm=T)
  Yraw <- apply(Yraw, 2, function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})
  
  
  iter.update <- 250
  
  nthin <- round(thinfac * nsave)
  ntot <- nburn + nsave
  thin.set <- floor(seq(nburn+1,ntot,length.out=nthin))
  in.thin <- 0
  
  p <- 12# lag order
  cons <- 0 # whether constant should be included
  
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
  PHI_draw <- PHI_approx <- A_OLS <- solve(crossprod(X)+diag(ncol(X))*100)%*%crossprod(X,Y)
  A_prior <- matrix(0, K, M)
  Sig_OLS <- crossprod(Y-X%*%A_OLS)/T
  
  # covariance related objects
  eta <- matrix(NA,T,M)
  H <- matrix(-3,T,M)
  A0_draw <- A0_approx <- diag(M)#t(chol(Sig_OLS))
  ASV_approx <- matrix(0,K+1,M)
  
  sv_draw <- list()
  sv_latent <- list()
  for (mm in seq_len(M)){
    sv_draw[[mm]] <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    sv_latent[[mm]] <- rep(0,T)
  }
  
  sv_priors <- list()
  if(sv == "SV"){
    for(mm in 1:M){
      sv_priors[[mm]] <- specify_priors(
        mu = sv_normal(mean = 0, sd = 10),
        phi = sv_beta(shape1 = 5, shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*0.01)),
        nu = sv_infinity(),
        rho = sv_constant(0)
      ) 
    }
  }
  
  # priors
  #A_prior <- matrix(0,K,M)
  theta_A <- matrix(10,K,M)
  theta_A0 <- matrix(1,M,M)
  
  a0sig <- 3
  b0sig <- 0.3
  
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
  
  E.hat <- Y - X%*%PHI_draw
  sampler.list <- list()
  sampler.hetero.list <- list() #XX This part is new for heteroBART
  
  svdraw.list <- list()
  for (jj in seq_len(M)){
    if (sv == "heteroBART" || sv == "SV"){ 
      prior.sig = c(10000^50, 0.5)
      sigma.init <- 1
    }else{
      sigma.init <- sqrt(Sig_OLS[jj,jj]) 
    }
    
    if (jj > 1){
      X.j <- cbind(X, E.hat[,1:(jj-1)])
      sampler.list[[jj]] <- dbarts(Y[,jj]~X.j, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=sigma.init, resid.prior = chisq(prior.sig[[1]], prior.sig[[2]]))
    }else{
      X.j <- X
      sampler.list[[jj]] <- dbarts(Y[,jj]~ X.j , control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=sigma.init, resid.prior = chisq(prior.sig[[1]], prior.sig[[2]]))
    }
    
   # This part is new for heteroBART
   if (sv == "heteroBART") sampler.hetero.list[[jj]] <-  dbarts(E.hat[, jj] ~ cbind(seq(1, T), X), control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=3)
  }
  
  sampler.run <- list()
  sampler.sv.run <- list() #XX new for heteroBART
  if (sv == "heteroBART"){
    # Stuff necessary for the 10-component mixture)
    m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
    v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
    q_st  = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)
  }
  
  sigma.mat <- matrix(1, M, 1)
  count.mat <- matrix(0, M*p, M)

  # Initialize HS prior on covariances (if any)
  lambda.A0 <- 1
  nu.A0 <- 1
  tau.A0 <- 1
  zeta.A0 <- 1
  prior.cov <- rep(1, M*(M-1)/2)
  
  # initialize HS prior on linear VAR coefs
  lambda.VAR <- 1
  nu.VAR <- 1
  tau.VAR <- 1
  zeta.VAR <- 1
  sv_parms_mat <- matrix(NA, M, 4)
  Sig_t <- array(0,dim=c(T,M,M))
  
  Y.fit.BART <- Y.fit.VAR <- Y.fit.full <- Y*0
  
  # storage objects
  if (!is.null(restr.var)){
    sl.quants <- seq(0, 1, by=0.05)
    range.conditional <- quantile(Y[,restr.var], sl.quants)
    R <- length(range.conditional)
    H_store <- array(NA, dim=c(nthin, T, M)) #Store stochastic volatilities
    fcst_store <- Hfcst_store <- array(NA,dim=c(nthin,fhorz,M, R))
  }else{
    H_store <- array(NA, dim=c(nthin, T, M)) #Store stochastic volatilities
    fcst_store <- Hfcst_store <- array(NA,dim=c(nthin,fhorz,M))
  }
  
  # -----------------------------------------------------------------------------
  # start Gibbs sampler
  # show progress
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = ntot, style = 3)
    start  <- Sys.time()
  }
  
  F.ginv <- MASS::ginv(cbind(seq(1, T), X))
  for(irep in 1:ntot){
    # 1) sample model coefficients (either linear VAR or BART)
      for (mm in seq_len(M)){
        sampler.list[[mm]]$setResponse(Y[,mm])
        if (mm > 1) X.j <-cbind(X, eta[, 1:(mm-1)]) else X.j <- X
        sampler.list[[mm]]$setPredictor(X.j)
        X.ginv <- MASS::ginv(X.j)
        
        # This part estimates the eq-specific BART model
        rep_mm <- sampler.list[[mm]]$run(0L, 1L)
        sampler.run[[mm]] <- rep_mm
        
        sig.draw <- rep_mm$sigma
        sigma.mat[mm,] <- sig.draw
        
        draw_tmp <- as.numeric(X.ginv %*% rep_mm$train)
        PHI_approx[,mm] <- draw_tmp[1:K]
        if(mm>1){
          A0_approx[mm,1:(mm-1)] <- draw_tmp[(K+1):length(draw_tmp)]
        }
        
        if (any(is.na(rep_mm$train))) break
        Y.fit.full[,mm] <- rep_mm$train
        eta[,mm] <- Y[,mm] - rep_mm$train
        beta.draw <- X.ginv %*% rep_mm$train
        PHI_draw[,mm] <- beta.draw[1:K,]
        if (mm > 1) A0_draw[mm,1:(mm-1)] <- beta.draw[(K+1):nrow(beta.draw),]
      }
      shock_norm <- eta

      for (mm in seq_len(M)){
        if (sv == "SV"){
          if (mm > 1) X.j <-cbind(X, eta[, 1:(mm-1)]) else X.j <- X
          svdraw_mm <- svsample_general_cpp(shock_norm[,mm]/sigma.mat[mm], startpara = sv_draw[[mm]], startlatent = sv_latent[[mm]], priorspec = sv_priors[[mm]])
          sv_draw[[mm]][c("mu", "phi", "sigma")] <- as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
          sv_latent[[mm]] <- svdraw_mm$latent
          sv_parms_mat[mm, ] <- c(svdraw_mm$para[, c("mu", "phi", "sigma")], svdraw_mm$latent[T])
          
          normalizer <- as.numeric(exp(-.5*svdraw_mm$latent))
          weights.new <- as.numeric(exp(-svdraw_mm$latent))
          if (mm > 1){
            dat <- dbartsData(formula = Y[,mm] ~ X.j, weights=weights.new)
            sampler.list[[mm]]$setData(dat)
            H[,mm] <- log(sigma.mat[mm]^2) + svdraw_mm$latent
          }else{
            H[, mm] <- svdraw_mm$latent
          }
        }else if (sv =="heteroBART"){
          eta.star <- shock_norm[ , mm]/sigma.mat[mm]
          eta.off <- 1e-15 
          eta.star <- log((eta.star - eta.off)^2)
          
          # Sample the mixture components (previous draw of ht)
          z <- sapply(eta.star-H[,mm], ncind, m_st, sqrt(v_st2), q_st)
          
          # Subset mean and variances to the sampled mixture components; (n x p) matrices
          m_st_all = matrix(m_st[z],T, 1)
          v_st2_all = matrix(v_st2[z], T, 1)
          eta.star <- matrix((eta.star - m_st_all))#/sqrt(v_st2_all))
          
          #Do BART for conditional variances
          sampler.hetero.list[[mm]]$setResponse(eta.star)
          sampler.hetero.list[[mm]]$setWeights(1/v_st2_all)
          
          rep.hetero.i <- sampler.hetero.list[[mm]]$run(0L, 1L)
          sampler.sv.run[[mm]] <- rep.hetero.i
          
          ht.i <- rep.hetero.i$train#*sqrt(v_st2_all) #Off-setting if we get a really weird draw
          ht.i[ht.i < -10] <- -10
          ASV_approx[,mm] <- F.ginv%*%ht.i
          
          weights.new <- as.numeric(exp(-ht.i))
          sampler.list[[mm]]$setWeights(weights.new)
          H[ , mm] <- ht.i + log(as.numeric(sigma.mat[mm])^2)
        }else{
          H[,mm] <- log(sigma.mat[mm]^2)
        }
      }
      
      for(tt in seq_len(T)){
        S_tmp <- exp(H[tt,])
        S.t <- A0_draw%*%tcrossprod(diag(S_tmp),(A0_draw))
        Sig_t[tt,,] <- S.t
      }
      
    if(irep %in% thin.set){
      in.thin <- in.thin+1
      H_store[in.thin,,] <- H
      
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
          HT <-  H[T, ] - log(as.numeric(sigma.mat)^2)#the final column of sv_parms_mat stores the final value of the log volatilitieslog(sigma.mat^2) +
        }else{
          HT <- log(as.numeric(sigma.mat)^2)
        }
        
        Sig_T <- Sig_t[T,,] # use final observation for Sigma
        tree.pred <- matrix(0, M)
        shocks.pred <- matrix(0, M, 1)
        for (hh in seq_len(fhorz)){
          if (sv == "SV"){
            HT  <- log(as.numeric(sigma.mat)^2) + (sv_parms_mat[, 1] + sv_parms_mat[ , 2] * (HT - sv_parms_mat[,1]) + sv_parms_mat[ , 3]*rnorm(M))   #("mu", "phi", "sigma")
          }else if (sv == "heteroBART"){
            # if(fc.approx=="exact"){
              vola.predict.tree <- matrix(0, M, 1)
              for (nn in seq_len(M)){
                vola.predict.tree[nn,] <- sampler.hetero.list[[nn]]$predict(c(T+hh, X.hat))#
              }
              HT <- log(as.numeric(sigma.mat)^2) + vola.predict.tree
            # }else if(fc.approx=="approx"){
            #   HT <- log(as.numeric(sigma.mat)^2) + as.numeric(c(T+hh, X.hat)%*%ASV_approx)
            # }
          }
          Hfc[hh,] <- exp(HT)
          
          Y.tp1 <- matrix(0, M, 1)
          for (j in seq_len(M)){
            if(fc.approx=="exact"){
              if (j > 1){
                shocks.pred[j, 1] <-  rnorm(1, 0, exp(.5 * HT[j]))
                Y.tp1[j, 1] <- sampler.list[[j]]$predict(c(X.hat, shocks.pred[1:(j-1)])) + shocks.pred[j, 1]
              }else{
                shocks.pred[j, 1] <-  rnorm(1, 0, exp(.5 * HT[1]))
                Y.tp1[j, 1] <- sampler.list[[j]]$predict(X.hat) + shocks.pred[j, 1] 
              }
            }else if(fc.approx=="approx"){
              if(j>1){
                shocks.pred[j, 1] <-  rnorm(1, 0, exp(.5 * HT[j]))
                Y.tp1[j, 1] <- sum(c(X.hat,shocks.pred[1:(j-1)]) * c(PHI_approx[,j],A0_approx[j,1:(j-1)])) + shocks.pred[j, 1]
              }else{
                shocks.pred[j, 1] <-  rnorm(1, 0, exp(.5 * HT[1]))
                Y.tp1[j, 1] <- sum(X.hat*PHI_approx[,j]) + shocks.pred[j, 1] 
              }
            }
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
            HT <-  H[T, ] - log(as.numeric(sigma.mat)^2)#the final column of sv_parms_mat stores the final value of the log volatilitieslog(sigma.mat^2) +
          }else{
            HT <- log(as.numeric(sigma.mat)^2)
          }
          
          Sig_T <- Sig_t[T,,] # use final observation for Sigma
          tree.pred <- matrix(0, M)
          shocks.pred <- matrix(0, M, 1)
          for (hh in seq_len(fhorz)){
            if (sv == "SV"){
              HT  <- log(as.numeric(sigma.mat)^2) + (sv_parms_mat[, 1] + sv_parms_mat[ , 2] * (HT - sv_parms_mat[,1]) + sv_parms_mat[ , 3]*rnorm(M))   #("mu", "phi", "sigma")
            }else if (sv == "heteroBART"){
              vola.predict.tree <- matrix(0, M, 1)
              for (nn in seq_len(M)){
                vola.predict.tree[nn,] <- sampler.hetero.list[[nn]]$predict(c(T+hh, X.hat))#
              }
              HT <- log(as.numeric(sigma.mat)^2) + vola.predict.tree
            }
            Hfc[hh,,r] <- exp(HT)*Ysd
            
            Y.tp1 <- matrix(0, M, 1)
            for (j in seq_len(M)){
              if (j > 1){
                shocks.pred[j, 1] <-  rnorm(1, 0, exp(.5 * HT[j]))
                Y.tp1[j, 1] <- sampler.list[[j]]$predict(c(X.hat, shocks.pred[1:(j-1)])) + shocks.pred[j, 1]
              }else{
                shocks.pred[j, 1] <-  rnorm(1, 0, exp(.5 * HT[1]))
                Y.tp1[j, 1] <- sampler.list[[j]]$predict(X.hat) + shocks.pred[j, 1] 
              }
              Y.tp1[sl.cond] <- cond.vals
            } 
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
  dimnames(fcst_store) <- dimnames(Hfcst_store) <- list(paste0("mcmc",1:nthin), paste0("fhorz", 1:fhorz), colnames(Y))
  return_obj <- list("fcst" = fcst_store, "Hfcst" = Hfcst_store)
  return(return_obj)
}
