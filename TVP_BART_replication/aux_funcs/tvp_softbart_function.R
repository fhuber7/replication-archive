#### --------------------------------------------------------------------- ####
#### -------------------- TVP-BART VAR: SoftBart sampler ----------------- ####
#### --------------------------------------------------------------------- ####
tvp.softbart.func <- function(y = y, X = X, Z = Z, Z.fsv = Z.fsv, Z.idio = Z.idio,
                              p = p, cons = cons, thin = thin,
                              nsave = nsave, nburn = nburn,
                              model.setup = model.setup, bart.setup = bart.setup,
                              off.setup = off.setup, irf.setup = irf.setup){

ids <- colnames(y)

list2env(model.setup, globalenv())
list2env(bart.setup,  globalenv())
list2env(off.setup,   globalenv())
list2env(irf.setup,   globalenv())

ntot     <- nsave*nthin + nburn
save.set <- seq(nthin, nsave*nthin, nthin) + nburn
save.len <- length(save.set)
save.ind <- 0

#### --------------------------------------------------------------------- ####
#### ------------------- Dimensions and data matrices -------------------- ####
#### --------------------------------------------------------------------- ####

K    <- ncol(X)
M    <- ncol(y)
M.pr <- length(prices)
N    <- nrow(y)
# Full coefficient dimension
KK <- K + K*Q + R

# Storage for the time-varying full regressor block
Xfull.i <- matrix(NA, N, KK)

# Kronecker expansion of X for Q-factor TVP mean: stacks Q copies of X column-wise
XX <- kronecker(t(rep(1, Q)), X)


#### --------------------------------------------------------------------- ####
#### --------------------------- Prior setup ----------------------------- ####
#### --------------------------------------------------------------------- ####

V0inv <- matrix(1, M, KK)

# 10-component Normal mixture approximation to log-chi-squared
m_st  <- c( 1.92677,  1.34744,  0.73504,  0.02266, -0.85173,
           -1.97278, -3.46788, -5.55246, -8.68384,-14.65000)
v_st2 <- c(0.11265, 0.17788, 0.26768, 0.40611, 0.62699,
           0.98583, 1.57469, 2.54498, 4.16591, 7.33342)
q_st  <- c(0.00609, 0.04775, 0.13057, 0.20674, 0.22715,
           0.18842, 0.12047, 0.05591, 0.01575, 0.00115)

#### --------------------------------------------------------------------- ####
#### -------------- BART forest initialisation: mean model --------------- ####
#### --------------------------------------------------------------------- ####

# M equations, each with Q trees; penalty on tree complexity increases with factor index
forest.list <- list()
for(ii in 1:M) forest.list[[ii]] <- list()
names(forest.list) <- colnames(y)

for(ii in seq_len(M)){
  for(jj in seq_len(Q)){
    if(Q.penalty){
      cgm.exp0   <- cgm.exp^jj   # tree depth penalty, increasingly stringent per factor
      cgm.level0 <- cgm.level^jj # tree complexity penalty
    }else{
      cgm.exp0   <- cgm.exp
      cgm.level0 <- cgm.level
    }

    hypers <- Hypers(X        = Z,
                     Y        = y[,ii],
                     alpha    = 1,
                     beta     = cgm.exp0,   # Chipman et al. (2010)
                     gamma    = cgm.level0, # Chipman et al. (2010)
                     k        = sd.mu,      # prior SD on leaf means
                     num_tree = num.trees,
                     sigma_hat = 1,
                     tau_rate  = tau.rate)
    opts <- Opts(update_sigma = FALSE)

    forest.list[[ii]][[jj]] <- MakeForest(hypers, opts)
  }
}

#### --------------------------------------------------------------------- ####
#### ----------- BART / SV initialisation: idiosyncratic variance -------- ####
#### --------------------------------------------------------------------- ####
if(idioSV == "SV"){
  # Standard stochastic volatility (stochvol) per equation
  sv_priors_idio <- specify_priors(
    mu     = sv_normal(mean = 0, sd = 10),
    phi    = sv_beta(shape1 = 25, shape2 = 1.5),
    sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*0.005)),
    nu     = sv_infinity(),
    rho    = sv_constant(0))

  svdraw_idio <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
  q_sv <- list()
  for(mm in 1:M) q_sv[[mm]] <- svdraw_idio
}else if(idioSV == "HOMO"){
  a.sig <- 3
  b.sig <- 0.03
}
#### --------------------------------------------------------------------- ####
#### --------------- BART / SV initialisation: factor variance ----------- ####
#### --------------------------------------------------------------------- ####

if(FSV == "FHB"){
  # BART-heteroscedasticity on R volatility factors
  forest.FSV.list <- list()
  eta.star <- c(0, diff(rowMeans(y)))
  for(jj in seq_len(R)){
    hypers <- Hypers(X        = Z.fsv,
                     Y        = eta.star,
                     alpha    = 1,
                     beta     = 2,            # Chipman et al. (2010)
                     gamma    = 0.95,         # Chipman et al. (2010)
                     k        = sd.mu.sv.idio,
                     num_tree = num.trees.sv,
                     sigma_hat = 1,
                     tau_rate  = tau.rate)
    opts <- Opts(update_sigma = FALSE)
    forest.FSV.list[[jj]] <- MakeForest(hypers, opts)
  }
  ht.i <- rep(-3, N)

}else if(FSV == "FSV"){
  # Standard factor stochastic volatility (stochvol) on R factors
  sv_priors_FSV <- specify_priors(
    mu     = sv_normal(mean = 0, sd = 100),
    phi    = sv_beta(shape1 = 25, shape2 = 1.5),
    sigma2 = sv_gamma(shape = 0.5, rate = 5),
    nu     = sv_infinity(),
    rho    = sv_constant(0))

  svdraw_FSV <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
  h_sv <- list()
  for(rr in 1:R) h_sv[[rr]] <- svdraw_FSV
}

#### --------------------------------------------------------------------- ####
#### -------------------------- Starting values -------------------------- ####
#### --------------------------------------------------------------------- ####
I_N    <- diag(1, N)
V      <- matrix(cut.off.V, M, K)
Lambda <- array(1e-15, c(M, K, Q))

Lambda.sv    <- lambda.hs.sv <- nu.sv <- matrix(0, M, R)
tau.sv       <- zeta.sv <- rep(1, R)
b.draw       <- matrix(0, M, K)
trees.fit.mat <- array(1, c(N, M, Q))

ht      <- matrix(-2, N, R)
sig2.V  <- matrix(NA, N, M)
for(i in seq_len(M)){
  for(t in seq_len(N)) sig2.V[t,i] <- X[t,] %*% diag(V[i,]) %*% X[t,]
}

fit.star <- matrix(0, N, M)
bt.draw  <- array(0, c(N, M, K))

qt       <- matrix(-2, N, M)
sig2.idio   <- matrix(1, N, M)
sig2.idio.i <- rep(1, N)

sig2.t <- Sigma <- array(0, c(N, M, M))
for(nn in 1:N) sig2.t[nn,,] <- Sigma[nn,,] <- diag(sig2.idio[nn,])

commonalities.sv <- matrix(0, N, M)
var.share.sv     <- array(0, c(N, M, R))
rec.share.sv     <- matrix(NA, M, R)
rownames(rec.share.sv) <- var.names

f_draw <- extract(y, R)[[1]]
ind    <- rep(1, KK)

# HS prior draws: constant part
phi.cons       <- lambda.hs.cons <- nu.cons <- matrix(1, M, K)
zeta.hs.cons   <- tau.cons <- rep(1, M)
# HS prior draws: TVP loadings
phitau.mat     <- nu.mat <- lambda.hs.mat <- array(1, c(M, K, Q))
tau.mat        <- zeta.mat <- matrix(1, M, Q)

# Tree split counters
count.trees            <- array(0, c(M, Q, ncol(Z)))
duplicate.trees        <- matrix(0, M, Q)
count.trees.fsv        <- matrix(0, R, ncol(Z))
duplicate.trees.fsv    <- matrix(0, R, 1)
count.trees.sv.idio    <- matrix(0, M, ncol(Z))
duplicate.trees.sv.idio <- matrix(0, M, 1)

#### --------------------------------------------------------------------- ####
#### ------------------------- Storage arrays ---------------------------- ####
#### --------------------------------------------------------------------- ####

# In-sample
ind.store          <- matrix(0, nsave, 1)
fit.store          <- array(NA, c(nsave, N, M))
Lambda.store       <- array(NA, c(nsave, M, K, Q))
Lambda.sv.store    <- array(NA, c(nsave, M, R))
f.sv.store         <- array(NA, c(nsave, N, R))
phi.cons.store     <- array(NA, c(nsave, M, K))
phitau.store       <- array(NA, c(nsave, M, K, Q))
bt.store           <- array(NA, c(nsave, N, M, K))
b.store            <- array(NA, c(nsave, M, K))
trees.fit.store    <- array(NA, c(nsave, N, M, Q))
trees.count.store  <- array(NA, c(nsave, M, Q, ncol(Z)))
trees.duplicate.store <- array(NA, c(nsave, M, Q))
eht.store          <- array(NA, c(nsave, N, R))
sig2.idio.store    <- array(0,  c(nsave, N, M))
trees.count.fsv.store        <- array(NA, c(nsave, R, ncol(Z)))
trees.duplicate.fsv.store    <- array(NA, c(nsave, R))
trees.count.sv.idio.store    <- array(NA, c(nsave, M, M, ncol(Z)))
trees.duplicate.sv.idio.store <- array(NA, c(nsave, M, M))
sig2t.store        <- array(0,  c(nsave, N, M, M))
commonalities.store    <- array(0, c(nsave, M, K))
commonalities.sv.store <- array(0, c(nsave, N, M))
var.share.store        <- array(0, c(nsave, M, K, Q))
var.share.sv.store     <- array(0, c(nsave, N, M, R))

# Forecasts and IRFs
girf_store         <- array(0, c(nsave, N, n.irf, nhor, R+1))
PC_store           <- array(0, c(nsave, N, M.pr, nhor, R+1))
rec.share.sv.store <- array(0, c(nsave, M, R))
eig_store          <- matrix(0, nsave, N)

dimnames(girf_store)         <- list(1:nsave, 1:N, var.names, 1:nhor, 1:(R+1))
dimnames(PC_store)           <- list(1:nsave, 1:N, prices,    1:nhor, 1:(R+1))
dimnames(rec.share.sv.store) <- list(1:nsave, var.names, 1:R)

to.plot <- 25
start <- Sys.time()
#### --------------------------------------------------------------------- ####
#### ---------------------------------- MCMC ----------------------------- ####
#### --------------------------------------------------------------------- ####
pb <- txtProgressBar(min = 0, max = ntot, style = 3) #start progress bar
irep <- 1
for (irep in seq_len(ntot)){
start <- Sys.time()
#### --------------------------------------------------------------------- ####
#### ------------ Step 1: Joint draw of a_m, Lambda_m, gamma_m ----------- ####
#### --------------------------------------------------------------------- ####


for (mm in 1:M){
    # get full X with X_i = (X_t, X_t g_1t, ..., X_t g_Qt, f_1t, ..., f_Rt)'
    Xfull.i[,1:(KK-R)]  <- cbind(X, XX*kronecker(trees.fit.mat[,mm,], t(rep(1, K))))
    Xfull.i[,(KK-R+1):KK] <- f_draw

    # Normalize X and y
    normalizer <- 1/sqrt(sig2.V[,mm] + sig2.idio[,mm])
    y.hat <- y[,mm]*normalizer
    X.hat <- Xfull.i*normalizer
    
    # Extract prior elements
    V0inv.vec <- V0inv[mm,]
    
    # Posterior moments
    if(N > KK){
      diagV0.inv <- diag(V0inv.vec)
      V1 <- try(solve(crossprod(X.hat) + diagV0.inv), silent=TRUE) 
      if (is(V1,"try-error")) V1 <- ginv(crossprod(X.hat) + diagV0.inv)
      b1 <- V1%*%(crossprod(X.hat,y.hat))
      
      obs.draw <- try(b1 + t(chol(V1))%*%rnorm(KK), silent=TRUE)
      if (is(obs.draw, "try-error")) obs.draw <- mvrnorm(1, b1, V1)
    }else{
      # 1.) Draw from prior
      V0.vec <- 1/V0inv.vec
      b0 <- rep(0,KK)
      uu <- rnorm(KK, b0, sqrt(V0.vec))
      #2.)  Draw tau ~ N(0, I)
      dd <- rnorm(N, 0, 1)
      #3.) 
      vv <- X.hat%*%uu + dd
      #4.) Mz = (zz*Omega0*zz' + I)
      XV0 <- t(X.hat)*V0.vec
      Mz <- X.hat%*%XV0 + I_N 
      #5.) solve(Mz)
      ww <- Matrix::solve(Mz, (y.hat - vv))
      #6.) Obtain a draw for TVPs
      obs.draw <- uu + XV0%*%ww
    }

    #Split coefficients and loadings
      b.draw[mm,] <- obs.draw[1:K]
      Lambda.vec <- obs.draw[(K+1):(KK-R)]
      Lambda[mm,,] <- matrix(Lambda.vec, K, Q)
      Lambda.sv.vec <- obs.draw[(KK-R+1):(KK)]
      Lambda.sv[mm,] <- as.numeric(Lambda.sv.vec)
}
  
#### --------------------------------------------------------------------- ####
#### ----------------- Step 2: Horseshoe hyperparameters ----------------- ####
#### --------------------------------------------------------------------- ####

  
#### --------------------------------------------------------------------- ####
#### ------------------- Step 2.1: HS prior for gamma_m ------------------ ####
#### --------------------------------------------------------------------- ####

hss.sv_draw <- get.hs.mat(bdraw.mat = Lambda.sv,lambda.hs.mat = lambda.hs.sv,nu.hs.mat = nu.sv,tau.hs.vec = tau.sv,zeta.hs.vec = zeta.sv)
phitau.sv <- hss.sv_draw$psi # M x R
tau.sv <- hss.sv_draw$tau  #R x 1
lambda.hs.sv <-  hss.sv_draw$lambda
nu.sv <- hss.sv_draw$nu
nu.sv[nu.sv > 1] <- 1
zeta.sv <- hss.sv_draw$zeta
zeta.sv[zeta.sv > 0.1] <- 0.1
rownames(phitau.sv) <- ids
phitau.sv[phitau.sv < off.set.hs] <- off.set.hs
phitau.sv[phitau.sv > 0.1] <- 0.1

for (mm in 1:M){
    #### --------------------------------------------------------------------- ####
    #### --------------------- Step 2.2: HS prior for a_m -------------------- ####
    #### --------------------------------------------------------------------- ####

        hss_cons <- get.hs(bdraw = b.draw[mm,],lambda.hs = lambda.hs.cons[mm,] ,nu.hs = nu.cons[mm,],tau.hs = tau.cons[[mm]] ,zeta.hs = zeta.hs.cons[[mm]])
        phi.cons.mm <- hss_cons$psi # K x Q
        
        lambda.hs.cons[mm,] <- hss_cons$lambda
        tau.cons[[mm]] <- hss_cons$tau # Q x 1
        
        nu.cons[mm,] <- hss_cons$nu
        zeta.hs.cons[[mm]] <- hss_cons$zeta
      

      phi.cons.mm[phi.cons.mm < off.set.hs] <- off.set.hs
      phi.cons.mm[phi.cons.mm > 1e4] <- 1e4
      phi.cons[mm,] <- phi.cons.mm
    
  
    #### --------------------------------------------------------------------- ####
    #### ------------------ Step 2.3: HS prior for Lambda_m ------------------ ####
    #### --------------------------------------------------------------------- ####


        hss_draw <- get.hs.mat(bdraw.mat = matrix(Lambda[mm,,], K, Q),lambda.hs.mat = matrix(lambda.hs.mat[mm,,], K, Q),nu.hs.mat = matrix(nu.mat[mm,,], K, Q),tau.hs.vec = tau.mat[mm,],zeta.hs.vec = zeta.mat[mm,])
        phitau.mm <- hss_draw$psi # K x Q
        tau.mm <- hss_draw$tau # Q x 1
        
        lambda.hs.mat[mm,,] <- hss_draw$lambda
        nu.mat[mm,,] <- hss_draw$nu
        zeta.mat[mm,] <- hss_draw$zeta
      
      
      phitau.mm[phitau.mm < off.set.hs] <- off.set.hs
      phitau.mm[phitau.mm > 1e4] <- 1e4
      phitau.mat[mm,,] <- phitau.mm
      tau.mat[mm,] <- tau.mm
      
    
   
    phi.mm <- c(phi.cons.mm, phitau.mm) 
    phi.mm <- c(phi.mm, phitau.sv[mm,])
    
    V0inv[mm,] <- 1/phi.mm#Diagonal(KK, 1/phi.mm)#diag(1/phi.mm)
}     
  
  
#### --------------------------------------------------------------------- ####
#### --------------------- Step 3: Sample BART trees --------------------- ####
#### --------------------------------------------------------------------- ####

 
fit.vals <- matrix(0, N, M) #Fit of TVPs x_t'b_t
state.mean <- matrix(NA, N, K) #Coefficients b_t

for (mm in 1:M){
    y.mm <- y[,mm] - f_draw%*%t(Lambda.sv)[,mm]
    b.draw.mm <- b.draw[mm,]  
    sig2.mm <- sig2.idio[,mm]
    sig2.V.mm <- sig2.V[,mm]  
  
      
    for (qq in 1:Q){
    #Construct lambda(qt) and lambda(jt) for j != q
      lambda.qt <- X%*%Lambda[mm,,qq]
      lambda.qstar.t <- X %*% Lambda[mm,, -qq]
      if (Q==1) trees.fit.mat.mm <- matrix(trees.fit.mat[,mm,], N, Q) else trees.fit.mat.mm <- trees.fit.mat[,mm,]
      fit.star <- rowSums(lambda.qstar.t*trees.fit.mat.mm[, -qq, drop = F])
      #Compute differences between y and the fitted values of the remaining trees and then normalize by lambda(qt)
      y.tilde <- (y.mm - X%*%b.draw.mm - fit.star) / lambda.qt # CHK CHKG 
      var.tilde <- as.numeric(sig2.mm + sig2.V.mm)/lambda.qt^2  
      
      # Run an independent BART model for the i-th conditional variance factor
      rep.i <- as.numeric(forest.list[[mm]][[qq]]$do_gibbs_weighted(Z, y.tilde, 1/var.tilde, Z, 1))
      trees.fit.mat[,mm,qq] <- rep.i   
      
      fit.vals[,mm] <- fit.vals[,mm] + lambda.qt * trees.fit.mat[,mm,qq]
      count.trees[mm,qq,]    <- as.numeric(forest.list[[mm]][[qq]]$get_counts())/num.trees
      duplicate.trees[mm,qq] <- length(trees.fit.mat[,mm,qq][!duplicated(trees.fit.mat[,mm,qq])])/num.trees
        
    }
    
      
    
    #### --------------------------------------------------------------------- ####
    #### ----------------------- Step 3.1: Sample TVPs ----------------------- ####
    #### --------------------------------------------------------------------- ####

    cond.mean.trees <- trees.fit.mat[,mm,]%*%t(Lambda[mm,,])
    state.mean <- matrix(NA, N, K)
    
    V.diag <- V[mm,]
    sig2.V[,mm] <- colSums(t(X^2)*V.diag)
  
    for (n in seq_len(N)){
      Qn <- 1/(sum(X[n, ]*V.diag*X[n, ]) + sig2.idio[n,mm]) #Qn <- 1/(tcrossprod(X[n, ], diag(V[mm,]))%*%X[n, ] + sig2.idio[n,mm])
      eps.n <- as.numeric(y.mm[n] - X[n,] %*% b.draw[mm,] - X[n, ]%*%Lambda[mm,,]%*%trees.fit.mat[n,mm,])
      eta.mean <- V.diag * X[n,] * eps.n * Qn #diag(V[mm,])%*%X[n, ]%*%Qn%*%eps.n
      state.mean[n,]  <- t(cond.mean.trees[n, ] + eta.mean) 
      
      VXn <- V.diag*X[n,,drop=F]
      Var.state <- diag(V[mm,]) - t(VXn)%*%Qn%*%VXn
      
      Var.state.diag <- diag(Var.state)
      Var.state[abs(Var.state) < 1e-6] <- 0 #Because off-diagonals are extremely close to zero --> zero it out to speed up sampling
      diag(Var.state) <- Var.state.diag
      
      if (isDiagonal(Var.state)){
        bt.draw[n,mm,] <- t(state.mean[n,]+ rnorm(K, 0, sqrt(diag(Var.state)))) + b.draw[mm,]
      }else{
        Var.chol <- try(t(chol(Var.state)), silent = TRUE)
        if (is(Var.chol, "try-error")) bt.draw[n,mm,] <- mvrnorm(1, state.mean[n,], Var.state)  + b.draw[mm,] else bt.draw[n,mm,] <- t(state.mean[n,]+ Var.chol%*%rnorm(K)) + b.draw[mm,]
      }
    }

    #### --------------------------------------------------------------------- ####
    #### --------- Step 3.2: Sample process innovation variances V_m --------- ####
    #### --------------------------------------------------------------------- ####

    shocks.states <- bt.draw[,mm,] -  state.mean
    SSE.states <- apply(shocks.states^2,2,sum)
    
    for (ll in seq_len(K)) V[mm,ll] <- GIGrvg::rgig(1, 1/2 - N/2, SSE.states[[ll]], 1/(2*1e-8))
    V[V > cut.off.V] <- cut.off.V
    
      
}
    

#### --------------------------------------------------------------------- ####
#### ----------- Step 4: Sample idiosyncratic log-volatilities ----------- ####
#### --------------------------------------------------------------------- ####

shocks.obs <- y - X%*%t(b.draw) - fit.vals
shocks.idio <- shocks.obs - f_draw%*%t(Lambda.sv)

# Step 9a: SV per equation (stochvol) 
if(idioSV == "SV"){
  for(mm in seq_len(M)){
    eta.star <- shocks.idio[,mm]  
    svdraw_idio <- svsample_fast_cpp(eta.star, startpara = q_sv[[mm]], startlatent = qt[,mm], priorspec = sv_priors_idio)
    svdraw_idio[c("mu", "phi", "sigma", "nu", "rho")] <- as.list(svdraw_idio$para[, c("mu", "phi", "sigma", "nu", "rho")])
    
    q_sv[[mm]] <- svdraw_idio
    sv_latent_idio <- svdraw_idio$latent
    qt.i <- as.numeric(sv_latent_idio)
    qt.i[qt.i > off.set.SV] <- off.set.SV
    qt.i[qt.i < -10] <- -10
    qt[,mm] <- qt.i
    
    sig2.idio[,mm] <-  exp(qt.i)
  }   
}else if(idioSV == "HOMO"){
   for (mm in 1:M){
      sig2.idio.i <- 1/rgamma(1, N/2 +  a.sig, sum(shocks.idio[,mm]^2)/2 + b.sig)
      sig2.idio[,mm] <- sig2.idio.i
      qt[,mm] <- log(sig2.idio.i)
   }    
}

  #### --------------------------------------------------------------------- ####
  #### ---------- Step 5: Sample factor log-volatilities via BART ---------- ####
  #### --------------------------------------------------------------------- ####
 
    for(rr in seq_len(R)){
      eta.star <- f_draw[,rr]  
      eta.off <- 1e-3
      
      eta.star <- log(eta.star^2 + eta.off)
      
      # Sample the mixture components (previous draw of ht)
      z <- sapply(eta.star-ht[,rr], ncind, m_st, sqrt(v_st2), q_st)
      
      # Subset mean and variances to the sampled mixture components; (n x p) matrices
      m_st_all = matrix(m_st[z],N, 1)
      v_st2_all = matrix(v_st2[z], N, 1)
      
      eta.star <- matrix((eta.star - m_st_all)/sqrt(v_st2_all))
      
      #Do BART for conditional variances
      ht.i <- as.numeric(forest.FSV.list[[rr]]$do_gibbs_weighted(Z.fsv, eta.star, 1/v_st2_all, Z.fsv, 1))
      ht.i.mean <- mean(ht.i)
      ht.i.sd   <- sd(ht.i)
      if(fsv.std && ht.i.sd != 0) ht.i <- (ht.i - ht.i.mean)/ht.i.sd else ht.i <- (ht.i - ht.i.mean) 
      ht.i[ht.i > off.set.SV] <- off.set.SV
      ht.i[ht.i < -20] <- -20
      ht[,rr] <- ht.i
    
      count.trees.fsv[rr,]     <- as.numeric(forest.FSV.list[[rr]]$get_counts())/num.trees.sv
      duplicate.trees.fsv[rr,] <- length(ht.i[!duplicated(ht.i)])/num.trees.sv
      
    }
    


  for (t in seq_len(N)){
      if(R == 1){
          sig2.t[t,,] <- exp(ht[t,])*Lambda.sv%*%t(Lambda.sv) + diag(sig2.idio[t,])   
      }else{
          sig2.t[t,,] <- Lambda.sv%*%diag(exp(ht[t,]))%*%t(Lambda.sv) + diag(sig2.idio[t,])  
          commonalities.sv[t,] <- diag(Lambda.sv%*%diag(exp(ht[t,]))%*%t(Lambda.sv))/diag(sig2.t[t,,])
          for(fj in 1:R){
            var.share.sv[t,,fj] <- diag(Lambda.sv[,fj,drop = F]%*%t(Lambda.sv[,fj, drop = F]))*exp(ht[t,fj])/diag(sig2.t[t,,])
          }
      }
        
  }
   
  if(R > 1){
      for(rr in 1:R) rec.share.sv[,rr] <- colMeans(var.share.sv[,,rr])
      permut.fac <- order(rec.share.sv[var.names[1],], decreasing = TRUE)
      Lambda.sv.perm <- Lambda.sv[,permut.fac]
      var.share.sv <- var.share.sv[,,permut.fac]
  }else{
    permut.fac <- 1
    Lambda.sv.perm <- Lambda.sv
  }  
   
  #### --------------------------------------------------------------------- ####
  #### ----------------- Step 6: Sample latent factors q_t ----------------- ####
  #### --------------------------------------------------------------------- ####
 
  err_f <- 0
  for (t in 1:N){
    if(M > 1){
      if (R>1){
      Ht <- exp(diag(ht[t,]))
      Qt <- Lambda.sv%*%Ht%*%t(Lambda.sv)+  diag(sig2.idio[t,])
      At <- Ht%*%t(Lambda.sv)%*%solve(Qt)
      f_mean <- At%*%t(shocks.obs[t,,drop = F])
      f_var <- Ht-At%*%Qt%*%t(At)
      f_chol <- try(t(chol(f_var)),silent=TRUE)
      if (is(f_chol,"try-error")) f_draw.t <- try(MASS::mvrnorm(1,f_mean,f_var, tol = 1e-2), silent = TRUE) else f_draw.t <- t(f_mean+f_chol%*%rnorm(R))
      if (is(f_draw.t,"try-error")) err_f <- err_f + 1 else f_draw[t,] <- f_draw.t
      
      }else if(R==1){
      Ht <- exp(ht[t,])
      Qt <- Lambda.sv%*%Ht%*%t(Lambda.sv)+  diag(sig2.idio[t,])
      At <- Ht%*%t(Lambda.sv)%*%solve(Qt)
      f_mean <- as.numeric(At%*%t(shocks.obs[t,,drop = F]))
      f_var <- Ht-At%*%Qt%*%t(At)
      f_chol <- sqrt(f_var)
      f_draw[t,] <- rnorm(1, f_mean, f_chol)  
      }
    
  }else{
    if (R>1){
      Ht <- exp(diag(ht[t,]))
      Qt <- Lambda.sv%*%Ht%*%t(Lambda.sv)+  sig2.idio[t,]
      At <- Ht%*%t(Lambda.sv)%*%solve(Qt)
      f_mean <- At%*%t(shocks.obs[t,,drop = F])
      f_var <- Ht-At%*%Qt%*%t(At)
      f_chol <- try(t(chol(f_var)),silent=TRUE)
      if (is(f_chol,"try-error")) f_draw.t <- try(MASS::mvrnorm(1,f_mean,f_var, tol = 1e-2), silent = TRUE) else f_draw.t <- t(f_mean+f_chol%*%rnorm(R))
      if (is(f_draw.t,"try-error")) err_f <- err_f + 1 else f_draw[t,] <- f_draw.t
      
    }else if(R==1){
      Ht <- exp(ht[t,])
      Qt <- Lambda.sv%*%Ht%*%t(Lambda.sv)+  sig2.idio[t,]
      At <- Ht%*%t(Lambda.sv)%*%solve(Qt)
      f_mean <- as.numeric(At%*%t(shocks.obs[t,,drop = F]))
      f_var <- Ht-At%*%Qt%*%t(At)
      f_chol <- sqrt(f_var)
      f_draw[t,] <- rnorm(1, f_mean, f_chol)  
    }
   }
  }
  
  f_draw <- apply(as.matrix(f_draw), 2, function(x) (x- mean(x)))


fit <- matrix(0, N, M)
for (n in seq_len(N)) fit[n,] <- X[n,]%*%t(bt.draw[n,,]) + f_draw[n,]%*%t(Lambda.sv)

if(irep %in% save.set)
{
    save.ind <- save.ind + 1
    
    ind.store[save.ind,] <- sum(ind)/KK
    fit.store[save.ind,,] <- fit 
    Lambda.store[save.ind,,,] <- Lambda
    Lambda.sv.store[save.ind,,] <- Lambda.sv
    f.sv.store[save.ind,,] <- f_draw
    sig2t.store[save.ind,,,] <- sig2.t 
    var.share.sv.store[save.ind,,,] <- var.share.sv
    rec.share.sv.store[save.ind,,] <- rec.share.sv
    commonalities.sv.store[save.ind,,] <- commonalities.sv
    
    phi.cons.store[save.ind,,] <- phi.cons
    phitau.store[save.ind,,,] <- phitau.mat
    bt.store[save.ind,,,] <- bt.draw 
    b.store[save.ind,,] <- b.draw 
    eht.store[save.ind,,] <- exp(ht)
    sig2.idio.store[save.ind,,] <- sig2.idio
    
    trees.fit.store[save.ind,,,] <- trees.fit.mat
    trees.count.store[save.ind,,,] <- count.trees
    trees.duplicate.store[save.ind,,] <- sort(duplicate.trees, decreasing = TRUE)
    trees.count.fsv.store[save.ind,,] <- count.trees.fsv
    trees.duplicate.fsv.store[save.ind,] <- sort(duplicate.trees.fsv, decreasing = TRUE)
    trees.count.sv.idio.store[save.ind,,,] <- count.trees.sv.idio
    trees.duplicate.sv.idio.store[save.ind,,] <- sort(duplicate.trees.sv.idio, decreasing = TRUE)

    if(Q > 1){
      for (jj in seq_len(M)){
        commonalities.store[save.ind, jj,] <- diag(Lambda[jj,,] %*% diag(apply(trees.fit.mat[,jj,],2,var)) %*% t(Lambda[jj,,])) / diag(Lambda[jj,,] %*% diag(apply(trees.fit.mat[,jj,],2,var)) %*% t(Lambda[jj,,]) +diag(V[jj,]))  
        for(kk in seq_len(Q)){
          var.share.store[save.ind, jj,,kk] <-  diag(as.matrix(Lambda[jj,,kk]) %*% apply(as.matrix(trees.fit.mat[,jj,kk]),2,var) %*% t(as.matrix(Lambda[jj,,kk])))
        }
         
      }
      
    }else{
      for (jj in seq_len(M)){
        commonalities.store[save.ind, jj,] <- diag(as.matrix(Lambda[jj,,]) %*% apply(as.matrix(trees.fit.mat[,jj,]),2,var) %*% t(as.matrix(Lambda[jj,,])) /(as.matrix(Lambda[jj,,]) %*% apply(as.matrix(trees.fit.mat[,jj,]),2,var) %*% t(Lambda[jj,,]) +diag(V[jj,])))
    
      }
    }
    

    # Eigenvalues — always computed (used for post-sampling stability filter)
    for(nn in 1:N){
      comp <- get.companion(A = t(bt.draw[nn,,1:(M*p)]), m = M, p = p, cons = FALSE)
      eig_store[save.ind,nn] <- max(abs(eigen(comp$MM)$values))
    }

    # Compute GIRFs 
    if(irf.on){
      girf <- array(0,dim=c(N, M, nhor, R +1))
      dimnames(girf) <- list(1:N, var.names, 1:nhor, 1:(R+1))
      PC <- array(0, c(N, M.pr, nhor,  R + 1))
      dimnames(PC) <- list(1:N, prices, 1:nhor, 1:(R+1))
      
      for(nn in 1:N){
        # Baseline initial conditions at time t
        Z0.girf0 <- y[nn,]/Z.iqr[nn,] - Z.min[nn,]/Z.iqr[nn,]
        Z0.girf  <- matrix(Z0.girf0, M, 1)
        X0.girf  <- matrix(c(y[nn,], X[nn,1:(M*(p-1))], 1), K, 1)
        
        # Impact response (h=0): construct shock matrix and apply scaling
        shk.mat <- matrix(0, M, R+1)
        shock.chol <- t(chol(sig2.t[nn,,]))
        colnames(shock.chol) <- rownames(shock.chol) <- rownames(shk.mat) <- rownames(Z0.girf) <- var.names
        if(shock.ident == "external") shock.chol[] <- bt.draw[nn,,K]
        rownames(Lambda.sv.perm) <- var.names
        if(shk.sc == "prd-wise"){
          girf[nn,,1,1:R] <- shk.mat[,1:R]   <- Lambda.sv.perm/t(matrix(Lambda.sv.perm[shk.max.var,], R, M))
          girf[nn,,1,R+1] <- shk.mat[,(R+1)] <- shock.chol[,shk.var]/shock.chol[shk.var, shk.var] 
        }else if(shk.sc == "avg"){
          girf[nn,,1,1:R] <- shk.mat[,1:R]   <- Lambda.sv.perm*colMeans(exp(ht/2))
          girf[nn,,1,R+1] <- shk.mat[,(R+1)] <- shock.chol[,shk.var]
        }else if(shk.sc == "none"){
          girf[nn,,1,1:R] <- shk.mat[,1:R]   <- Lambda.sv.perm*exp(ht[nn,]/2) 
          girf[nn,,1,R+1] <- shk.mat[,(R+1)] <- shock.chol[,shk.var]
        }
        PC[nn,,1,] <- girf[nn,prices,1,]/girf[nn,"UNRATE",1,]
        
        # Shocked initial conditions for each of the R+1 shocks
        Z1.girf <- matrix(Z0.girf, M, R+1) + shk.mat/matrix(Z.iqr[nn,], M, R+1)
        X1.girf <- list()
        for(ss in 1:(R+1)) X1.girf[[ss]] <- X0.girf + c(shk.mat[,ss], rep(0, M*(p-1)), 0) 
            
        # Propagate baseline and shocked paths forward horizon by horizon
        for(hh in 2:nhor){
          path0.girf      <- matrix(0, M, 1)
          trees0.girf.mat <- matrix(0, M, Q)
          bt0.girf        <- matrix(0, M, K)
          
          path1.girf <- matrix(0, M, R+1)
          trees1.girf.mat <- array(0, c(R+1, M, Q))
          bt1.girf <- array(0, c(R+1, M, K))
          
          for(mm in 1:M){
            for(qq in 1:Q){
              trees0.girf.mat[mm,qq] <-  forest.list[[mm]][[qq]]$do_predict(t(Z0.girf))
              trees1.girf.mat[,mm,qq] <- forest.list[[mm]][[qq]]$do_predict(t(Z1.girf))
            }
            # Predict TVP coefficients for baseline and shocked paths
            bt0.girf[mm,]  <- trees0.girf.mat[mm,]%*%t(Lambda[mm,,]) + b.draw[mm,]
            bt1.girf[,mm,] <- trees1.girf.mat[,mm,]%*%t(Lambda[mm,,]) + t(matrix(b.draw[mm,], K, R+1))  
          }
          
          # Advance baseline one step
          path0.girf <- bt0.girf%*%X0.girf
          X0.girf <- matrix(c(path0.girf, X0.girf[1:(M*(p-1))], 1))
          
          # Advance each shocked path; GIRF = shocked minus baseline
          for(ss in 1:(R+1)){
            bt1.ss <- bt1.girf[ss,,]
            path1.girf[,ss] <- bt1.ss%*%X1.girf[[ss]]
            girf[nn,,hh,ss] <- path1.girf[,ss] - path0.girf
            X1.girf[[ss]] <- matrix(c(path1.girf[,ss],X1.girf[[ss]][1:(M*(p-1))], 1))
          }
          # Re-normalise Z paths for BART prediction at next horizon
          Z0.girf <- path0.girf/Z.iqr[nn,] - Z.min[nn,]/Z.iqr[nn,]
          Z1.girf <- path1.girf/matrix(Z.iqr[nn,], M, R+1) - matrix(Z.min[nn,]/Z.iqr[nn,], M, R+1)
          
          PC[nn,,hh,] <- girf[nn,prices,hh,]/girf[nn,"UNRATE",hh,]
        }
      }
      girf_store[save.ind,,,,] <- girf
      PC_store[save.ind,,,,] <- PC

    }
    
}
  
end <- Sys.time()


setTxtProgressBar(pb, irep)

}

end <- Sys.time()
# Time in minutes
time.min <- (ts(end)-ts(start))/60

return(list(fit.store = fit.store, Lambda.store = Lambda.store, Lambda.sv.store = Lambda.sv.store, f.sv.store = f.sv.store, sig2t.store = sig2t.store, sig2.idio.store = sig2.idio.store, commonalities.sv.store = commonalities.sv.store, var.share.sv.store = var.share.sv.store, phi.cons.store = phi.cons.store,  phitau.store = phitau.store, bt.store = bt.store, b.store = b.store, eht.store = eht.store, trees.fit.store = trees.fit.store, trees.count.store = trees.count.store, trees.duplicate.store = trees.duplicate.store, trees.count.fsv.store = trees.count.fsv.store, trees.duplicate.fsv.store = trees.duplicate.fsv.store, trees.count.sv.idio.store = trees.count.sv.idio.store, trees.duplicate.sv.idio.store = trees.duplicate.sv.idio.store,M = M, KK = KK, time.min = time.min, ntot = ntot, commonalities.store = commonalities.store, ind.store = ind.store, var.share.store = var.share.store, girf_store = girf_store, PC_store = PC_store, rec.share.sv.store = rec.share.sv.store, eig_store = eig_store))
  
}


