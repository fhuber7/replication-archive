###----------------------------------- Functions ----------------------------------------
#z = draw.indicators(res = ystar-h_prev, nmix = list(m = m_st, v = v_st2, p = q))
ncind <- function(y,mu,sig,q){
  sample(1:length(q),
         size = 1,
         prob = q*dnorm(y,mu,sig))
}

var.spec <- function(B.comp,p, M, cons = cons, Sig = H) {
  K <- M*p  + cons
  Id <- diag(1,nrow=K,ncol=K) # identity matrix 
  J <- matrix(0,K,M)
  J[1:M,1:M] <- diag(1,M)  
  A <- (Id - B.comp)
  sp <- try(t(J)%*%(solve(A)%*%Sig%*%t(solve(A)))%*%J, silent=TRUE) 
  if (is(sp,"try-error")) sp <- t(J)%*%(ginv(A)%*%Sig%*%t(ginv(A)))%*%J
  #sp <- t(J)%*%(ginv(A)%*%tcrossprod(Sig)%*%t(ginv(A)))%*%J
  return(sp) 
} 

get.post <- function(store,dims){
  mean <- apply(store,2:dims,median)
  lo <- apply(store,2:dims,quantile,0.16)
  hi <- apply(store,2:dims,quantile,0.84)
  lo1 <- apply(store,2:dims,quantile,0.05)
  hi1 <- apply(store,2:dims,quantile,0.95)
  return(list(mean,lo,hi,lo1,hi1))
}   

mlag <- function(X,lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)  
}

mlag0 <- function(X,lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[p:Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii+1),(1:N)]
  }
  return(Xlag)  
}


gammacoef <- function(mode, sd){
  k.shape <- (2+mode^2/sd^2+sqrt((4+mode^2/sd^2)*mode^2/sd^2))/2
  theta.scale <- sqrt(sd^2/k.shape)
  return(data.frame(shape = k.shape, rate = 1/theta.scale))
}

get_V <- function(a_bar_1,a_bar_2,ind1,sigma_sq1,p1,set.cons,M1,K1,M2=1){
  V_i <- matrix(0,K1,M2)
  
  #this double loop fills the prior covariance matrix
  for (i in 1:M2){ #for each i equation
    for (j in 1:K1){ #for each variable on the rhs
      if (set.cons==1) {
        if (j==K1) {
          V_i[j,i] <- a_bar_2*sigma_sq1[i,1] #variance on constant,trend, dummies and on ex variables         
        }    
        else if (any(j==ind1[i,])) {
          ll <- which(ind1[i,]==j)
          
          V_i[j,i] <- a_bar_1/(ll^2) #variance on own lags
        }else{
          for (kj in 1:M1){
            if (any(j==ind1[kj,]))
            {
              ll <- kj
            }  
          } 
          V_i[j,i] <- (a_bar_2*sigma_sq1[i,1])/((ll^2)*sigma_sq1[ll,1])
        }
      }       
      else{
        if (any(j==ind1[i,])) {
          ll <- which(ind1[i,]==j)
          
          V_i[j,i] <- a_bar_1/(ll^2) #variance on own lags
        } 
        else
        {
          for (kj in 1:M1){
            if (any(j==ind1[kj,]))
            {
              ll <- kj
            }  
          } 
          V_i[j,i] <- (a_bar_2*sigma_sq1[i,1])/((ll^2)*sigma_sq1[ll,1])
        }
      }
    }
  }  
  return(V_i)
}


# get_V_gprior without lag penalisation and a_bar_1 = a_bar_2
get_V_gprior <- function(a_bar_1, ind1,sigma_sq1,p1,set.cons,M1,K1,M2=1){
  V_i <- matrix(0,K1,M2)
  
  #this double loop fills the prior covariance matrix
  for (i in 1:M2){ #for each i equation
    for (j in 1:K1){ #for each variable on the rhs
      if (set.cons==1) {
        if (j==K1) {
          V_i[j,i] <- a_bar_1*sigma_sq1[i,1] #variance on constant,trend, dummies and on ex variables         
        }    
        else if (any(j==ind1[i,])) {
          ll <- which(ind1[i,]==j)
          
          V_i[j,i] <- a_bar_1 #variance on own lags
        }else{
          for (kj in 1:M1){
            if (any(j==ind1[kj,]))
            {
              ll <- kj
            }  
          } 
          V_i[j,i] <- (a_bar_1*sigma_sq1[i,1])/(sigma_sq1[ll,1])
        }
      }       
      else{
        if (any(j==ind1[i,])) {
          ll <- which(ind1[i,]==j)
          
          V_i[j,i] <- a_bar_1#variance on own lags
        } 
        else
        {
          for (kj in 1:M1){
            if (any(j==ind1[kj,]))
            {
              ll <- kj
            }  
          } 
          V_i[j,i] <- (a_bar_1*sigma_sq1[i,1])/(sigma_sq1[ll,1])
        }
      }
    }
  }  
  return(V_i)
}


hamiltonfilter <- function(beta1,beta2,sigma1,sigma2,p,q,y,x){
  T <- nrow(y)
  isig1 <- solve(sigma1)
  isig2 <- solve(sigma2)
  dsig1 <- sigma1^2
  dsig2 <- sigma2^2
  #Initialise the Filter
  pr_tr <- rbind(c(p[1],1-q[1]),c(1-p[1],q[1]))
  A <- rbind(diag(2)-pr_tr,matrix(1,1,2))
  EN <- matrix(c(0,0,1),3,1)
  ett <- ginv(crossprod(A))%*%t(A)%*%EN
  if (any(is.na(ett))){
    ett <- c(0.5,0.5)
  }
  #Forward Filter
  lik <- 0
  fprob <- matrix(0,T,2)
  for (it in 1:T){
    pr_tr <- rbind(c(p[it],1-q[it]),c(1-p[it],q[it]))
    em1 <- (y[it,] - x[it,]%*%beta1)
    em2 <- (y[it,] - x[it,]%*%beta2)
    
    neta1 <- (1/sqrt(dsig1))*exp(-0.5*(em1%*%isig1%*%t(em1))) 
    neta2 <- (1/sqrt(dsig2))*exp(-0.5*(em2%*%isig2%*%t(em2))) 
    
    neta1[neta1 < 1e-30] <- 1e-30
    neta2[neta2 < 1e-30] <- 1e-30
    neta1[neta1 > 1e30] <- 1e30
    neta2[neta2 > 1e30] <- 1e30
    
    #  if ((neta1 && neta2)==0) neta1 <- 0.1
    
    #Kim and Nelson Algorithm
    #  if (any(is.na(ett))) ett <- c(0.5,0.5)
    
    ett1 <- ett*rbind(neta1,neta2)
    
    fit <- sum(ett1)
    ett <- (pr_tr%*%ett1)/fit
    
    fprob[it,] <- t(ett1/fit)
    
    if (fit>0 && !is.na(fit)){
      lik <- lik+log(fit)
    }else{
      lik <- lik-10
    }
  }
  
  return(list(fprob=fprob,lik=lik))
}

getS <- function(fprob,p,q,Ncrit,maxdraws){
  PROBLEM <- 0
  j <- 1
  chck <- -1
  while(j<maxdraws && chck<0){
    ST1 <- getst(fprob,p,q)
    ST <- ST1[[1]]
    
    check1 <- length(ST[!duplicated(ST)])==2
    T1 <- sum(ST==0)
    T2 <- sum(ST==1)
    check2 <- T1>=Ncrit
    check3 <-  T2>=Ncrit# T2 <30#CHGGGGGGGG
    if (check1+check2+check3==3){
      chck=10
    }else{
      j=j+1
    }
  }
  if (check1+check2+check3<3){
    PROBLEM=1
  }
  return(list(ST=ST,prob=PROBLEM,smooth=ST1[[2]]))
}
getst <- function(fprob,p,q){
  T <- nrow(fprob)
  ST <- matrix(0,T,1)
  filtered <- matrix(0,T,1)
  
  p00 <- fprob[T,1]
  p01 <- fprob[T,2]
  r <- runif(1,0,1)
  ST[T,1] <- (r>=(p00/(p00+p01)))
  #go backwards
  for (it in (T-1):1){
    pr_tr <- rbind(c(p[it],1-q[it]),c(1-p[it],q[it]))
    if (ST[it+1]==0){
      p00 <- pr_tr[1,1]*fprob[it,1]
      p01 <- pr_tr[1,2]*fprob[it,2]
    }else if (ST[it+1]==1){
      p00 <- pr_tr[2,1]*fprob[it,1]
      p01 <- pr_tr[2,2]*fprob[it,2]
    }
    #sample regime numbers
    r <- runif(1,0,1)
    alph <- (p00/(p00+p01))
    if (r<alph){
      ST[it] <- 0
    }else{
      ST[it] <- 1
    }
    filtered[it] <- alph
  }
  
  return(list(ST,filtered))
}


KF_R <- function(y, Z,Ht,Qtt, m, p, t, B0, V0){
  #Define everything calculated down there
  # B0 <- as.vector(B0)
  # V0 <- as.matrix(V0)
  # Qtt <- as.matrix(Qtt)
  # bt <- matrix(0, t,m)
  # Vt <- matrix(0, m^2, t)
  # R <- matrix(0, p,m)
  # H <- matrix(0, t*m, p)
  # cfe <- matrix(0, p,1)
  # yt <- matrix(0, p, 1)
  # f <- matrix(0, p,p)
  # inv_f <- matrix(0, p,p)
  # btt <- matrix(0, m, 1)
  # Vtt <- matrix(0, m, m)
  
  bp <- B0 #the prediction at time t=0 is the initial state
  Vp <- V0 #Same for the variance
  bt <- matrix(0,t,m) #Create vector that stores bt conditional on information up to time t
  Vt <- matrix(0,m^2,t) #Same for variances
  
  for (i in 1:t){
    R <- Ht[i,] #CHK LATER NOISE INNOVATION VARIANCE TO MEASUREMENT ERRORS
    if(ncol(Qtt) > 1) Qt <- diag(Qtt[i,]) else Qt <- Qtt[i,]
    H <- Z#[i,,drop = F]
    
    cfe <- y[i] - H%*%bp   # conditional forecast error
    f <- H%*%Vp%*%t(H) + R    # variance of the conditional forecast error
    inv_f <- try(t(H)%*%solve(f), silent = T)
    if(is(inv_f, "try-error")) inv_f <- t(H)%*%ginv(f)
    btt <- bp + Vp%*%inv_f%*%cfe  #updated mean estimate for btt Vp * inv_F is the Kalman gain
    Vtt <- Vp - Vp%*%inv_f%*%H%*%Vp #updated variance estimate for btt
    if (i < t){
      bp <- btt
      Vp <- Vtt + Qt
    }
    bt[i,] <- t(btt)
    Vt[,i] <- matrix(Vtt,m^2,1)
  }
  
  # draw the final value of the latent states using the moments obtained from the KF filters' terminal state
  bdraw <- matrix(0,t,m)
  
  bdraw.temp <- try(btt+t(chol(Vtt))%*%rnorm(nrow(Vtt)), silent=T)
  if (is(bdraw.temp, "try-error")) bdraw.temp <- mvrnorm(1, btt, Vtt+diag(1e-6,m))
  bdraw[t,] <- bdraw.temp
  
  #Now do the backward recurssions
  for (i in 1:(t-1)){
    if(ncol(Qtt) > 1) Qt <- diag(Qtt[t-1,]) else Qt <- Qtt[t-1,]
    bf <- t(bdraw[t-i+1,])
    btt <- t(bt[t-i,])
    Vtt <- matrix(Vt[,t-i,drop=FALSE],m,m)
    f <- Vtt + Qt
    
    inv_f <- try(Vtt%*%solve(f), silent = T)
    if(is(inv_f, "try-error")) inv_f <- Vtt%*%ginv(f)
    
    cfe <- bf - btt
    bmean <- t(btt) + inv_f%*%t(cfe)
    bvar <- Vtt - inv_f%*%Vtt
    
    bdraw.temp <- try(bmean+t(chol(bvar))%*%rnorm(nrow(bvar)), silent=T)
    if (is(bdraw.temp, "try-error")) bdraw.temp <- mvrnorm(1, bmean, bvar+diag(1e-6,m))
    
    bdraw[t-i,] <- bdraw.temp
  }
  
  return(bdraw)
}




extract <- function(data,k){
  t <- nrow(data);n <- ncol(data)
  xx <- crossprod(data)
  eigs <- eigen(xx)
  evec <- eigs$vectors;eval <- eigs$values
  
  lam <- sqrt(n)*evec[,1:k]
  fac <- data%*%lam/n
  
  return(list(fac,lam))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

get.mcg <- function(Lambda, tau, df = 3, a1 = 2.1, a2 = 3.1){
  # Lambda is K x Q
  # tau is Q x 1
  K <- nrow(Lambda)
  Q <- ncol(Lambda)
  phi <- phi.mat <- matrix(NA, K, Q)
  delta <- rep(1, Q)
  
  for (k in seq_len(Q)) {
    lambda2tau <- Lambda[,k]^2 * tau[k]
    phi.jh <- rgamma(K, df/2 + 0.5, df/2 + lambda2tau/2)
    phi[, k] <- phi.jh
    phi.mat[,k] <- phi.jh * Lambda[,k]^2
  }
  
  if(Q > 1){
    xi.sum <- apply(phi.mat, 2, sum)
    delta.1 <- rgamma(1, a1 + K * Q/2, 1 + 1/2 * sum(prod(delta[-1]) * xi.sum))
    delta[1] <- delta.1
    
    for (h in 2:Q) {
      delta[h] <- rgamma(1, a2 + K * (Q - h + 1)/2, 1 + 1/2 * sum(prod(delta[-h]) * xi.sum))
    } 
  }else{
    delta <- 1
  }
  
  tau <- cumsum(delta)
  
  ret <- list("psi"=(t(t(1/phi)/tau)),"tau"=tau,"phi"=phi)
  return(ret)
}


get.ng <- function(bdraw,psi,pr=0,d0=0.01,d1=0.01,th=0.1){
  k <- NROW(bdraw)
  
  # global shrinkage parameter
  d0_po <- d0 + th*k
  d1_po <- d1 + (th * mean(psi)*k)/2
  lambda2 <- rgamma(1,d0_po,d1_po)
  
  # local shrinkage scalings
  psi <- matrix(0,k,1)
  for(kk in 1:k){
    scale <- ifelse(bdraw[kk]^2<=1e-15, 1e-15, bdraw[kk]^2) # for stability
    psi[kk] <- GIGrvg::rgig(n=1, lambda=th-0.5, chi=scale, psi=lambda2*th)
  }
  return(psi)
}


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

get.hs.mat <- function(bdraw.mat,lambda.hs.mat,nu.hs.mat,tau.hs.vec,zeta.hs.vec){
  psi.hs.mat <- matrix(0, nrow(bdraw.mat), ncol(bdraw.mat))
  for (jj in 1:ncol(bdraw.mat)) {
    bdraw <- bdraw.mat[,jj]
    lambda.hs <- lambda.hs.mat[,jj]
    nu.hs <- nu.hs.mat[,jj]
    tau.hs <- tau.hs.vec[jj]
    zeta.hs <- zeta.hs.vec[jj]
 
    # standard HS function
    k <- length(bdraw)
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
    lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
    nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
    zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
    
    psi.hs <- (lambda.hs*tau.hs)
    
    psi.hs.mat[,jj] <-  psi.hs
    lambda.hs.mat[,jj] <- lambda.hs
    nu.hs.mat[,jj] <- nu.hs
    tau.hs.vec[jj] <- tau.hs
    zeta.hs.vec[jj] <- zeta.hs
  }

  ret <- list("psi"= psi.hs.mat,"lambda"=lambda.hs.mat,"tau"=tau.hs.vec,"nu"=nu.hs.mat,"zeta"=zeta.hs.vec)
  return(ret)
}


############## VAR functions

# get companion matrix
get.companion <- function(A,m,p,cons){
  nn <- m
  nd <- cons
  nl <- p
  
  nkk <- nn*nl+nd
  Jm <- matrix(0,nkk,nn)
  MM <- matrix(0,nkk,nkk)
  Jm[1:nn,1:nn] <- diag(nn)
  if (nd==1){
    MM <- rbind(t(A), cbind(diag((nl-1)*nn),matrix(0,(nl-1)*nn,nn+nd)),c(matrix(0,nd,nn*nl),nd))
  }else{
    MM <- rbind(t(A), cbind(diag((nl-1)*nn),matrix(0,(nl-1)*nn,nn+nd)))
  }
  return(list(MM=MM,Jm=Jm))
}



do.fcst <- function(Y,X, A.draw,Sig.t,pred.map = pred.map, fhorz, var.names, forc.var, T, M, M.map, K,p,cons){
  #Predict Sigma_t multi-steps ahead
  A.draw <- A.draw
  Sig.t <- Sig.t
  
  get.comp <- get.companion(A=A.draw,m, p, cons)
  MM <- get.comp$MM
  Jm <- get.comp$Jm
  eig.max <- try(max(Re(eigen(MM)$values)), silent = TRUE)
  if(is(eig.max, "try-error")) eig.max <- 100
  
  if(p > 1){
    if(cons) X0 <- c(Y[T,],X[T,(1:(K-M-1))],1) else X0 <- c(Y[T,],X[T,(1:(K-M))])
  }else{
    if(cons) X0 <- c(Y[T,],1) else X0 <- c(Y[T,])
  }  
  
  X0 <- as.matrix(X0)
  SS <- matrix(0,K,K)
  Sigma00 <- as.matrix(Sig.t)
  
  
  #Predict Sigma_t multi-steps ahead
  pred.cov <- array(NA,c(M.map,M.map,fhorz))
  pred.mean <- matrix(NA, M.map, fhorz)
  pred.draw <- matrix(NA,M.map,fhorz)
  rownames(pred.draw) <- rownames(pred.mean) <- var.names
  dimnames(pred.cov) <- list(var.names,var.names,1:fhorz)
  
  for (jj in 1:fhorz){
    # Mean
    X0 <- MM%*%X0 
    
    # Variance
    SS <- MM%*%SS%*%t(MM)+Jm%*%Sigma00%*%t(Jm)
    cholSig.0 <- try(t(chol(SS[1:M,1:M])),silent=TRUE)
    
    if (is(cholSig.0,"try-error")){
      y.forecast <- matrix(rmvnorm(1,mean=X0[1:M],sigma=SS[1:M,1:M]))
    }else{
      y.forecast <- X0[1:M]+cholSig.0%*%rnorm(M,0,1)
    }  
    pred.draw[,jj] <- pred.map%*%y.forecast
    pred.cov[,,jj] <- pred.map%*%SS[1:M,1:M]%*%t(pred.map)
    pred.mean[,jj] <- pred.map%*%X0[1:M]
  }
  
  pred.draw <- pred.draw[forc.var,]
  pred.cov <- pred.cov[forc.var,forc.var,]
  pred.mean <- pred.mean[forc.var,]
  
  return(list(pred.draw, pred.mean, pred.cov, eig.max))
}


# get IRF
impulsdtrf <- function(B,smat,nstep,time=FALSE,B1=NULL){
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
    if (time){
      if (it==8){
        B <- B1[, , seq(from=lags, to=1, by=-1)]  #reverse time index to allow matrix mult instead of loop
        B <- matrix(B1,nrow=nvar)
      }
    }
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




