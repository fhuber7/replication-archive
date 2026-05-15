require(compiler)

fcst.tvp <- function(Xraw.full,p=2,nsave=15000,nburn=15000,thin=1,h=80,fhorz=4,TVS){
  
  Xraw <- Xraw.full[1:(nrow(Xraw.full)-h),]
  
  class(Xraw) <- "numeric"
  
  M <- ncol(Xraw)
  Y <- Xraw
  #Y <- apply(Y,2,function(x) x-mean(x)) #demean
  #Create lagged Y matrix
  Xlag <- mlag(Y,p)
  Y <- Y[(p+1):nrow(Y),,drop=FALSE]
  X <- cbind(1,Xlag[(p+1):nrow(Xlag),])
  T <- nrow(X)
  K <- ncol(X)
  start <- Sys.time()
  model_VAR <- estimate_tvp(Y,X,save=nsave,burn=nburn,p=p,sv_on = TRUE,thin = thin,priorbtheta = list(B_1=3,B_2=0.03,kappa0=1e-8),priormu=c(0,10),h0prior="stationary", grid.length = 150, thrsh.pct = 0.2,thrsh.pct.high = 1,TVS=TVS,CPU=1) 
  end <- Sys.time()
  print(end-start)
  #Compute draw from the predictive density
  pred_store <- array(0,c(nsave,M,fhorz))
  for (irep in 1:nsave){
    #Specific for a draw from the VAR posterior
    #irep denotes the irep'th draw 
    SIGMA <- model_VAR$VAR_coeff$S_post[T,,,irep]
    A_draw <- model_VAR$VAR_coeff$A_post[T,,,irep]
    zt <- X[T,]#Xraw[nrow(Xraw)-h,]
    
    z1 <- zt
    Mean00 <- zt
    Sigma00 <- matrix(0,K,K)
    Sigma00[1:M,1:M] <- SIGMA
    Cm <- diag(M)
    varndxv <- c(M,1,p)
    SIGMAt <- SIGMA
    companion <- get_companion(A_draw,varndxv)
    Jm <- companion$Jm;Mm <- companion$MM
    # if (max(abs(Re(eigen(Mm)))>1.05)) next
    for (ih in 1:fhorz){
      z1 <- Mm%*%z1
      cholSig <- try(t(chol(Sigma00[1:M,1:M])),silent=TRUE)
      if (is(cholSig,"try-error")){
        yf <- rmvnorm(1,z1[1:M],Sigma00[1:M,1:M])
      }else{
        yf <- z1[1:M]+cholSig%*%rnorm(M,0,1)
      }
      #first and second moments
      Mean00 <- Mm%*%Mean00
      Sigma00 <- Mm%*%Sigma00%*%t(Mm)+Jm%*%SIGMAt%*%t(Jm)
      meanyt <- Cm%*%Mean00[1:M]
      varyt <- Cm%*%Sigma00[1:M,1:M]%*%t(Cm)
      pred_store[irep,,ih] <- yf
    }
  }
  return(pred_store)
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
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

int_posterior <- function(d=d_prop,Achg1=Achg[jj,],sqrttheta11=sqrttheta1[jj,jj],sqrttheta21=sqrttheta2[jj,jj]){
  
  #d=d_prop[kk];Achg1=Achg[jj,];sqrttheta11=sqrttheta1[jj,jj];sqrttheta21=sqrttheta2[jj,jj]
  #if (d<1e-08) d <- 1e-08
  D_i <- (abs(Achg1)>d)*1
  
  b_01 <- Achg1[D_i==1]
  b_02 <- Achg1[D_i==0]
  A_1 <- sum(dnorm(b_01,0,sqrt(sqrttheta11),log=TRUE))
  A_2 <- sum(dnorm(b_02,0,sqrt(sqrttheta21),log=TRUE))
  
  post_1 <- A_1+A_2#+dinvgamma(d,1.5,1)
  return(post_1)
}


Mu1b <- Vectorize(int_posterior, "d") 


get_companion <- function(Beta_,varndxv){
   #Beta_ <- A_draw
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  if (nd==1){
  MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn+1)),c(matrix(0,1,nn*nl),1))
  }else{
    MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn)))
  }
  return(list(MM=MM,Jm=Jm))
}

