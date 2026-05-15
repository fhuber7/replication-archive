#setwd("!inflation_GP_BNP/Quarterly")
#rm(list = ls())

#set.seed(123)

set.var <- "CPILFESL"
set.var <- "CPIAUCSL"

foldername <- paste0("Results_", set.var)
dir.create(foldername, showWarnings = FALSE)


# adapted mlag function
mlag <- function(X,lag)
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

Rcpp::sourceCpp("BAKRGibbs.cpp") # This needs Rcpp installed
source("gpsubspace_function.R") # Sources the main function
source("ucsv.R") # sources the UCSV version
source("BART_function.R") # Sources BART versions

require(bvarsv)
require(Rcpp)
require(stochvol)
require(flexmix)
# Preliminaries ------------------------------------------------------------------
length.hold.out <- 185 # end of in-sample
h <- c(1, 4)#,3, 12                                 # 1-month, 1-quarter, 1-year horizon
p <- 1                                 # p lags
cons <- TRUE                           # constant in X
nburn <- 2000
ntot <- 4000

# Create grid ------------------------------------------------------------------
model.est <-  c("BART","reg", "gp", "subspace","UCSV")#, "UCSV"
tau.grid <- c(1e-4, 1e+4)
G.grid <- c("mixSV")#, "25", "1",  "SV"
data.type <- c("small","full","AR")
time.grid <- seq(1, length.hold.out)

combi <- expand.grid(h = h, time = time.grid, model.est, G.grid, data=data.type, stringsAsFactors = FALSE)
run<- 1
#run <-  as.integer(Sys.getenv("SGE_TASK_ID"))

dir.create("errors", showWarnings = FALSE)
# Hold-out/horizions combination
slct.run <- as.character(combi[run,])
h.slct <- combi[run, "time"]
h <- combi[run, "h"]
data.type <- combi[run, "data"]
mod.slct <- combi[run, 3]
G.slct <- combi[run, 4]

dir.create("errors", showWarnings =FALSE)
filename <- paste0(foldername, "/", "res_time_",h.slct, "_fhorz_", h, "_model_", mod.slct, "_G_", G.slct,"_data_",data.type,".RData")
outputFile <- file(paste0("errors/", "res_time_",h.slct, "_fhorz_", h, "_model_", mod.slct, "_G_", G.slct,"_data_",data.type,".txt"))

if(!file.exists(filename)){
  tryCatch({
# Define X and y------------------------------------------------------------------------
load(paste0("Data/makroUS_Q_h_", h,"_",set.var,".rda"))  
  
y <- as.matrix(Xmat[,1])
y <- (y-mean(y))/sd(y)

X <- Xmat[, -1]
X <- apply(X, 2, function(x) (x-mean(x))/sd(x))

X.full <- embed(X, p) #care here: X is already first lag of the covariates, hence embed yields p+1 lags
y.full <- y[p:nrow(y), , drop=F]

X.in <- X.full[1:(nrow(X.full)-length.hold.out+h.slct-1),]
X.out <- X.full[(nrow(X.full)-length.hold.out+h.slct):nrow(X.full),,drop=FALSE]

y.in <- y.full[1:(length(y.full)-length.hold.out+h.slct-1)]
y.out <- y.full[(length(y.full)-length.hold.out+h.slct):length(y.full)]

if (data.type == "small"){
  sl.vars <- which(colnames(X) %in% info.sets$M)#which(mc.cracken$SMALL == "x")
  PCA <- FALSE
}else if (data.type == "full"){
  sl.vars <- which(colnames(X) %in% info.sets$XL)
  PCA <- TRUE
}else if (data.type == "AR"){
  sl.vars <- which(colnames(X) %in% "CPIAUCSL")
  PCA <- FALSE
}

Xho <- X.out[h ,sl.vars, drop=F]
X <- X.in[,sl.vars, drop=F]
y <- y.in
yho <- c(0,as.numeric(y.out[h]))

if(cons)
{
  X <- cbind(X,1)
  Xho <- c(Xho, 1)
}

N <- nrow(y)
K <- ncol(X)
v <- N-K  


if (!(G.slct == "SV") && !(G.slct == "mixSV")){
  G.slct <- as.numeric(G.slct)
}

  if (G.slct == "mixSV"){
    G.slct <- 25
    mix.sv <- TRUE
    source("gp_svdlp.R")
  }else{
    mix.sv <- FALSE
  }
  
  if (mod.slct == "reg"){
    mod.j <- gp_bnp(y, X, yho, Xho,nburn=nburn, ntot=ntot, thinning=2, PDP=FALSE, sample.omega=FALSE, tau2.start=0.01, G.max=G.slct, fcst=TRUE, a0=1, a1=1, PCA=PCA, mix.sv=mix.sv)
  }else if (mod.slct == "gp"){
    mod.j <- gp_bnp(y, X, yho, Xho,nburn=nburn, ntot=ntot, thinning=2, PDP=FALSE, sample.omega=FALSE, tau2.start=1e+3, G.max=G.slct, fcst=TRUE, a0=1, a1=1, PCA=PCA, mix.sv=mix.sv)
  }else if (mod.slct =="subspace"){
    mod.j <- gp_bnp(y, X, yho, Xho,nburn=nburn, ntot=ntot, thinning=2, PDP=FALSE, sample.omega=TRUE, tau2.start=0.5, G.max=G.slct, fcst=TRUE, a0=1, a1=1, PCA=PCA, mix.sv=mix.sv)
  }else if (mod.slct =="UCSV"){
    if (data.type=="small"){
     mod.j <- ucsv(y, yho,h, G.max = G.slct, mix.SV = mix.sv)
    }else{
     break()
    }
  }else if (mod.slct == "BART"){
    if (G.slct == "mixSV"){
      mod.j <- BART_mixSV(y, X, yho, Xho,nburn=nburn, ntot=ntot, thinning=2, PDP=FALSE, sample.omega=FALSE, tau2.start=1e+3, G.max=G.slct, fcst=TRUE, a0=1, a1=1, PCA=PCA, mix.sv=mix.sv)
    }else{
      mod.j <- BART_bnp(y, X, yho, Xho,nburn=nburn, ntot=ntot, thinning=2, PDP=FALSE, sample.omega=FALSE, tau2.start=1e+3, G.max=G.slct, fcst=TRUE, a0=1, a1=1, PCA=PCA, mix.sv=mix.sv)
    }
  }
  save(mod.j, file=filename)
  }, error = function(e) {
    writeLines(as.character(e), outputFile)
  })
}


  

