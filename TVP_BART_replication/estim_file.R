#### --------------------------------------------------------------------- ####
#### --------                Main estimation script              --------- ####
#### -------- Bayesian Modeling of Time-Varying Parameter Vector --------- ####
#### --------       Autoregressions Using Regression Trees       --------- ####
#### --------         Hauzenberger, Huber, Koop & Mitchell       --------- ####
#### --------            Annals of Applied Statistics            --------- ####
#### --------------------------------------------------------------------- ####
rm(list = ls())
w.dir <- ""

#### --------------------------------------------------------------------- ####
#### ---------------------------- Libraries ------------------------------ ####
#### --------------------------------------------------------------------- ####

# MCMC / estimation
library(SoftBart)
library(MASS)
library(RcppArmadillo)
library(Rcpp)
library(coda)
library(stochvol)
library(zoo)
library(invgamma)
library(sparseMVN)
library(GIGrvg)
library(scales)
library(shrinkTVP)

# Helper 
library(Matrix)
library(mvtnorm)
library(BVAR)

# Plots
library(withr)
library(ggplot2)
library(crayon)
library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)

source(paste0(w.dir, "aux_funcs/aux_files.R"))
source(paste0(w.dir, "aux_funcs/tvp_softbart_function.R"))

#### --------------------------------------------------------------------- ####
#### ----------------------- Model Specification ------------------------- ####
#### --------------------------------------------------------------------- ####

# MCMC settings
nsave    <- 5000
nsave.2  <- 5000
nburn    <- 5000
nthin    <- 2

# Data / model identifiers
data.set   <- "USMACRO"
info.set   <- "woEBP"
prices     <- c("PCECTPI", "PCEPILFE", "GDPCTPI", "CPILFESL", "CPIAUCSL")
prices.labs <- c("PCECTPI", "PCEPI core", "GDPCTPI", "CPI core", "CPI")
dFreq      <- "Q"
data.freq <- 4
trans      <- "mix_yoy"
p        <- 5
ts.start <- 1975
end.in   <- 2019+3/4
irf.on   <- TRUE
nhor     <- data.freq*4 + 1
cons      <- TRUE
stdz      <- FALSE

shock.ident <- "internal"

if(info.set == "woEBP"){
  shk.var     <- "UNRATE" # Cholesky
  shk.max.var <- "UNRATE" # Maximum variance explained identification
}else if(info.set == "wEBP"){
  shk.var     <- "EBP"    # Cholesky
  shk.max.var <- "UNRATE" # Maximum variance explained identification
}
shk.sc <- "prd-wise"

# Mean specification
mean.type  <- c(25, 1, TRUE)   # Q, num.trees, Q.penalty
Q          <- as.numeric(mean.type[[1]])
num.trees  <- as.numeric(mean.type[[2]])
Q.penalty  <- as.logical(mean.type[[3]])

# Variance specification
var.type     <- c("FHB", "SV", 3, 250)
FSV          <- var.type[[1]]
idioSV       <- var.type[[2]]
R            <- as.numeric(var.type[[3]])
num.trees.sv <- as.numeric(var.type[[4]])

ts.st <- ts.start - p/data.freq

#### --------------------------------------------------------------------- ####
#### -------------------------- Data set-up ------------------------------ ####
#### --------------------------------------------------------------------- ####

load(paste0(w.dir, "data/makroUS_Q_16.rda"))
var.set <- info.sets[["S"]]

Xraw.mix[,"S.P.500"] <- Xraw.stat[,"S.P.500"]
Xraw      <- Xraw.mix[,var.set]
Xraw      <- window(Xraw,      start = ts.st, end = end.in)
var.names <- colnames(Xraw)

# OPTIONAL: Add EBP as endogenous variable or construct external instrument residual
load(paste0(w.dir, "data/EBPUS_Q.rda"))
ebp <- window(ebp[,"ebp"], start = ts.st, end = end.in)
if(info.set == "wEBP"){
  Xraw      <- ts.intersect("EBP" = ebp, Xraw)
  colnames(Xraw) <- c("EBP", var.names)
  var.names <- colnames(Xraw)
}else if(info.set == "woEBP"){
  ebp.instr <- ar.ols(ebp, p = p)
  ebp.instr <- ebp.instr$resid
}

# Standardise or demean
if(stdz){
  sc.mean <- apply(Xraw, 2, mean)
  sc.sd   <- apply(Xraw, 2, sd)
  Xraw    <- apply(Xraw, 2, function(x) (x-mean(x))/sd(x))
  Xraw    <- ts(Xraw, start = ts.st, frequency = data.freq)
}else{
  sc.mean <- apply(Xraw, 2, mean)
  sc.sd   <- rep(1, ncol(Xraw))
  Xraw    <- apply(Xraw, 2, function(x) (x-mean(x)))
  names(sc.mean) <- names(sc.sd) <- colnames(Xraw)
  Xraw    <- ts(Xraw, start = ts.st, frequency = data.freq)
}

# Build regressor matrix X and BART covariate matrix Z
Xraw <- window(Xraw, end = end.in)
TT   <- nrow(Xraw)
X    <- mlag(Xraw, p)

if(cons){
  X <- cbind(X, 1)
  colnames(X) <- c(paste0(rep(var.names, p), "_l", rep(1:p, each = length(var.names))), "iota")
}else{
  colnames(X) <- paste0(rep(var.names, p), "_l", rep(1:p, each = length(var.names)))
}

Z <- X[,paste0(var.names, "_l1")]
colnames(Z) <- paste0(var.names, "_l1")

y <- Xraw[(p+1):TT,,drop = F]
X <- X[(p+1):TT,,drop = F]
Z <- Z[(p+1):TT,,drop = F]

K <- ncol(X)
M <- ncol(y)
N <- nrow(y)

# Normalise Z to [0,1] for BART splitting rules
Z.max <- t(matrix(apply(Z, 2, max), M, N))
Z.min <- t(matrix(apply(Z, 2, min), M, N))
Z.iqr <- Z.max - Z.min

Z      <- (Z - Z.min)/Z.iqr
Z.fsv  <- Z
Z.idio <- matrix((1:N - 1)/(N-1), N, 1)

#### --------------------------------------------------------------------- ####
#### ---------------------- Model and prior lists ------------------------ ####
#### --------------------------------------------------------------------- ####

model.setup <- list()

model.setup$Q         <- Q
model.setup$R         <- R
model.setup$FSV       <- FSV
model.setup$idioSV    <- idioSV

bart.setup <- list(num.trees = num.trees, num.trees.sv = num.trees.sv, cgm.exp = 2, cgm.level = 0.95, Q.penalty = Q.penalty,
                   sd.mu = 2, sd.mu.fsv = 10, sd.mu.sv.idio = 50,
                   fsv.std = TRUE, tau.rate = 100000)
off.setup  <- list(cut.off = 10, off.set.SV = 8, cut.off.V = 1e-4, cut.off.idio = 10, off.set.hs = 1e-15)

irf.var   <- var.names
irf.setup <- list(irf.on = irf.on, nhor = nhor, var.names = var.names, irf.var = irf.var, n.irf = length(irf.var), shk.var = shk.var, shk.max.var = shk.max.var, shk.sc = shk.sc, shock.ident = shock.ident, prices = prices, price.labs = prices.labs, sc.mean = sc.mean, sc.sd = sc.sd, Z.min = Z.min, Z.max = Z.max, Z.iqr = Z.iqr)

#### --------------------------------------------------------------------- ####
#### ------------------------ Output directory --------------------------- ####
#### --------------------------------------------------------------------- ####
dir.create(paste0(w.dir, data.set), showWarnings = FALSE)
dir <- gsub(paste0(w.dir, data.set, "/", dFreq, "_", trans), pattern = " ", replacement = "_")
dir.create(dir, showWarnings = FALSE)
file.name <- gsub(paste0(c(info.set, shock.ident, mean.type, var.type), collapse = "_"), pattern = " ", replacement = "_")
foldername <- paste0(dir, "/", file.name)
dir.create(foldername, showWarnings = FALSE)

#### --------------------------------------------------------------------- ####
#### -------------------- Step 1: Estimation ----------------------------- ####
#### --------------------------------------------------------------------- ####
Regr.obj <- tvp.softbart.func(y = y, X = X, Z = Z, Z.fsv = Z.fsv, Z.idio = Z.idio, p = p, cons = cons, thin = nthin, nsave = nsave, nburn = nburn, model.setup = model.setup, bart.setup = bart.setup, off.setup = off.setup, irf.setup = irf.setup)

list2env(x = Regr.obj, envir = .GlobalEnv)

#### --------------------------------------------------------------------- ####
#### ------------------- Step 2: Post-processing ------------------------- ####
#### --------------------------------------------------------------------- ####
girf_store             <- girf_store[order(apply(eig_store, 1, max)),,,,, drop = F]
PC_store               <- PC_store[order(apply(eig_store, 1, max)),,,,, drop = F]
bt.store               <- bt.store[order(apply(eig_store, 1, max)),,,,drop = F]
sig2t.store            <- sig2t.store[order(apply(eig_store, 1, max)),,,,drop = F]
var.share.sv.store     <- var.share.sv.store[order(apply(eig_store, 1, max)),,,,drop = F]
commonalities.sv.store <- commonalities.sv.store[order(apply(eig_store, 1, max)),,,drop = F]
eht.store              <- eht.store[order(apply(eig_store, 1, max)),,,drop = F]
trees.fit.store        <- trees.fit.store[order(apply(eig_store, 1, max)),,,,drop = F]
Lambda.store           <- Lambda.store[order(apply(eig_store, 1, max)),,,,drop = F]
Lambda.sv.store        <- Lambda.sv.store[order(apply(eig_store, 1, max)),,,drop = F]
var.share.store        <- var.share.store[order(apply(eig_store, 1, max)),,,,drop = F]
commonalities.store    <- commonalities.store[order(apply(eig_store, 1, max)),,,drop = F]
trees.count.store      <- trees.count.store[order(apply(eig_store, 1, max)),,,,drop = F]
trees.count.fsv.store  <- trees.count.fsv.store[order(apply(eig_store, 1, max)),,,drop = F]
rec.share.sv.store     <- rec.share.sv.store[order(apply(eig_store, 1, max)),,, drop = F]
sig2.idio.store        <- sig2.idio.store[order(apply(eig_store, 1, max)),,,drop = F]
eig_store              <- eig_store[order(apply(eig_store, 1, max)),,drop = F]

girf_store             <- girf_store[1:nsave.2,,,,,drop = F]
PC_store               <- PC_store[1:nsave.2,,,,,drop = F]
bt.store               <- bt.store[1:nsave.2,,,,drop = F]
sig2t.store            <- sig2t.store[1:nsave.2,,,,drop = F]
var.share.sv.store     <- var.share.sv.store[1:nsave.2,,,,drop = F]
commonalities.sv.store <- commonalities.sv.store[1:nsave.2,,,drop = F]
rec.share.sv.store     <- rec.share.sv.store[1:nsave.2,,, drop = F]
eht.store              <- eht.store[1:nsave.2,,,drop = F]
trees.fit.store        <- trees.fit.store[1:nsave.2,,,,drop = F]
Lambda.store           <- Lambda.store[1:nsave.2,,,,drop = F]
Lambda.sv.store        <- Lambda.sv.store[1:nsave.2,,,drop = F]
var.share.store        <- var.share.store[1:nsave.2,,,,drop = F]
commonalities.store    <- commonalities.store[1:nsave.2,,,drop = F]
trees.count.store      <- trees.count.store[1:nsave.2,,,,drop = F]
trees.count.fsv.store  <- trees.count.fsv.store[1:nsave.2,,,drop = F]
sig2.idio.store        <- sig2.idio.store[1:nsave.2,,,drop = F]
eig_store              <- eig_store[1:nsave.2,,drop = F]

#### --------------------------------------------------------------------- ####
#### ----------------------- Step 3: Analysis ---------------------------- ####
#### --------------------------------------------------------------------- ####

# --- Time axis, recession dates, sub-samples, and variable labels ---
load(paste0(w.dir, "data/uncUS_Q.rda"))
rec.time   <- exo.df[,"REC"]
rec.time   <- time(rec.time)[rec.time == 1]
rec.time   <- zoo::as.yearqtr(rec.time[c(6,10,11,13,14,19,20,22,23,25,26,31)])
rec_dates  <- data.frame(Peak = 1:6, Trough = 1:6)
rec_dates$Peak   <- zoo::as.Date.yearqtr(rec.time)[seq(1,12,2)]
rec_dates$Trough <- zoo::as.Date.yearqtr(rec.time)[seq(2,12,2)]

time.line <- as.character(seq(as.Date("1975-01-01"), as.Date("2019-12-31"), by = "quarter"))
sub.smp1  <- as.character(seq(as.Date("1975-01-01"), as.Date("1989-12-31"), by = "quarter"))
sub.smp2  <- as.character(seq(as.Date("1990-01-01"), as.Date("2019-12-31"), by = "quarter"))

if(info.set == "wEBP"){
  clean.labs <- data.frame("var"   = c("EBP", "RGDP", "EMPL" , "UNRATE", "AWH", "PCECTPI", "PCEPI core", "GDPCTPI", "CPI"  , "CPI core", "AHE", "FFR"  , "GS10"),
                           "group" = c("SMALL", "v1" , "v1"   , "SMALL" , "v1" , "PRICES" , "PRICES"    , "PRICES" , "SMALL", "PRICES"  , "v1" , "SMALL", "v1"))
}else if(info.set == "woEBP"){
  clean.labs <- data.frame("var"   = c("RGDP", "EMPL" , "UNRATE", "AWH", "PCECTPI", "PCEPI core", "GDPCTPI", "CPI"  , "CPI core", "AHE", "FFR"  , "GS10"),
                           "group" = c("v1"  , "v1"   , "SMALL" , "v1" , "PRICES" , "PRICES"    , "PRICES" , "SMALL", "PRICES"  , "v1" , "SMALL", "v1"))
}

clean.labs.X <- colnames(X)

# --- Posterior summaries: Conditional mean

# Communalities
comm_post  <- apply(commonalities.store, c(2,3), mean)
var.share_post <- apply(var.share.store, c(2,3,4), mean)
rownames(comm_post) <- clean.labs[,1]
colnames(comm_post) <- clean.labs.X
dimnames(var.share_post) <- list("y" = clean.labs[,1], "X" = clean.labs.X, "factor" = paste0("Fac", 1:Q))

comm.sv_post     <- apply(commonalities.sv.store, c(2,3), mean)
var.share.sv_post <- apply(var.share.sv.store, c(2,3,4), mean)
dimnames(comm.sv_post)     <- list("date" = as.character(time.line), "y" = clean.labs[,1])
dimnames(var.share.sv_post) <- list("date" = as.character(time.line), "y" = clean.labs[,1], "factor" = paste0("Fac", 1:R))

si_post <- apply(trees.fit.store, c(2,3,4), function(x) c(quantile(x, 0.16), mean(x), quantile(x, 0.84)))
dimnames(si_post) <- list("moment" = c("low","med","high"), "date" = as.character(time.line), "variable" = clean.labs[,1], "factor" = 1:Q)

mu_post <- apply(Lambda.store, c(2,3,4), function(x) c(quantile(x, 0.16), mean(x), quantile(x, 0.84)))
dimnames(mu_post) <- list("moment" = c("low","med","high"), "y" = clean.labs[,1], "X" = clean.labs.X, "factor" = 1:Q)

# --- Posterior summaries: Variances

f.var_post <- apply(f.sv.store, c(2,3), function(x) c(quantile(x, 0.16), mean(x), quantile(x, 0.84)))
dimnames(f.var_post) <- list("moment" = c("low","med","high"), "date" = as.character(time.line), "factor" = 1:R)

si.var_post <- apply(sqrt(eht.store), c(2,3), function(x) c(quantile(x, 0.16), mean(x), quantile(x, 0.84)))
dimnames(si.var_post) <- list("moment" = c("low","med","high"), "date" = as.character(time.line), "factor" = 1:R)

mu.var_post <- apply(Lambda.sv.store, c(2,3), function(x) c(quantile(x, 0.16), mean(x), quantile(x, 0.84)))
dimnames(mu.var_post) <- list("moment" = c("low","med","high"), "variable" = clean.labs[,1], "factor" = 1:R)

trees.count_post <- apply(trees.count.store, c(2,3,4), function(x) c(quantile(x, 0.16), mean(x), quantile(x, 0.84)))
dimnames(trees.count_post) <- list("moment" = c("low","med","high"), "y" = clean.labs[,1], "factor" = 1:Q, "Z" = colnames(Z))

trees.count.fsv_post <- apply(trees.count.fsv.store, c(2,3), function(x) c(quantile(x, 0.16), mean(x), quantile(x, 0.84)))
dimnames(trees.count.fsv_post) <- list("moment" = c("low","med","high"), "factor" = 1:R, "Z" = colnames(Z))

# --- IRFs: Stability check

# Set draws with max eigenvalue > eigen.crit to NA
eigen.crit <- 1.05

irf_stab       <- girf_store
PC_stab        <- PC_store
sig2.t_stab    <- sig2t.store
bt_stab        <- bt.store
sig2.idio_stab <- sig2.idio.store
for(rr in 1:N){
  irf_stab[(eig_store[,rr]>eigen.crit),rr,,,]   <- NA
  PC_stab[(eig_store[,rr]>eigen.crit),rr,,,]    <- NA
  sig2.t_stab[(eig_store[,rr]>eigen.crit),rr,,] <- NA
  bt_stab[(eig_store[,rr]>eigen.crit),rr,,]     <- NA
  sig2.idio_stab[(eig_store[,rr]>eigen.crit),rr,] <- NA
}

dimnames(irf_stab) <- list("draws" = 1:nsave.2, "date" = as.character(time.line), "variable" = clean.labs[,1], "horizon" = 0:(nhor-1), "shock" = c(paste0("fac", 1:R), "chol"))
dimnames(PC_stab)  <- list("draws" = 1:nsave.2, "date" = as.character(time.line), "variable" = prices.labs,    "horizon" = 0:(nhor-1), "shock" = c(paste0("fac", 1:R), "chol"))


# --- Posterior summaries: IRFs

# Full-sample time-varying
irf_post <- apply(irf_stab, c(2,3,4,5), quantile, probs = c(.16,.5,.84), na.rm = T)
dimnames(irf_post) <- list("moment" = c("low","med","high"), "date" = as.character(time.line), "variable" = clean.labs[,1], "horizon" = 0:(nhor-1), "shock" = c(paste0("fac", 1:R), "chol"))

PC_post <- apply(PC_stab, c(2,3,4,5), quantile, probs = c(.16,.5,.84), na.rm = T)
dimnames(PC_post) <- list("moment" = c("low","med","high"), "date" = as.character(time.line), "variable" = prices.labs, "horizon" = 0:(nhor-1), "shock" = c(paste0("fac", 1:R), "chol"))

# Average over full sample
irf_avg <- apply(irf_stab, c(1,3,4,5), mean, na.rm = T)
PC_avg  <- apply(PC_stab,  c(1,3,4,5), mean, na.rm = T)

irf_avg_post <- apply(irf_avg, c(2,3,4), quantile, probs = c(.16,.5,.84), na.rm = T)
dimnames(irf_avg_post) <- list("moment" = c("low","med","high"), "variable" = clean.labs[,1], "horizon" = 0:(nhor-1), "shock" = c(paste0("fac", 1:R), "chol"))

PC_avg_post <- apply(PC_avg, c(2,3,4), quantile, probs = c(.16,.5,.84), na.rm = T)
dimnames(PC_avg_post) <- list("moment" = c("low","med","high"), "variable" = prices.labs, "horizon" = 0:(nhor-1), "shock" = c(paste0("fac", 1:R), "chol"))

# Volatility and coefficients
sig2.t_post <- apply(sig2.t_stab, c(2,3,4), quantile, probs = c(.16,.5,.84), na.rm = T)
dimnames(sig2.t_post) <- list("moment" = c("low","med","high"), "date" = as.character(time.line), "var1" = clean.labs[,1], "var2" = clean.labs[,1])

sig2.idio_post <- apply(sig2.idio_stab, c(2,3), quantile, probs = c(.16,.5,.84), na.rm = T)
dimnames(sig2.idio_post) <- list("moment" = c("low","med","high"), "date" = as.character(time.line), "var1" = clean.labs[,1])

bt_post <- apply(bt_stab, c(2,3,4), quantile, probs = c(.16,.5,.84), na.rm = T)
dimnames(bt_post) <- list("moment" = c("low","med","high"), "date" = as.character(time.line), "var1" = clean.labs[,1], "var2" = clean.labs.X)


save(file = paste0(foldername, "_", shk.sc, "_draws.Rda"),
     list = c("irf_post", "PC_post",
              "irf_avg_post", "PC_avg_post",
              "comm_post", "comm.sv_post",
              "var.share_post", "var.share.sv_post",
              "sig2.t_post", "sig2.idio_post", "bt_post",
              "si_post", "mu_post", "si.var_post",
              "f.var_post", "mu.var_post",
              "trees.count_post", "trees.count.fsv_post",
              "Z", "M", "p", "cons", "y", "K", "N"))
