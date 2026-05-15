require(stochvol); require(ggplot2); library(reshape2)
set.seed(123) #set seed 
source("Data_QD/prepare_dat.R") # Careful here: this script always downloads the latest vintage of the McCracken Ng dataset
source("aux.R")

var.S.true <- c("GDPC1", "UNRATE", "CPIAUCSL","FEDFUNDS", "GS10TB3Mx","BAA10YM","S.P.500") # Select a few variables
Xraw <- Xraw[,var.S.true]

X.mean <- apply(Xraw, 2, mean) # Gets the mean of the series
X.sd <- apply(Xraw, 2, sd) # Gets the standard deviations of the series
Xraw <- apply(Xraw, 2, function(x) (x-mean(x))/sd(x))  #normalizes the series 

M <- ncol(Xraw) # Total number of series used in our VAR

p <- 5 # Number of lags
nsave <- 5000 # Number of saved draws
nburn <- 5000 # Number of burn-ins
G.max <- 20 # Number of total clusters (will typically be smaller when we estimate the model)
SV.idi <- FALSE # In case we have SV on the idiosyncratic shocks

# Start by estimating the model with BNP shocks first
bvar.G <- bvar.mix(Xraw, nsave=nsave, nburn=nburn, G.max=G.max, p=p, SV.idi = SV.idi, fhorz = 1, irf.hor = 24,prior.cov = "AR", sample.w = "conditional")

# Estimate the model with G=1 
bvar.G1 <- bvar.mix(Xraw, nsave=nsave, nburn=nburn, G.max=1, p=p, SV.idi = SV.idi, fhorz = 1, irf.hor = 24, prior.cov = "AR", sample.w = "conditional")

# Compute posterior moments of the IRFs
IRF.post.G <- apply(bvar.G$IRF, c(2,3,4,5), function(x) quantile(x,c(0.16, 0.5, 0.84),na.rm=TRUE))
IRF.post.G1 <- apply(bvar.G1$IRF, c(2,3,4,5), function(x) quantile(x,c(0.16, 0.5, 0.84),na.rm=TRUE))
dimnames(IRF.post.G1) <- dimnames(IRF.post.G) <- list(c("0.16", "0.5", "0.84"), colnames(Xraw), colnames(Xraw),NULL, NULL)

foldername <- "Insample_results" #Creates a folder that stores all the results
dir.create(foldername, showWarnings = FALSE)

#Create a simple table that shows the number of non-zero components
G.prob <- matrix(NA, G.max, 1)
for (jj in 1:G.max){
  G.prob[jj,] <- sum(bvar.G$G[,2]==jj)/nsave
}

G.tab <- t(G.prob[1:9])

#Start by plotting the log(determinant of Sigma)
logdet.store <- array(NA, c(nsave, G.max))
for (j in seq_len(nsave)){
  for (g in seq_len(G.max)){
    logdet.store[j,g] <- determinant(bvar.G$Sigma_G[j,,,g], log=TRUE)$modulus
  }
}
colnames(logdet.store) <- paste0("Cluster ",seq(1,G.max))
logdet.vec <- melt(logdet.store[, 1:5])
det.prior <- median(logdet.store[,8])# Selects 

pdf(paste0(foldername,"/detSIGMA.pdf"), width=14, height=3)
c <- ggplot(logdet.vec, aes(x=Var2, y=value)) + 
  geom_boxplot(outlier.shape = NULL, outlier.size = 0) + ylab("") + xlab("") + theme_bw() + geom_hline(yintercept = det.prior, col = "red")
print(c)
dev.off()

#Time series plots of the regime allocation
n.length <- 4
S.median <- apply(bvar.G$S,2,median)
S.prob <- matrix(0, length(S.median), G.max)
for (j in seq_len(G.max)){
  S.prob[,j] <- apply(bvar.G$S==j,2,mean)
}

S.prob <- S.prob[, sort(apply(S.prob,2,mean),decreasing = TRUE, index.return=T)$ix]
S.prob <- apply(S.prob, 2, function(x) rollmean(x, n.length))

date <- seq(1, nrow(S.prob))
time.line <- as.character(seq(as.Date("1961-01-01"),as.Date("2022-12-30"),by="quarter"))
#time.line <- time.line[(p+1):length(time.line)]
#time.line <- time.line[(n.length-3) : length(time.line)]
time.line <- time.line[2:length(time.line)]
rownames(S.prob) <- time.line
colnames(S.prob) <- paste0("probability of being in cluster ",seq(1, G.max))

ind.mat <- t(matrix(seq(1:10), 2, 5))

jj <- 1

list.prob <- list()

for (jj in seq_len(5)){

S.mat <- melt(S.prob[, ind.mat[jj,]])
S.mat$Var1 <- as.Date(S.mat$Var1)
S.mat$Var2 <- as.factor(S.mat$Var2)
prob.ot <- ggplot(S.mat, aes(x=Var1, y=value, col=Var2))+ 
  geom_line()+ theme_bw()+ ylim(c(0,1)) + ylab("") + xlab("")+
  theme(legend.title=element_blank(), legend.position = "bottom")# + nberShade()

pdf(paste0(foldername,"/prob_regime_", jj, ".pdf"), width=14, height=3)
print(prob.ot)
dev.off()
}

eta.mean <- apply(bvar.G$eta,2,mean)

#Plot impulse responses for all variables and across different values of G
library(gridExtra)
sl.var <- "GDPC1"
for (sl.var in colnames(Xraw)){

  y.range <- range(X.sd[sl.var]*IRF.post.G[,sl.var,"FEDFUNDS",,],X.sd[sl.var]*IRF.post.G1[,sl.var,"FEDFUNDS",,], na.rm=T)
  
  plot.G <- list()
  for (sl.G in seq_len(10)){
    irf.mat <- X.sd[sl.var]*cbind(t(IRF.post.G[,sl.var, "FEDFUNDS",,sl.G]), t(IRF.post.G1[,sl.var,"FEDFUNDS",,1]))
    plot.var <- as.data.frame(cbind(1:24,irf.mat))
    colnames(plot.var) <- c("Quarters","BVARBNP_016", "BVARBNP_050", "BVARBNP_084", "BVAR_016", "BVAR_050", "BVAR_084")
    
    
    h <- ggplot(plot.var, aes(x=Quarters, y=BVARBNP_050))
    h <- h + ylab("") + ylim(y.range)
    # Add aesthetic mappings 
    h <- h + geom_ribbon(
      aes(ymin = BVARBNP_016, ymax = BVARBNP_084), fill = "grey70") + 
      geom_line(aes(y = BVARBNP_050))+
      geom_line(aes(y = BVAR_016), col="orange")+
      geom_line(aes(y = BVAR_084), col="orange")+
      ggtitle(paste0("Cluster ", sl.G))+
      geom_hline(yintercept = 0, col="red", lty="dashed")+
      theme_bw()
    
    plot.G[[sl.G]] <- h
  }
  
  pdf(paste0(foldername,"/", sl.var, ".pdf"), width=14, height=3)
  grid.arrange(plot.G[[1]], plot.G[[2]],plot.G[[3]],plot.G[[4]], plot.G[[5]], ncol=5)
  dev.off()
  
  
  #Now plot only posterior medians across clusters
  plot.var <- IRF.post.G["0.5",sl.var,"FEDFUNDS",,]*X.sd[sl.var]
  pdf(paste0(foldername,"/", sl.var, "medians.pdf"), width=4, height=4)
  par(mar=c(2,2,2,2))
  ts.plot(plot.var, lty=rep(1, 20), xlab="Quarters", ylab=""); abline(h=0, col="red"); grid()
  for (jj in 1:ncol(plot.var)) points(plot.var[,jj], pch=jj)
  if (sl.var == "CPIAUCSL") legend("topright",legend=paste0("J=",1:5), lty=1, pch = 1:5)
  dev.off()
  
  
  #Plot impulse responses weighted by eta
  IRF.avg <- array(NA, c(3, M, 24)); dimnames(IRF.avg) <- list(NULL, colnames(Xraw), NULL)
  IRF.avg.high <- IRF.avg.median <- IRF.avg.low <- matrix(0, M, 24)
  for (g in 1:5){
    IRF.avg.low <- IRF.avg.low + eta.mean[g]*IRF.post.G[1,,"FEDFUNDS",,g]
    IRF.avg.median <- IRF.avg.median + eta.mean[g]*IRF.post.G[2,,"FEDFUNDS",,g]
    IRF.avg.high <- IRF.avg.high + eta.mean[g]*IRF.post.G[3,,"FEDFUNDS",,g]
  }
  IRF.avg[1,,] <- IRF.avg.low
  IRF.avg[2,,] <- IRF.avg.median
  IRF.avg[3,,] <- IRF.avg.high
  
  IRF.avg.mat <- IRF.avg[,sl.var,]
  plot.var <- as.data.frame(cbind(1:24,t(IRF.avg.mat), t(IRF.post.G1[,sl.var,"FEDFUNDS",,1])))
  colnames(plot.var) <- c("Quarters","BVARBNP_016", "BVARBNP_050", "BVARBNP_084", "BVAR_016", "BVAR_050", "BVAR_084")
  
  y.range <- range(cbind(IRF.avg.mat,IRF.post.G1[,sl.var,"FEDFUNDS",,1]))
  h <- ggplot(plot.var, aes(x=Quarters, y=BVARBNP_050))
  h <- h + ylab("") + ylim(y.range)
  # Add aesthetic mappings 
  h <- h + geom_ribbon(
    aes(ymin = BVARBNP_016, ymax = BVARBNP_084), fill = "grey70") + 
    geom_line(aes(y = BVARBNP_050))+
    geom_line(aes(y = BVAR_016), col="orange")+
    geom_line(aes(y = BVAR_084), col="orange")+
    geom_hline(yintercept = 0, col="red", lty="dashed")+
    theme_bw()
  
  pdf(paste0(foldername,"/", sl.var, "_averaged.pdf"), width=4, height=4)
  print(h)
  dev.off()
}

#Make IRFs for selected time horizons over time
date.new <- as.character(seq(as.Date("1960-04-01"),as.Date("2023-12-31"),by="quarter"))
date.new <- date.new[(p+1):length(date.new)]
sl.horz <- c(2, 4, 8)
sl.quants <- c(0.16, 0.5, 0.84)

IRF.post <- bvar.G$IRF; dimnames(IRF.post) <- list(NULL, colnames(Xraw), colnames(Xraw), NULL, NULL)
IRF.mat <- array(NA, c(nsave, length(date.new), M, length(sl.horz))); dimnames(IRF.mat) <- list(NULL,date.new, colnames(Xraw), sl.horz)
for (t in seq_len(length(date.new))){
  count <- 0
  for (fhorz in sl.horz){
    count <- count+1
    for (irep in 1:nsave){
      delta.t <- bvar.G$S[irep, t]
      irf.t <- IRF.post[irep,,"FEDFUNDS", fhorz, delta.t] * X.sd
      
      IRF.mat[irep, t, , count] <- irf.t
    }
  }
  print(t)
}

# sl.dates <- c("1981-01-01","2008-10-01")
IRF.quants <- apply(IRF.mat, c(2,3,4), function(x) quantile(x, c(0.16, 0.84)))
# 
# 
# 
for (var.slct in colnames(Xraw)){
  plot.var <- as.data.frame(cbind(1:length(date.new),t(IRF.quants[,,var.slct, "4"]), t(IRF.quants[,,var.slct, "8"])))
  plot.var$V1 <- as.Date(date.new)
  colnames(plot.var) <- c("Quarters","BVARBNP_016", "BVARBNP_084","BVARBNP_016_8",  "BVARBNP_084_8")


  h <- ggplot(plot.var, aes(x=Quarters, y=BVARBNP_016))
  h <- h + ylab("")
  # Add aesthetic mappings
  h <- h + geom_ribbon(
    aes(ymin = BVARBNP_016, ymax = BVARBNP_084), fill = "grey70") +
    #geom_line(aes(y = BVARBNP_050))+
    geom_line(aes(y = BVARBNP_016_8), col="darkblue", lty="dashed")+
    geom_line(aes(y = BVARBNP_084_8), col="darkblue", lty="dashed")+
    geom_hline(yintercept = 0, col="red", lty="dotted")+
    theme_bw()

  pdf(paste0(foldername,"/irfot_", var.slct, ".pdf"), width=14, height=3)
  print(h)
  dev.off()

}



