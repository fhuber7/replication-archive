remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
setwd("!inflation_GP_BNP/Quarterly")
model.est <- c("reg", "gp", "subspace", "UCSV", "BART")
sl.var <- "CPIAUCSL"
#sl.var <- "CPILFESL"

tau.grid <- c(1e-4, 1e+4)
G.grid <- c(1, 25, "SV", "mixSV")
length.hold.out <- 185
h <- c(1, 4, 8)
time.grid <- seq(1, length.hold.out)

#data.grid <- "small"

sl.quants <- seq(0.05, 0.95, by=0.05)
loss <- c("lpl", "se", paste0("qs", sl.quants))
data.grid <- c("AR","small", "full")
pit.array <- store.array <- array(NA, c(length(time.grid), length(model.est), length(G.grid), length(h),length(loss), length(data.grid))); dimnames(pit.array) <- dimnames(store.array) <- list(time.grid, model.est, G.grid, h, loss,data.grid)
pred.array <- array(NA, c(length(time.grid), length(model.est), length(G.grid), length(h),length(sl.quants)+1, length(data.grid))); dimnames(pred.array) <- list(time.grid, model.est, G.grid, h, c(sl.quants,"outcome"), data.grid)



for (mod.type in data.grid){
  for (month in time.grid){
    for (h.sl in h){
      for (mod.slct in model.est){
        for (G.slct in G.grid){
        
          outputFile <- paste0("Results_", sl.var,"/", "res_time_",month, "_fhorz_", h.sl, "_model_", mod.slct, "_G_", G.slct,"_data_",mod.type,".RData")
          
          if (file.exists(outputFile)){
            load(outputFile)
            #This part computes the various loss functions to be defined above
            for (j in loss){
              if (j == "lpl"){
                y.median <- median(mod.j$predictions, na.rm=T)
                sd.median <- sd(mod.j$predictions, na.rm=T)
                
                lpl.i <- remove_outliers(mod.j$lpl)
                
                store.array[(month), mod.slct, as.character(G.slct), as.character(h.sl), j,mod.type] <- log(mean(exp(lpl.i-7), na.rm=T))+7#dnorm(mod.j$yho, y.median, sd.median, log=TRUE)#log(mean(exp(mod.j$lpl-7), na.rm=T))+7
                pit.array[(month), mod.slct, as.character(G.slct), as.character(h.sl), j,mod.type] <- as.numeric(qnorm(pnorm(mod.j$yho, y.median, sd.median)))
              }else if (j == "se"){
                store.array[(month), mod.slct, as.character(G.slct), as.character(h.sl), j, mod.type] <- (median(mod.j$predictions, na.rm=T) - mod.j$yho)^2
              }else if (substr(j,1,2) == "qs"){
                q.pred <- quantile(mod.j$predictions, sl.quants, na.rm=T)
                count <- 0
                y.h <- mod.j$yho
                for (q in sl.quants){
                  count <- count+1  
                  QS.q <- (y.h - q.pred[[count]])*(q - (y.h <= q.pred[[count]])*1)
                  #print(QS.q)
                  store.array[(month), mod.slct, as.character(G.slct), as.character(h.sl), paste0("qs", q), mod.type] <- QS.q
                }
              }
            }
            
            pred.array[(month), mod.slct, as.character(G.slct), as.character(h.sl), ,mod.type] <-  c(quantile(mod.j$predictions, sl.quants, na.rm=T),mod.j$yho)
            
          }else{
            next
          }
        }
      }
    }
  }
}
# 
# pdf("test.pdf")
# ts.plot(pred.array[,3,1,1,], col=c(rep(1,length(sl.quants)),2))
# dev.off()

res.list <- list("forecasts"=pred.array, "res"=store.array, "pit"=pit.array)
save(res.list, file=paste0("results_",sl.var,".RData"))

res.mean <- apply(store.array, c(2,3,4,5,6), mean, na.rm=T)

