model.est <- c("UCSV")
tau.grid <- c(1e-4, 1e+4)
G.grid <- c("mixSV", "1","25", "SV")
eval_type <- "RT"                      # "pseudo" versus "RT" (real-time)
final.vint <- 2021+10/12               # final vintage date (November 2021)
roll <- TRUE  
length.hold.out <- 185
h <- c(1, 4, 8)

#data.grid <- "small"
sl.var <- "CPIAUCSL"
sl.var <- "CPILFESL"

sl.quants <- seq(0.05, 0.95, by=0.05)
loss <- c("lpl", "se", paste0("qs", sl.quants))
data.grid <- c("small")
pit.array <- store.array <- array(NA, c((length.hold.out), length(model.est), length(G.grid), length(h),length(loss), length(data.grid))); dimnames(store.array)<- dimnames(pit.array) <- list(1:length.hold.out, model.est, G.grid, h, loss,data.grid)
pred.array <- array(NA, c((length.hold.out), length(model.est), length(G.grid), length(h),length(sl.quants)+1, length(data.grid))); dimnames(pred.array)  <- list(1:length.hold.out, model.est, G.grid, h, c(sl.quants,"outcome"), data.grid)
#pit.array <- array(NA, c(length.hold.out, length(model.est), length(G.grid), length(h),length(sl.quants)+1, length(data.grid))); dimnames(pit.array) <- list(1:length.hold.out, model.est, G.grid, h, c(sl.quants,"outcome"), data.grid)

for (mod.type in data.grid){
    for (month in 1:length.hold.out){
          for (h.sl in h){
      for (mod.slct in model.est){
        for (G.slct in G.grid){
          if (mod.type == "small"){
            outputFile <- paste0("Results_", sl.var,"/", "res_time_",month, "_fhorz_", h.sl, "_model_", mod.slct, "_G_", G.slct,"_data_",mod.type,".RData")
          }else{
            outputFile <- paste0("Results_", sl.var,"/", "res_time_",month, "_fhorz_", h.sl, "_model_", mod.slct, "_G_", G.slct,".RData")
          }
          if (file.exists(outputFile)){
            load(outputFile)
            #This part computes the various loss functions to be defined above
            for (j in loss){
              if (j == "lpl"){
                y.median <- median(mod.j$predictions, na.rm=T)
                sd.median <- sd(mod.j$predictions, na.rm=T)
                store.array[(month), mod.slct, as.character(G.slct), as.character(h.sl), j,mod.type] <- dnorm(mod.j$yho, y.median, sd.median, log=TRUE)#log(mean(exp(mod.j$lpl-7), na.rm=T))+7
                pit.array[(month), mod.slct, as.character(G.slct), as.character(h.sl), j, mod.type] <- as.numeric(qnorm(pnorm(mod.j$yho, y.median, sd.median)))
                #store.array[(month), mod.slct, as.character(G.slct), as.character(h.sl), j,mod.type] <- dnorm(mod.j$yho, median(mod.j$predictions, na.rm=T), sd(mod.j$predictions, na.rm=T), log=TRUE)#log(mean(exp(mod.j$lpl-7), na.rm=T))+7
              }else if (j == "se"){
                store.array[(month), mod.slct, as.character(G.slct), as.character(h.sl), j, mod.type] <- (median(mod.j$predictions, na.rm=T) - mod.j$yho)
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

res.uc.list <- list("forecasts"=pred.array, "res"=store.array, "pit"=pit.array)
save(res.uc.list, file=paste0("ucresults_", sl.var, ".RData"))
