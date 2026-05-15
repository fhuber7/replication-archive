model.est <- c("gp", "subspace")
G.grid <- c( 25,  "mixSV")
length.hold.out <- 185
h <- c(1, 4)
time.grid <- seq(1, length.hold.out)


features <- c("fevds", "G", "vola", "mean")
data.grid <- c("small", "full")

features.store <- array(NA, c(length.hold.out, length(h), length(model.est), length(G.grid), length(features), length(data.grid)))
dimnames(features.store) <- list(time.grid, h, model.est, G.grid, features, data.grid)

for (month in time.grid){
  for (h.sl in h){
    for (mod.slct in model.est){
      for (G.slct in G.grid){
        for (mod.type in data.grid){
         outputFile <- paste0("Results_CPIAUCSL/", "res_time_",month, "_fhorz_", h.sl, "_model_", mod.slct, "_G_", G.slct,"_data_",mod.type,".RData")  
         
         if (!file.exists(outputFile)) next
         
         load(outputFile) 
         
         for (sl.feature in features){
           if (sl.feature == "fevds"){
             #Compute FEVDs
             features.store[as.character(month), as.character(h.sl), mod.slct, as.character(G.slct), sl.feature,mod.type] <- median(mod.j$fevds, na.rm=T)
           }else if (sl.feature == "G"){
             features.store[as.character(month), as.character(h.sl), mod.slct, as.character(G.slct), sl.feature,mod.type] <- mod.j$S[length(mod.j$S)]
           }else if (sl.feature == "vola"){
             features.store[as.character(month), as.character(h.sl), mod.slct, as.character(G.slct), sl.feature, mod.type] <- mod.j$sd.mix[length(mod.j$sd.mix)]
           }else if (sl.feature == "mean"){
             features.store[as.character(month), as.character(h.sl), mod.slct, as.character(G.slct), sl.feature, mod.type] <- mod.j$mu.mix[length(mod.j$mu.mix)]
           }
         }
        }
      }
    }
  }
}

save(features.store, file="features.RData")
