sl.var <- "CPIAUCSL"
model.est <- c("subspace", "gp","UCSV")
G.grid <- c("1","25","mixSV", "SV")
data.grid <- c("small","full")
length.hold.out <- 185
h <- c(1, 4)
time.grid <- seq(1, length.hold.out)

nsave <- 1000
pred.store <- array(NA, c(nsave,length.hold.out, length(h), length(model.est), length(G.grid),  length(data.grid)))

dimnames(pred.store) <- list(NULL,time.grid, h, model.est, G.grid, data.grid)

for (month in time.grid){
  for (h.sl in h){
    for (mod.slct in model.est){
      for (G.slct in G.grid){
        for (mod.type in data.grid){
          outputFile <- paste0("Results_", sl.var,"/", "res_time_",month, "_fhorz_", h.sl, "_model_", mod.slct, "_G_", G.slct,"_data_",mod.type,".RData")  
          if (!file.exists(outputFile)) next
          load(outputFile) 
          sl.samples <- sample(1:nrow(mod.j$predictions), nsave)
        
          pred.store[,as.character(month), as.character(h.sl), mod.slct, as.character(G.slct),mod.type] <- mod.j$predictions[sl.samples,]
        
        }
      }
    }
  }
}

save(pred.store, file="densities.RData")
