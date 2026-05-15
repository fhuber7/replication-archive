rm(list = ls())
require(readxl)
require(zoo)
library(httr)
library(ustyc)
library(dplyr)
library(quantmod)
library(lubridate)

sl.inf <- "C"

data.transform.func <- function(data = data, transform.code = rep(4,ncol(data)), freq = "Q")
{
  if(length(transform.code) != ncol(data)){stop("Transformation code does not equal number of variables!")}
  
  if(freq == "M")
  {
    n.freq <- 12
    scale <- 1200
  }else{
    scale <- 400
    n.freq <- 4
  } 
  
  for(ii in 1:ncol(data))
  {  
    
    series.temp <- as.numeric(data[,ii])
    trans.temp <- as.numeric(transform.code[[ii]])
    if(trans.temp == 4) # log-transformation
    {
      series.temp <- log(series.temp)*100
    }else if(trans.temp == 5){ # log.diff-transformation
      series.temp <- c(NA, diff(log(series.temp)))*100
    }else if(trans.temp == 50){  # yoy-transformation
      series.temp <- c(log(series.temp) - dplyr::lag(log(series.temp),n.freq))*100
    }else if(trans.temp == 6){ # log.diff^2-transformation
      series.temp <- c(NA, NA, diff(diff(log(series.temp))))*100
    }else if(trans.temp == 7){
      series.temp <- c(NA,diff(series.temp))/lag(series.temp) 
    }else if(trans.temp == 2){ #differences in levels
      series.temp <- c(NA, diff(series.temp))
    }
    data[,ii] <- series.temp # levels
  }
  return(data)
}



freq <- "Q"

if(freq == "Q"){
  url_fred <- "https://files.stlouisfed.org/files/htdocs/fred-md/quarterly/current.csv"
  GET(url_fred, write_disk(tf <- tempfile(fileext = "current.csv")))   
  data <- read.csv(tf)
  var.info <- read.csv2(file = "Data_QD/US_QDvariables.csv")
}else if(freq == "M"){
  url_fred <- "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv"
  GET(url_fred, write_disk(tf <- tempfile(fileext = "current.csv")))   
  data <- read.csv(tf)
  var.info <- read.csv2(file = "US_MDvariables.csv")
}

var.info <- var.info[,c("ID", "FRED.Mnemonic", "type", "XS", "S", "M", "L", "XL")]
var.info <- var.info[order(var.info[,"type"], decreasing = TRUE),]

var.XS <- var.info[var.info[,"XS"] == "x","FRED.Mnemonic"]
var.S <- var.info[var.info[,"S"] == "x","FRED.Mnemonic"]
var.M <-var.info[var.info[,"M"] == "x","FRED.Mnemonic"]
var.L <-var.info[var.info[,"L"] == "x","FRED.Mnemonic"]
var.XL <-var.info[var.info[,"XL"] == "x","FRED.Mnemonic"]


var.M <- var.M[!(var.M %in% c("HOUST", "PERMIT", "NONREVSLx", "TOTRESNS", "NONBORRES", "AAAFFM"))]
var.M <- c(var.M, "PPIACO", "WPSID61", "WPSID62", "COMPRNFB", "ULCNFB","BAA10YM")
var.XL <- c(var.XL, "ULCNFB")
#Specific suggestions for dropping:   "HOUST"         "PERMIT"       "NONREVSLx"     "TOTRESNS"      "NONBORRES”   “AAAFFM"
#Adding:  PPIACO, WPSID61, WPSID62, COMPRNFB, ULCNFB

info.sets <- list(XS = var.XS, S = var.S, M = var.M, L = var.L, XL = var.XL)
no.vars   <- list(XS = length(var.XS), S = length(var.S), M = length(var.M), L = length(var.L), XL = length(var.XL))

var.choice <- var.XL

if(freq == "Q"){
  trans.code.sugg <- data[2,-1]
  data <- ts(data[-c(1,2),-1], start = 1959, frequency = 4)
}else if(freq == "M"){
  trans.code.sugg <- data[1,-1]
  data <- ts(data[-1,-1], start = 1959, frequency = 12)
  
}
names(trans.code.sugg) <- colnames(data)

missing.var <- setdiff(var.choice, colnames(data))
na.ind <- apply(data, 2, function(x)sum(is.na(x)))
na.var <- colnames(data)[as.numeric(which(na.ind > 5))]
missing.var <- c(missing.var, na.var)

var.choice <- setdiff(var.choice, missing.var)
var.XL <- setdiff(var.XL, missing.var)
var.L <- setdiff(var.L, missing.var)

data <- data[,var.choice]
trans.code.sugg <- trans.code.sugg[,var.choice]
names(trans.code.sugg) <- colnames(data)

# trans code stationary
trans.code.stat <- trans.code.sugg
# trans code integrated 
levels.var <- names(trans.code.sugg)[which(trans.code.sugg[1,] < 4)]
log.diff.var <- names(trans.code.sugg)[which(trans.code.sugg == 6)]
perc.var <- names(trans.code.sugg)[which(trans.code.sugg == 7)]
trans.code.int <- trans.code.sugg
trans.code.int[1,] <- rep(4, ncol(data)) 
trans.code.int[1,levels.var] <- 1
trans.code.int[1,log.diff.var] <- 5
trans.code.int[1,perc.var] <- 7

Xraw.lev <- data.transform.func(data = data, transform.code = trans.code.int, freq = freq)

#trans.code.stat[, sl.inf] <- 0 
#trans.code.stat <- trans.code.stat[-length(trans.code.stat)]
Xraw.stat <- data.transform.func(data = data, transform.code = trans.code.stat, freq = freq)
Xraw.stat <- window(Xraw.stat, start = c(1960,1))
Xraw <- Xraw.stat[!is.na(apply(Xraw.stat,1,mean)),]

lvl.dat <- list("data"=Xraw.lev, "s"=var.XS, "m"=var.S, "l"=var.M)


