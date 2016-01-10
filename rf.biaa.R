#!/usr/bin/Rscript

#### Load the data
library(randomForest)

options(stringsAsFactors=FALSE)
#setwd("../output_allMetal")
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

rawdata <- read.table("rf.smallest.fe.txt", header = FALSE)
id <- rawdata[,1]
orderid <- order(id)
data <- cbind(rawdata[orderid,2:5], (rawdata[orderid,7] + rawdata[orderid,8])/2)

getGroup <- function (x) {
  if (as.numeric(x[1]) * 0.04 + as.numeric(x[5]) > 4.7)  {"normal"}
  else {"compressed"}
}
#################### prepare the data for rf ####################

## Training data
prediction.all <- apply(data, 1, getGroup) 

save(prediction.all, file="rf.results.RData")

