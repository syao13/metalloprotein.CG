#!/usr/bin/Rscript

#### Load the data
library(randomForest)

options(stringsAsFactors=FALSE)
#setwd("../output_allMetal")
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

load("angle.correction.RData")
rawdata <- read.table("rf.smallest.txt", header = FALSE)
id <- rawdata[,1]
orderid <- order(id)
data <- rawdata[orderid,]
ligCombos <- names(table(angle.correction$lig_combo))

getGroup <- function (x) {
  if (paste(x[3], x[4], sep = ".") %in% ligCombos) {"compressed"}
  else {"normal"}
}
#################### prepare the data for rf ####################

## Training data
prediction.all <- apply(data, 1, getGroup) 

save(prediction.all, file="rf.results.RData")

