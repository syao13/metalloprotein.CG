#!/usr/bin/Rscript

############################################################################################
##
##   Written by Sen Yao, 07/20/2016
##   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/20/2016. All rights reserved.
##
###########################################################################################

#### Load the data
library(randomForest)

options(stringsAsFactors=FALSE)
#setwd("../output_allMetal")
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

rawdata <- read.table("rf.smallest.txt", header = FALSE)
data <- rawdata[, 2:5]
table(rawdata[,6])

########## ligand Encoding ########### 
ligs <- c(data[,2], data[,3])
length(table(ligs))
sort(table(ligs),decreasing=TRUE)[1:20]
#sum(sort(table(ligs),decreasing=TRUE)[9:108])

ligRank <- names(sort(table(ligs),decreasing=TRUE))

encodeCols <- function(item, rank) {
  ind <- which(rank==item)
  
  col <- 0
  if (ind < 31) {
    col <- item
  }
  else if (ind < 100) {
    col <- "RARE"
  }
  else {
    col <- "VERYRARE"
  }
  
  col 
}


## Keep only information about the smallest angle
names(data) <- c("smallestAngle", "ligand1", "ligand2", "bidentate")
data[,2] <- sapply(data[,2], function(x) encodeCols(x, ligRank))
data[,3] <- sapply(data[,3], function(x) encodeCols(x, ligRank))
data[,2:4] <- lapply(data[,2:4] , factor)
str(data)


#### Group the data into normal, compressed, super-compressed, and leave out minangle between 58 and 68.
minAngle <- data[,1]
compress <- function (x) {
  comp <- "normal"
  if (x < as.numeric(args[3]) && x > as.numeric(args[2])) {comp <- "leaveout"}
  else if (x <= as.numeric(args[2]) && x >= 38) {comp <- "compressed"}
  else if (x < 38) {comp <- "supercompressed"}
  
  comp
}

#################### prepare the data for rf ####################

## Training data
groups <- sapply(minAngle, compress)
ind.train <- groups != "leaveout"

data.train <- data[ind.train,]
yy <- groups[ind.train]

## Testing data
ind.test <- groups == "leaveout"
data.test <- data[ind.test,]
data.all <- data

#################### execute random forest ######################

###### sorted middle four ordering
## Train forest
rf.sorted <- randomForest(factor(yy) ~ .,data=data.train)

## Prediction
prediction.test <- predict(rf.sorted, data.test)
prediction.all <- predict(rf.sorted, data.all)

rf.sorted$confusion
save(rf.sorted, prediction.all, file="rf.results.RData")

###### Set simple 63 degree cutoff on test data to compare with the prediction results
compress63 <- function (minAngle) {
  comp <- "normal"
  if (minAngle < 40) {comp <- "supercompressed"}
  else if (minAngle < 63) {comp <- "compressed"}

  comp
}

test.cutoff <- sapply(data[ind.test,1], compress63)
all.cutoff <- sapply(data[,1], compress63)

table(test.cutoff, prediction.test)
table(all.cutoff, prediction.all)



