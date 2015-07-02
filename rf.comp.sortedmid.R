#!/usr/bin/Rscript

#### Load the data
library(randomForest)

options(stringsAsFactors=FALSE)
setwd("../output")
rawdata <- read.table("four.chi.txt", header = TRUE)

id <- rawdata[,1]
orderid <- order(id)
data <- rawdata[orderid,]

angles <- data[,2:7]
bidentates <- data[,16:21]
ligands <- data[,8:11]

########## ligand Encoding ########### 
ligs <- c(ligands[,1], ligands[,2], ligands[,3], ligands[,4])
#length(table(ligs))
sort(table(ligs),decreasing=TRUE)[1:100]
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

########## sorted middle four ordering ########## 
lig1cols <- sapply(ligands[,1], function(x) encodeCols(x, ligRank))
lig2cols <- sapply(ligands[,2], function(x) encodeCols(x, ligRank))
lig3cols <- sapply(ligands[,3], function(x) encodeCols(x, ligRank))
lig4cols <- sapply(ligands[,4], function(x) encodeCols(x, ligRank))

lig12 <- cbind(lig1cols, lig2cols)
lig34 <- cbind(lig3cols, lig4cols)
lig12sorted <- t(apply(lig12, 1, sort))
lig34sorted <- t(apply(lig34, 1, sort))

## Sort bi status according to angle orders
sortedBi <- t(sapply(1:nrow(bidentates), 
                     function(x) c(bidentates[x,1], bidentates[x,(order(angles[x,2:5])+1)], bidentates[x,6])))
biInOne<- apply(sortedBi, 1, function(x) paste(x[1], x[2],x[3],x[4],x[5],x[6], sep=""))
sortedAngles <- t(apply(angles, 1, function(x) c(x[1], sort(x[2:5]), x[6])))


data.all <- as.data.frame(cbind(sortedAngles, biInOne, lig12sorted, lig34sorted))
data.all[,1:6] <- lapply(data.all[,1:6] , as.numeric)
data.all[,8:11] <- lapply(data.all[,8:11] , factor)
names(data.all) <- c("angle1", "angle2", "angle3", "angle4", "angle5", "angle6", "bidentate",
                        "ligand1", "ligand2", "ligand3", "ligand4")
str(data.all)


#### Group the data into normal, compressed, super-compressed, and leave out minangle between 58 and 68.
minAngle <- apply(angles, 1, min)
compress <- function (x) {
  comp <- "normal"
  if (x < 75 && x > 58) {comp <- "leaveout"}
  else if (x <= 58 && x >= 38) {comp <- "compressed"}
  else if (x < 38) {comp <- "supercompressed"}
  
  comp
}

#################### prepare the data for rf ####################

## Training data
groups <- sapply(minAngle, compress)
ind.train <- groups != "leaveout"
data.train.sorted <- data.all[ind.train,]

yy <- groups[ind.train]

## Testing data
ind.test <- groups == "leaveout"
data.test <- data.all[ind.test,]

#################### execute random forest ######################

###### sorted middle four ordering
## Train forest
rf.sorted <- randomForest(factor(yy) ~ .,data=data.train.sorted)

## Prediction
prediction.test <- predict(rf.sorted, data.test)
prediction.all <- predict(rf.sorted, data.all)

rf.sorted$confusion
save(rf.sorted, prediction.all, file="rf.sorted.RData")

###### Set simple 63 degree cutoff on test data to compare with the prediction results
compress63 <- function (x) {
  comp <- "normal"
  minAngle <- min(x)
  if (minAngle < 40) {comp <- "supercompressed"}
  else if (minAngle < 63) {comp <- "compressed"}

  comp
}

test.cutoff <- apply(data.test[,1:6], 1, compress63)
all.cutoff <- apply(data.all[,1:6], 1, compress63)

table(test.cutoff, prediction.test)
table(all.cutoff, prediction.all)



