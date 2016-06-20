#!/usr/bin/Rscript

###########################################################################################
#
##   Written by Sen Yao, 07/09/2015
##   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/09/2015. All rights reserved.
##
##   Usage: ./center.prob.metallist.R directory normal_k compressed_k combined_k angle_space(6 for individual ligNum, 7 for combined)
##
###########################################################################################

##################    load data  ####################  
#setwd("../output")
options(stringsAsFactors=FALSE)
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

load("normal_cluster_assg.RData")
load("compressed_cluster_assg.RData")
load("combined_cluster_assg.RData")

Tet <- c(109.5, 109.5, 109.5, 109.5, 109.5, 109.5)
Bva <- c(120, 90, 90,120, 120, 90)
Bvp <- c(180, 90, 90, 90, 90, 120)
Pyv <- c(180, 90, 90, 90, 90, 90)
Spl <- c(180, 90, 90, 90, 90, 180)

Tbp <- c(180, 90, 90, 90, 90, 90, 90, 120, 120, 120)
Spy <- c(180, 90, 90, 90, 90, 90, 90, 90, 180, 90)
Tpv <- c(131.8, 70.6, 90, 90, 90, 90, 131.8, 131.8, 131.8, 70.6)
#Ppl <- c(144, 72, 72, 72, 72, 144, 144, 144, 144, 72)

Oct <- c(180, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 180, 180, 90)
Pva <- c(144, 72, 72, 72, 72, 90, 90, 90, 90, 90, 144, 144, 144, 144, 72)
Pvp <- c(180, 72, 72, 90, 90, 90, 90, 90, 90, 90, 90, 144, 144, 144, 72)
Tpr <- c(131.8, 70.6, 70.6, 90, 90, 90, 90, 90, 90, 131.8, 131.8, 131.8, 131.8, 131.8, 70.6)

Pbp <- c(180, 72, 72, 72, 72, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 144, 144, 144, 144, 144, 72)
Hva <- c(180, 60, 60, 60, 60, 60, 90, 90, 90, 90, 90, 90, 120, 120, 120, 120, 120, 120, 180, 180, 60)
Hvp <- c(180, 60, 60, 60, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 120, 120, 120, 120, 180, 180, 60)
Cuv <- c(180, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 180, 180, 70.5)
Sav <- c(143.6, 70.5, 70.5, 70.5, 70.5, 70.5, 82, 82, 82, 82, 82, 82, 109.5, 109.5, 109.5, 143.6, 143.6, 143.6, 143.6, 143.6, 70.5)

Hbp <- c(180, 60, 60, 60, 60, 60, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 120, 120, 120, 120, 120, 120, 120, 180, 180, 60)
Cub <- c(180, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 109.5, 180, 180, 180, 70.5)
Sqa <- c(143.6, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 82, 82, 82, 82, 82, 82, 82, 82, 109.5, 109.5, 109.5, 109.5, 143.6, 143.6, 143.6, 143.6, 143.6, 143.6, 143.6, 70.5)

rawdata <- read.table("r.allLig.txt", header = FALSE, comment.char = "")
colnames(rawdata) <- c("znID", "method", "year", "resolution", "angleCombo", "ligandCombo", "bondlengthCombo", "biStatusCombo", "bfactorCombo", "biLigs", "chainIDCombo", "residueCombo", "atomCombo", "extra")
row.names(rawdata) <- rawdata$znID
fullAngleSpace <- 7
############## cluster centers ###################
angleSapce <- function(angleCombo, num) {
  angles <- as.numeric(strsplit(angleCombo, ",")[[1]])
  anglesSort <- sort(angles[2:(length(angles)-1)])
  if (num == 5) { c(angles[1], anglesSort[c(1, floor((length(anglesSort) + 1)/2),length(anglesSort))], angles[length(angles)]) }
  else if (num == 6) { c(angles[1], anglesSort[c(1, floor(quantile(1:length(anglesSort), 0.34)), floor(quantile(1:length(anglesSort), 0.67)), length(anglesSort))], angles[length(angles)]) }
  else if (num == 7) { c(angles[1], anglesSort, angles[length(angles)]) }
}

rmsdCG <- function(cc) {
  CGs <- NULL
  if (ncol(cc) == 6) {
    CGs <- c("Tet", "Bva", "Bvp", "Pyv", "Spl")
  } else if (ncol(cc) == 10) {
    CGs <- c("Tbp", "Spy", "Tpv")
  } else if (ncol(cc) == 15) {
    CGs <- c("Oct", "Pva", "Pvp", "Tpr")
  } else if (ncol(cc) == 21) {
    CGs <- c("Pbp", "Hva", "Hvp", "Cuv", "Sav")
  } else if (ncol(cc) == 28) {
    CGs <- c("Hbp", "Cub", "Sqa")
  }

  sapply(CGs, function(x) sqrt(rowMeans((sweep(cc, 2, get(x)))^2)))
}

genCSVtable <- function(rf, k) {
  cluster <- get(paste(rf, ".", k, ".clusters", sep=""))
  cluster[,2] <- as.numeric(cluster[,2])

  angles <- rawdata[cluster[,1],]$angleCombo
  #length(angles)
  sortedAngles <- t(sapply(angles, function(x) angleSapce(x, args[5])))
  rownames(sortedAngles) <- NULL
  #dim(sortedAngles)
  #dim(cluster)

  centers <- round(apply(sortedAngles, 2, function(x) tapply(x, cluster[,2], mean)), digit=1)
  #centers
  #round(rmsdCG(clusterCenters),1)
  std <- round(apply(sortedAngles, 2, function(x) tapply(x, cluster[,2], sd)), digit=1)
  size <- round(apply(sortedAngles, 2, function(x) tapply(x, cluster[,2], length))[,1], digit=1)
  meanstd <- matrix(paste(centers, "+/-", std, sep=""), nrow = k)
 
  nAngle <- ncol(meanstd)
  star <- rep("", nAngle-2) 
  star[c(1, floor(quantile(1:(nAngle-2), 0.34)), floor(quantile(1:(nAngle-2), 0.67)), nAngle-2)] <- "*"
  colnames(meanstd) <- c("largest_angle*", paste("middle_", 1:(nAngle-2), star, sep=""), "smallest_opposite_angle*")

  #cbind(size, meanstd)
  numLig <- rawdata[cluster[1,1],17]
  idx <- NULL
  if (args[5] == 6) {
    idx <- 1:18
  } else if (numLig == "four") { 
    idx <- c(1,3,4,7,8)
  } else if (numLig == "five") {
    idx <- c(2,6,10)
  } else if (numLig == "six") {
    idx <- c(5,9,12,13)
  } else if (numLig == "seven") {
    idx <- c(11,15,17,18)
  } else if (numLig == "eight") {
    idx <- c(14,16)
  }

  probs <- read.table("../ia.cutoff/probs.txt", header = TRUE)
  if (rf == "compressed") {probs <- read.table("../ia.cutoff/probs.leaveOut.txt", header = TRUE)}

  probs <- probs[! duplicated(probs[,1]), ]
  rownames(probs) <- probs[,1]
  probs <- probs[cluster[,1],]
  avgprobs <- round(apply(probs[2:19], 2, function(x) tapply(x, cluster[,2], mean)), digit=3)[,idx]

  probTable <- cbind(size, meanstd, avgprobs)
  write.csv(probTable[order(as.numeric(rownames(probTable))), ], file=paste("average_prob.", rf, ".csv", sep=""))
#}

#getMembers <- function(rf, k) {
#  cluster <- get(paste(rf, ".", k, ".clusters", sep=""))
#  cluster[,2] <- as.numeric(cluster[,2])

  members <- NULL
  for (i in 1:k) {
    members <- rbind(members, c(i, toString(cluster[cluster[,2] == i, 1])))
  }
  colnames(members) <- c("cluster", "members")
  write.table(members, file = paste("cluster_members.", rf, ".txt", sep=""))
}

#if (substring(args[1], 14, 19) != "allLig") {
  genCSVtable("normal", as.numeric(args[2]))
  genCSVtable("compressed", as.numeric(args[3]))
  genCSVtable("combined", as.numeric(args[4]))
#} 

#getMembers("normal", as.numeric(args[2]))
#getMembers("compressed", as.numeric(args[3]))
#getMembers("combined", as.numeric(args[4]))

############## Representative ###################
#### Find the best examples for clusters
## the ids are reordered, and this is only normal or compressed group, so be careful about the ids and row numbers.
#centers
#cen <- centers[8,] ## change the cluster number
#distTOcenter <- apply(abs(t(t(selectAngles)-cen)), 1, sum)
#ord <- order(distTOcenter)
#ord[1:5]
#selectAngles[ord[1:5],]
#data[row.names(selectAngles[ord[1:5],]),]
#cen

## Print out results
#paste(centers, "+/-", centers.sd, sep="")
#matrix( data = paste(centers, "+/-", centers.sd, sep=""), ncol=6, byrow=TRUE)

