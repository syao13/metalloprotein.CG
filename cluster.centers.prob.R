#!/usr/bin/Rscript

####################    load data  ####################  
#setwd("../output")
options(stringsAsFactors=FALSE)
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

load("normal_cluster_assg.RData")
load("compressed_cluster_assg.RData")
load("combined_cluster_assg.RData")

rawdata1 <- read.table("r.allLig.txt", header = FALSE)
colnames(rawdata1) <- c("znID", "method", "year", "resolution", "angleCombo", "ligandCombo", "bondlengthCombo", "biStatusCombo", "bfactorCombo", "biLigs", "chainIDCombo", "residueCombo", "atomCombo", "extra")
#rawdata2 <- rawdata1 #read.table("four.chi.leaveOut.txt", header = TRUE)

## Set the number k normal, compressed, and combined
normal.k <- args[3]
normal.cluster <- get(paste("normal.", normal.k, ".clusters", sep=""))
normal.cluster[,2] <- as.numeric(normal.cluster[,2])

compressed.k <- args[4]
compressed.cluster <- get(paste("compressed.", compressed.k, ".clusters", sep=""))
compressed.cluster[,2] <- as.numeric(compressed.cluster[,2])

combined.k <- args[5]
combined.cluster <- get(paste("combined.", combined.k, ".clusters", sep=""))
combined.cluster[,2] <- as.numeric(combined.cluster[,2])


## Remove duplicate zinc ids
dupId <- rawdata1$znID[duplicated(rawdata1$znID)]
rawdata1 <- rawdata1[! rawdata1$znID %in% dupId,]
#dupId2 <- rawdata2$znID[duplicated(rawdata2$zn_ID)]
#rawdata2 <- rawdata2[! rawdata2$znID %in% dupId2,]
dim(rawdata1)
#dim(rawdata2)

normal.cluster <- normal.cluster[! normal.cluster[,1] %in% dupId,]
dim(normal.cluster)
compressed.cluster <- compressed.cluster[! compressed.cluster[,1] %in% dupId,]
dim(compressed.cluster)
combined.cluster <- combined.cluster[! combined.cluster[,1] %in% dupId,]
dim(combined.cluster)

row.names(rawdata1) <- rawdata1$znID
#row.names(rawdata2) <- rawdata2$znID

############## cluster centers ###################
angleSapce <- function(angleCombo, num) {
  angles <- as.numeric(strsplit(angleCombo, ",")[[1]])
  anglesSort <- sort(angles[2:(length(angles)-1)])
  if (num == 5) { c(angles[1], anglesSort[c(1, floor((length(anglesSort) + 1)/2),length(anglesSort))], angles[length(angles)]) }
  else if (num == 6) { c(angles[1], anglesSort[c(1, floor(quantile(1:length(anglesSort), 0.34)), floor(quantile(1:length(anglesSort), 0.67)), length(anglesSort))], angles[length(angles)]) }
  else if (num == 7) { c(angles[1], anglesSort, angles[length(angles)]) }
}

#### normal
angles.norm <- rawdata1[normal.cluster[,1],]$angleCombo
length(angles.norm)
sortedAngles.normal <- t(sapply(angles.norm, function(x) angleSapce(x, args[2])))
rownames(sortedAngles.normal) <- NULL
dim(sortedAngles.normal)
dim(normal.cluster)

print("Table 5")
round(apply(sortedAngles.normal, 2, function(x) tapply(x, normal.cluster[,2], mean)), digit=1)
round(apply(sortedAngles.normal, 2, function(x) tapply(x, normal.cluster[,2], sd)), digit=1)
round(apply(sortedAngles.normal, 2, function(x) tapply(x, normal.cluster[,2], length))[,1], digit=1)

#### compressed
angles.comp <- rawdata1[compressed.cluster[,1],]$angleCombo
length(angles.comp)
sortedAngles.compressed <- t(sapply(angles.comp, function(x) angleSapce(x, args[2])))
rownames(sortedAngles.compressed) <- NULL
dim(sortedAngles.compressed)
dim(compressed.cluster)

print("Table 6")
round(apply(sortedAngles.compressed, 2, function(x) tapply(x, compressed.cluster[,2], mean)), digit=1)
round(apply(sortedAngles.compressed, 2, function(x) tapply(x, compressed.cluster[,2], sd)), digit=1)
round(apply(sortedAngles.compressed, 2, function(x) tapply(x, compressed.cluster[,2], length))[,1], digit=1)

#### combined
angles.all <- rawdata1[combined.cluster[,1],]$angleCombo
length(angles.all)
sortedAngles.combined <- t(sapply(angles.all, function(x) angleSapce(x, args[2])))
rownames(sortedAngles.combined) <- NULL
dim(sortedAngles.combined)
dim(combined.cluster)

print("Table S7")
round(apply(sortedAngles.combined, 2, function(x) tapply(x, combined.cluster[,2], mean)), digit=1)
round(apply(sortedAngles.combined, 2, function(x) tapply(x, combined.cluster[,2], sd)), digit=1)
round(apply(sortedAngles.combined, 2, function(x) tapply(x, combined.cluster[,2], length))[,1], digit=1)

############## Chi-squared probabilities ###################
#### normal
#probs.normal <- data.normal[,34:38]
#dim(probs.normal)

#print("Table 7")
#round(apply(probs.normal, 2, function(x) tapply(x, normal.cluster[,2], mean)), digit=3)
#round(apply(probs.normal, 2, function(x) tapply(x, normal.cluster[,2], sd)), digit=3)

#### compressed
#probs.compressed <- data.compressed[,34:38]
#dim(probs.compressed)

#round(apply(probs.compressed, 2, function(x) tapply(x, compressed.cluster[,2], mean)), digit=3)
#round(apply(probs.compressed, 2, function(x) tapply(x, compressed.cluster[,2], sd)), digit=3)

#### compressed and leaving out the smallest angle
#data.compressed.leaveout <- rawdata2[compressed.cluster[,1],]
#probs.compressed.leaveout <- data.compressed.leaveout[,34:38]
#dim(probs.compressed.leaveout)

#print("Table 8")
#round(apply(probs.compressed.leaveout, 2, function(x) tapply(x, compressed.cluster[,2], mean)), digit=3)
#round(apply(probs.compressed.leaveout, 2, function(x) tapply(x, compressed.cluster[,2], sd)), digit=3)


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

