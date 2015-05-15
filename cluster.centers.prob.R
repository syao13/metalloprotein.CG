####################    load data  ####################  
setwd("~/Desktop/zinc.CG.2015")
load("normal_cluster_assg.RData")
load("compressed_cluster_assg.RData")
load("combined_cluster_assg.RData")

rawdata1 <- read.table("four.chi.txt", header = TRUE)
rawdata2 <- read.table("four.chi.leaveOut.txt", header = TRUE)

## Set the number k normal, compressed, and combined
normal.k <- 10
normal.cluster <- get(paste("normal.", normal.k, ".clusters", sep=""))
normal.cluster[,2] <- as.numeric(normal.cluster[,2])

compressed.k <- 8
compressed.cluster <- get(paste("compressed.", compressed.k, ".clusters", sep=""))
compressed.cluster[,2] <- as.numeric(compressed.cluster[,2])

combined.k <- 14
combined.cluster <- get(paste("combined.", combined.k, ".clusters", sep=""))
combined.cluster[,2] <- as.numeric(combined.cluster[,2])


## Remove duplicate zinc ids
dupId <- rawdata1$Zn_ID[duplicated(rawdata1$Zn_ID)]
rawdata1 <- rawdata1[! rawdata1$Zn_ID %in% dupId,]
dupId2 <- rawdata2$Zn_ID[duplicated(rawdata2$Zn_ID)]
rawdata2 <- rawdata2[! rawdata2$Zn_ID %in% dupId2,]
dim(rawdata1)
dim(rawdata2)

normal.cluster <- normal.cluster[! normal.cluster[,1] %in% dupId,]
dim(normal.cluster)
compressed.cluster <- compressed.cluster[! compressed.cluster[,1] %in% dupId,]
dim(compressed.cluster)
combined.cluster <- combined.cluster[! combined.cluster[,1] %in% dupId,]
dim(combined.cluster)

row.names(rawdata1) <- rawdata1$Zn_ID
row.names(rawdata2) <- rawdata2$Zn_ID

############## cluster centers ###################
#### normal
data.normal <- rawdata1[normal.cluster[,1],]
sortedAngles.normal <- t(apply(data.normal[,2:7], 1, function(x) c(x[1], sort(x[2:5]), x[6])))
colnames(sortedAngles.normal) <- c("angle1", "mid1", "mid2", "mid3", "mid4","angle6")
dim(sortedAngles.normal)
dim(normal.cluster)

round(apply(sortedAngles.normal, 2, function(x) tapply(x, normal.cluster[,2], mean)), digit=1)
round(apply(sortedAngles.normal, 2, function(x) tapply(x, normal.cluster[,2], sd)), digit=1)
round(apply(sortedAngles.normal, 2, function(x) tapply(x, normal.cluster[,2], length))[,1], digit=1)

#### compressed
data.compressed <- rawdata1[compressed.cluster[,1],]
sortedAngles.compressed <- t(apply(data.compressed[,2:7], 1, function(x) c(x[1], sort(x[2:5]), x[6])))
colnames(sortedAngles.compressed) <- c("angle1", "mid1", "mid2", "mid3", "mid4","angle6")
dim(sortedAngles.compressed)
dim(compressed.cluster)

round(apply(sortedAngles.compressed, 2, function(x) tapply(x, compressed.cluster[,2], mean)), digit=1)
round(apply(sortedAngles.compressed, 2, function(x) tapply(x, compressed.cluster[,2], sd)), digit=1)
round(apply(sortedAngles.compressed, 2, function(x) tapply(x, compressed.cluster[,2], length))[,1], digit=1)

#### combined
data.combined <- rawdata1[combined.cluster[,1],]
sortedAngles.combined <- t(apply(data.combined[,2:7], 1, function(x) c(x[1], sort(x[2:5]), x[6])))
colnames(sortedAngles.combined) <- c("angle1", "mid1", "mid2", "mid3", "mid4","angle6")
dim(sortedAngles.combined)
dim(combined.cluster)

round(apply(sortedAngles.combined, 2, function(x) tapply(x, combined.cluster[,2], mean)), digit=1)
round(apply(sortedAngles.combined, 2, function(x) tapply(x, combined.cluster[,2], sd)), digit=1)
round(apply(sortedAngles.combined, 2, function(x) tapply(x, combined.cluster[,2], length))[,1], digit=1)

############## Chi-squared probabilities ###################
#### normal
probs.normal <- data.normal[,34:38]
dim(probs.normal)

round(apply(probs.normal, 2, function(x) tapply(x, normal.cluster[,2], mean)), digit=3)
round(apply(probs.normal, 2, function(x) tapply(x, normal.cluster[,2], sd)), digit=3)

#### compressed
probs.compressed <- data.compressed[,34:38]
dim(probs.compressed)

round(apply(probs.compressed, 2, function(x) tapply(x, compressed.cluster[,2], mean)), digit=3)
round(apply(probs.compressed, 2, function(x) tapply(x, compressed.cluster[,2], sd)), digit=3)

#### compressed and leaving out the smallest angle
data.compressed.leaveout <- rawdata2[compressed.cluster[,1],]
probs.compressed.leaveout <- data.compressed.leaveout[,34:38]
dim(probs.compressed.leaveout)

round(apply(probs.compressed.leaveout, 2, function(x) tapply(x, compressed.cluster[,2], mean)), digit=3)
round(apply(probs.compressed.leaveout, 2, function(x) tapply(x, compressed.cluster[,2], sd)), digit=3)


############## Representative ###################
#### Find the best examples for clusters
## the ids are reordered, and this is only normal or compressed group, so be careful about the ids and row numbers.
centers
cen <- centers[8,] ## change the cluster number
distTOcenter <- apply(abs(t(t(selectAngles)-cen)), 1, sum)
ord <- order(distTOcenter)
ord[1:5]
selectAngles[ord[1:5],]
data[row.names(selectAngles[ord[1:5],]),]
cen

## Print out results
paste(centers, "+/-", centers.sd, sep="")
matrix( data = paste(centers, "+/-", centers.sd, sep=""), ncol=6, byrow=TRUE)

