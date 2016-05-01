#!/mlab/data/software/R-3.2.1-F22/bin/Rscript

###############################################
## Original on all number of ligands together
## Now can handle all number of ligands, given number of ligands and no-heme fe
## args: directory, #lig/all/nonheme, angle space
###############################################

###!/usr/bin/Rscript
####################    load data  ####################  
options(stringsAsFactors=FALSE)
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

load("rf.results.RData")
load("finalMetalList.RData")
data <- read.table("r.allLig.txt", header=FALSE, comment.char = "")
colnames(data) <- c("metalID", "method", "year", "resolution", "angleCombo", "ligandCombo", "bondlengthCombo", "biStatusCombo", "bfactorCombo", "biLigs", "chainIDCombo", "residueCombo", "atomCombo", "amaineN", "occupancy", "solvent")

znList <- data[data[,1] %in% finalZnList & data[,4] < 3, 1]
ligNum <- sapply(data$ligandCombo, function(x) length(strsplit(x, ",")[[1]]))
heme <- sapply(data$ligandCombo, function(x) "HEM" %in% sort(matrix(unlist(strsplit(strsplit(x, ",")[[1]], "[.]")), byrow=TRUE, ncol=3)[,1]))
table(ligNum)

group <- NULL
if (args[2] == "all") {
  group <- rep(1, length(prediction.all))
} else if (args[2] == "noheme") {
  group <- ! heme
} else {
  group <- ligNum == args[2]
}

#### Define the data into normal and compressed from rf prediciton on 58-68 angles.
ind.normal <- (prediction.all=="normal" | prediction.all==2 ) & group 
ind.compress <- (prediction.all=="compressed" | prediction.all==1 )& group 

normal <- data[ind.normal,]
compressed <- data[ind.compress,]
all <- data[ind.normal | ind.compress,]
normal.nr <- normal[normal[,1] %in% znList, ]
compressed.nr <- compressed[compressed[,1] %in% znList, ]
all.nr <- all[all[,1] %in% znList,]
dim(normal.nr) 
dim(compressed.nr) 
dim(all.nr) 

############ define normal vs compressed ############

##### only one compressed
secondComp <- function(angleCombo) {
  angles <- as.numeric(strsplit(angleCombo, ",")[[1]])
  anglesSort <- sort(angles)

  if (anglesSort[2] <= 60) {1}
  else {0}
}

ind.comp2 <- sapply(compressed.nr$angleCombo, secondComp)
compressed.nr <- compressed.nr[!ind.comp2, ]

#### reduced angle space
angleSapce <- function(angleCombo, num, mode="median") {
  angles <- as.numeric(strsplit(angleCombo, ",")[[1]])
  anglesSort <- sort(angles[2:(length(angles)-1)])
  if (num == 5) { c(angles[1], anglesSort[c(1, floor((length(anglesSort) + 1)/2),length(anglesSort))], angles[length(angles)]) }
  else if (num == 6) { c(angles[1], anglesSort[c(1, floor(quantile(1:length(anglesSort), 0.34)), floor(quantile(1:length(anglesSort), 0.67)), length(anglesSort))], angles[length(angles)]) }
  else if (num == 4) { c(angles[1], anglesSort[c(1, 2, length(anglesSort)-1, length(anglesSort))], angles[length(angles)]) }
  else if (num == 3) { c(angles[1], anglesSort[c(1, length(anglesSort)-2, length(anglesSort)-1, length(anglesSort))], angles[length(angles)]) }
  else if (num == 7) { c(angles[1], anglesSort, angles[length(angles)]) }
}

## normal
angles.norm <- normal.nr$angleCombo
length(angles.norm)
angle.sorted.norm <- t(sapply(angles.norm, function(x) angleSapce(x, args[3])))
rownames(angle.sorted.norm) <- NULL

## compressed
angles.comp <- compressed.nr$angleCombo
length(angles.comp)
angle.sorted.comp <- t(sapply(angles.comp, function(x) angleSapce(x, args[3])))
rownames(angle.sorted.comp) <- NULL

## combined
angles.all <- all.nr$angleCombo
length(angles.all)
angle.sorted.all <- t(sapply(angles.all, function(x) angleSapce(x, args[3])))
rownames(angle.sorted.all) <- NULL

####################################################

####################    k-means   ####################    
library(cluster) 
library(clue)

## Test stability with different number of clusters.
## Distance is calculated as sum of (the best match centers' absolute differences)
## Jaccard is calulated as sum of intersection/(i+j-intersetion)
kmeansStab <- function(data, k, nrep){
  fits <- lapply(1:nrep, function(x) kmeans(data, k, iter.max = 30))
  dmat <- matrix(-1, nrow=nrep, ncol=nrep)
  jmat <- matrix(-1, nrow=nrep, ncol=nrep)
  
  for (i in 1:(nrep-1)) {
    for (j in (i+1):nrep) {
      cen1 <- fits[[i]]$centers
      cen2 <- fits[[j]]$centers
      
      clu1 <- fits[[i]]$cluster
      clu2 <- fits[[j]]$cluster
      
      diss <- matrix(0, nrow=k, ncol=k) ## distances of all pair-wise clusters' centers
      for (p in 1:k) {
        for (q in 1:k) {
          diss[p,q] <- sum(abs(cen1[p,] - cen2[q,]))
        }
      }
      
      assign <- solve_LSAP(diss)
      totdis <- sum(sapply(1:k, function(x) diss[x,assign[x]]))
      
      jtab <- table(clu1,clu2)
      jassign <- solve_LSAP(jtab, maximum=TRUE)
      totjac <- sum(sapply(1:k, function(x) 
        jtab[x,jassign[x]]/(sum(clu1==x) + sum(clu2==jassign[x]) - jtab[x,jassign[x]]) ))
      
      dmat[i,j] <- dmat[j,i] <- round(totdis, digits=3)
      jmat[i,j] <- jmat[j,i] <- round(totjac, digits=3)
    }
  }
  ind <- which(rowsum(dmat,rep(1,nrep)) == min(rowsum(dmat,rep(1,nrep))))[1]
  
  avg <- round(apply(data, 2, function(x) tapply(x, fits[[ind]]$cluster, mean)), digit=1)
  std <- round(apply(data, 2, function(x) tapply(x, fits[[ind]]$cluster, sd)), digit=1)
  size <- fits[[ind]]$size
  cluster <- fits[[ind]]$cluster
  
  list(distance=dmat, jaccard=jmat, cluster=cluster, mean=avg, std=std, size=size)
}


#### Better run on lab machine with 500 repetitions.
library(parallel)
sumdiff.norm <- sumdiff.comp <- sumdiff.all <- 0
jaccard.norm <- jaccard.comp <- jaccard.all <- 0
nangle <- dim(angle.sorted.all)[2]
nrep <- 500

date()
resultList <- mclapply(1:30, function(x) kmeansStab(angle.sorted.norm, x, nrep), mc.cores=10)
for (i in 1:length(resultList)){
  k <- length(resultList[[i]]$mean)/nangle

  sumdiff.norm[k] <- mean(resultList[[i]]$distance[resultList[[i]]$distance!=-1])/k
  jaccard.norm[k] <- mean(resultList[[i]]$jaccard[resultList[[i]]$jaccard!=-1])/k
  cluster <- resultList[[i]]$cluster
  clusterAssg <- cbind(normal.nr[,1], cluster)
  assign(paste("normal.", k, ".clusters",  sep=""), clusterAssg)
}

date()
resultList <- mclapply(1:30, function(x) kmeansStab(angle.sorted.comp, x, nrep), mc.cores=10)
for (i in 1:length(resultList)){
  k <- length(resultList[[i]]$mean)/nangle

  sumdiff.comp[k] <- mean(resultList[[i]]$distance[resultList[[i]]$distance!=-1])/k
  jaccard.comp[k] <- mean(resultList[[i]]$jaccard[resultList[[i]]$jaccard!=-1])/k
  cluster <- resultList[[i]]$cluster
  clusterAssg <- cbind(compressed.nr[,1], cluster)
  assign(paste("compressed.", k, ".clusters",  sep=""), clusterAssg)
}

date()
resultList <- mclapply(1:30, function(x) kmeansStab(angle.sorted.all, x, nrep), mc.cores=10)
for (i in 1:length(resultList)){
  k <- length(resultList[[i]]$mean)/nangle

  sumdiff.all[k] <- mean(resultList[[i]]$distance[resultList[[i]]$distance!=-1])/k
  jaccard.all[k] <- mean(resultList[[i]]$jaccard[resultList[[i]]$jaccard!=-1])/k
  cluster <- resultList[[i]]$cluster
  clusterAssg <- cbind(all.nr[,1], cluster)
  assign(paste("combined.", k, ".clusters",  sep=""), clusterAssg)
}
date()

save(list=sapply(1:30, function(x) paste0("normal.", x, ".clusters",  sep="")), file="normal_cluster_assg.RData")
save(list=sapply(1:30, function(x) paste0("compressed.", x, ".clusters",  sep="")), file="compressed_cluster_assg.RData")
save(list=sapply(1:30, function(x) paste0("combined.", x, ".clusters",  sep="")), file="combined_cluster_assg.RData")

save(list=c("sumdiff.norm", "jaccard.norm", "sumdiff.comp", "jaccard.comp", "sumdiff.all", "jaccard.all"), file="two_measures_over_k.RData")
# load("~/Desktop/zinc.CG.2015/two_measures_over_k.RData")

#sumdiff.norm
#plot(sumdiff.norm,type="b", main="Normal, sum of absolute differences", xlab="k")
#jaccard.norm
#plot(jaccard.norm,type="b", main="Normal, Jaccard metric", xlab="k")


#sumdiff.comp
#plot(sumdiff.comp,type="b", main="Compressed, sum of absolute differences", xlab="k")
#jaccard.comp
#plot(jaccard.comp,type="b", main="Compressed, Jaccard metric", xlab="k")

#sumdiff.all
#plot(sumdiff.all,type="b", main="Combined, sum of absolute differences", xlab="k")
#jaccard.all
#plot(jaccard.all,type="b", main="Combined, Jaccard metric", xlab="k")




