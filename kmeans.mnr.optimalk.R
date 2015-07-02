#!/usr/bin/Rscript

date()
## set anglesU to change between normal and compressed angles group
####################    load data  ####################  
setwd("../firstManuscriptOutput")
options(stringsAsFactors=FALSE)

load("rf.sorted.RData")
load("finalZnList.RData")
rawdata <- read.table("four.chi.txt", header = FALSE)
znList <- rawdata[rawdata[,1] %in% finalZnList & rawdata[,24] < 3, 1]

id <- rawdata[,1]
orderid <- order(id)
data <- rawdata[orderid,]

angles <- data[,2:7]
bidentates <- data[,16:21]
ligands <- data[,8:11]

resolution <- data[,24]
anglesU <- angles

#### Define the data into normal and compressed from rf prediciton on 58-68 angles.
minAngle <- apply(angles, 1, min)

ind.normal <- names(prediction.all)[prediction.all=="normal"]
ind.compress <- names(prediction.all)[prediction.all=="compressed"]

angles.normal <- angles[ind.normal,]
angles.compress <- angles[ind.compress,]

normal <- data[ind.normal,]
compressed <- data[ind.compress,]
normal.nr <- normal[normal[,1] %in% znList, ]
compressed.nr <- compressed[compressed[,1] %in% znList, ]
all.nr <- data[data[,1] %in% znList,]
dim(normal.nr) 
dim(compressed.nr) 
dim(all.nr) 

############ define normal vs compressed ############

##### only one compressed
sortedA <- t(apply(compressed.nr[,2:7], 1, sort))
ind.comp2 <- sortedA[,2] <= 63
compressed.nr <- compressed.nr[!ind.comp2, ]
angles.comp <- compressed.nr[,2:7]
dim(angles.comp)
angle.sorted.comp <- t(apply(angles.comp, 1, function(x) c(x[1], sort(x[2:5]), x[6])))
colnames(angle.sorted.comp) <- c("angle1", "mid1", "mid2", "mid3", "mid4","angle6")

##### normal
angles.norm <- normal.nr[,2:7]
dim(angles.norm)
angle.sorted.norm <- t(apply(angles.norm, 1, function(x) c(x[1], sort(x[2:5]), x[6])))
colnames(angle.sorted.norm) <- c("angle1", "mid1", "mid2", "mid3", "mid4","angle6")

##### combined
angles.all <- all.nr[,2:7]
dim(angles.all)
angle.sorted.all <- t(apply(angles.all, 1, function(x) c(x[1], sort(x[2:5]), x[6])))
colnames(angle.sorted.all) <- c("angle1", "mid1", "mid2", "mid3", "mid4","angle6")
#####################################################

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
      
      dmat[i,j] <- round(totdis, digits=3)
      jmat[i,j] <- round(totjac, digits=3)
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
nrep <- 500 
sumdiff.norm <- sumdiff.comp <- sumdiff.all <- 0
jaccard.norm <- jaccard.comp <- jaccard.all <- 0

for (k in 1:30) {
  kmat.norm <- kmeansStab(angle.sorted.norm, k, nrep)
  sumdiff.norm[k] <- mean(kmat.norm$distance[kmat.norm$distance!=-1])/k
  jaccard.norm[k] <- mean(kmat.norm$jaccard[kmat.norm$jaccard!=-1])/k
  cluster <- kmat.norm$cluster
  clusterAssg <- cbind(normal.nr[,1], cluster)
  assign(paste("normal.", k, ".clusters",  sep=""), clusterAssg)
  
  kmat.comp <- kmeansStab(angle.sorted.comp, k, nrep)
  sumdiff.comp[k] <- mean(kmat.comp$distance[kmat.comp$distance!=-1])/k
  jaccard.comp[k] <- mean(kmat.comp$jaccard[kmat.comp$jaccard!=-1])/k
  cluster <- kmat.comp$cluster
  clusterAssg <- cbind(compressed.nr[,1], cluster)
  assign(paste("compressed.", k, ".clusters",  sep=""), clusterAssg)
  
  kmat.all <- kmeansStab(angle.sorted.all, k, nrep)
  sumdiff.all[k] <- mean(kmat.all$distance[kmat.all$distance!=-1])/k
  jaccard.all[k] <- mean(kmat.all$jaccard[kmat.all$jaccard!=-1])/k
  cluster <- kmat.all$cluster
  clusterAssg <- cbind(all.nr[,1], cluster)
  assign(paste("combined.", k, ".clusters",  sep=""), clusterAssg)
}

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

date()



