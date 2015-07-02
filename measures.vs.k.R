#!/usr/bin/Rscript

### Load the data
library(cluster) 
library(clue)

options(stringsAsFactors=FALSE)
setwd("../firstManuscriptOutput")

rawdata <- read.table("four.chi.txt", header = TRUE)

############ normal vs. compressed#####################
### choose one to operate 
## normal
load("./normal_cluster_assg.RData")
data.normal <- rawdata[sapply(normal.1.clusters[,1], function(x) which(rawdata[,1] == x)[1]),]
angles.normal <- t(apply(data.normal[2:7], 1, function(x) c(x[1], sort(x[2:5]), x[6])))
colnames(angles.normal) <- c("angle1", "mid1", "mid2", "mid3", "mid4","angle6")


## compressed
load("compressed_cluster_assg.RData")
data.compressed <- rawdata[sapply(compressed.1.clusters[,1], function(x) which(rawdata[,1] == x)[1]),]
angles.compressed <- t(apply(data.compressed[2:7], 1, function(x) c(x[1], sort(x[2:5]), x[6])))
colnames(angles.compressed) <- c("angle1", "mid1", "mid2", "mid3", "mid4","angle6")


## combined
load("combined_cluster_assg.RData")
data.combined <- rawdata[sapply(combined.1.clusters[,1], function(x) which(rawdata[,1] == x)[1]),]
angles.combined <- t(apply(data.combined[2:7], 1, function(x) c(x[1], sort(x[2:5]), x[6])))
colnames(angles.combined) <- c("angle1", "mid1", "mid2", "mid3", "mid4","angle6")

###### Compute all distance matrices for all k
getDistM <- function(selectAngles, clusters) {
  clusterNum <- length(unique(clusters))
  dist.mat <- matrix(0, nrow=clusterNum, ncol=clusterNum) 
  
  for (i in 1:clusterNum) {
    for (j in i:clusterNum) {
      cluster1 <- selectAngles[clusters == i,]
      cluster2 <- selectAngles[clusters == j,]
      
      rmsd <- 0
      for (m in 1: nrow(cluster1)) {
        for (n in 1: nrow(cluster2)) {
          angles1 <- cluster1[m,]
          angles2 <- cluster2[n,]
          
          rmsd <- rmsd + sqrt(sum((angles1-angles2)^2))
        }
      }
      dist.mat[i,j] <- dist.mat[j,i] <- rmsd / nrow(cluster1) / nrow(cluster2)      
    }
  }
  dist.mat
}

## This takes an hour or so on my laptop

for (i in 1:30) {
  clusters <- get(paste("normal.",i,".clusters", sep=""))
  assign(paste("dist.struct.norm.", i,  sep=""), getDistM(angles.normal, as.numeric(clusters[,2])))
  
  clusters <- get(paste("compressed.",i,".clusters", sep=""))
  assign(paste("dist.struct.comp.", i,  sep=""), getDistM(angles.compressed, as.numeric(clusters[,2])))

  clusters <- get(paste("combined.",i,".clusters", sep=""))
  assign(paste("dist.struct.comb.", i,  sep=""), getDistM(angles.combined, as.numeric(clusters[,2])))
}

save(list=sapply(1:30, function(x) paste("dist.struct.norm.", x,  sep="")), file="dists_struct_normal.RData")
save(list=sapply(1:30, function(x) paste("dist.struct.comp.", x,  sep="")), file="dists_struct_compressed.RData")
save(list=sapply(1:30, function(x) paste("dist.struct.comb.", x,  sep="")), file="dists_struct_combined.RData")

# for (x in 1:30) {assign(paste("dist.struct.norm.", x, sep=""), get(paste("dist.norm.", x, sep="")))}
# for (x in 1:30) {assign(paste("dist.struct.comp.", x, sep=""), get(paste("dist.comp.", x, sep="")))}
# save(list=sapply(1:30, function(x) paste("dist.struct.norm.", x,  sep="")), file="dists_struct_normal.RData")
# save(list=sapply(1:30, function(x) paste("dist.struct.comp.", x,  sep="")), file="dists_struct_compressed.RData")


####################  rho/p-value vs k #################### 
# load("dists_struct_normal.RData")
# load("dists_struct_compressed.RData")
# load("dists_struct_combined.RData")

load("normal_funct_dist.RData")
load("compressed_funct_dist.RData")
load("combined_funct_dist.RData")

just_dist <- function(dist_list){lapply(dist_list, function(x){x$dist})}
compare_list_distances <- function(dist_a, dist_b, method = "spearman"){
  out_cor <- lapply(names(dist_a), function(in_dist){
    s_dist <- as.vector(dist_a[[in_dist]])
    
    pos_funct <- grep(in_dist, names(dist_b), value = TRUE)
    s_f_cor <- sapply(dist_b[pos_funct], function(in_funct){
      tmp <- cor.test(s_dist, as.vector(in_funct), alternative = "two.sided", method = method)
      c(tmp$estimate, tmp$p.value)
    })
  })
  
  names(out_cor) <- names(dist_a)
  return(out_cor)
}

dist.func.norm <- just_dist(normal_funct_dist)
dist.func.comp <- just_dist(compressed_funct_dist)
dist.func.comb <- just_dist(combined_funct_dist)

compareDist <- function(dist.a, dist.b, method="spearman") {
  vec.a <- as.vector(dist.a)
  vec.b <- as.vector(dist.b)
  
  ind <- vec.b!=0
  tmp <- cor.test(vec.a[ind], vec.b[ind], alternative="two.sided", method = method)
  c(tmp$estimate, tmp$p.value)
}


rhopVsk <- function(nLig, nc="norm", kVec) {
  rho <- NULL
  pval <- NULL
  for (i in kVec) {
    dist.struc <- as.dist(get(paste("dist.struct.", nc, ".", i, sep="")))
    
    if (nc == "norm") {name <- paste(nLig, ",normal.", i, ".clusters", sep="") }
    if (nc == "comp") {name <- paste(nLig, ",compressed.", i, ".clusters", sep="") }
    if (nc == "comb") {name <- paste(nLig, ",combined.", i, ".clusters", sep="") }
    
    dist.func <- get(paste("dist.func.", nc, sep=""))
    dist.func <- dist.func[names(dist.func)==name][[1]]
    
    temp <- compareDist(dist.struc, dist.func)
    rho <- c(rho, temp[1])
    pval <- c(pval, temp[2])
    
    #plot(dist.struc, dist.func)
    }
  
  list(rho=rho, p_value=pval)
}

kVec <- 3:30
rhop.1.norm <- rhopVsk(1, "norm", kVec)
rhop.2.norm <- rhopVsk(2, "norm", kVec)
rhop.3.norm <- rhopVsk(3, "norm", kVec)
rhop.4.norm <- rhopVsk(4, "norm", kVec)

rhop.1.comp <- rhopVsk(1, "comp", kVec)
rhop.2.comp <- rhopVsk(2, "comp", kVec)
rhop.3.comp <- rhopVsk(3, "comp", kVec)
rhop.4.comp <- rhopVsk(4, "comp", kVec)

rhop.1.all <- rhopVsk(1, "comb", kVec)
rhop.2.all <- rhopVsk(2, "comb", kVec)
rhop.3.all <- rhopVsk(3, "comb", kVec)
rhop.4.all <- rhopVsk(4, "comb", kVec)


load("two_measures_over_k.RData")
save(list=c("sumdiff.norm", "jaccard.norm", "rhop.1.norm", "rhop.2.norm", "rhop.3.norm", "rhop.4.norm",
            "sumdiff.comp", "jaccard.comp", "rhop.1.comp", "rhop.2.comp", "rhop.3.comp", "rhop.4.comp",
            "sumdiff.all", "jaccard.all", "rhop.1.all", "rhop.2.all", "rhop.3.all", "rhop.4.all"), 
     file="four_measures_over_k.RData")





