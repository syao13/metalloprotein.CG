#!/usr/bin/Rscript

###########################################################################################
#
##   Written by Sen Yao, 07/09/2015
##   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/09/2015. All rights reserved.
##
##   Usage: ./hierarchical.clustering.R directory normal_k normal_n_lig compressed_k compressed_n_lig
###########################################################################################

library(ggplot2)
options(stringsAsFactors=FALSE)
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

load("dists_struct_normal.RData")
load("dists_struct_compressed.RData")
load("normal_funct_dist.RData")
load("compressed_funct_dist.RData")


##########################
## normal group, k=10
##########################
normal.k <- args[2] 
struc.dist.norm <- get(paste("dist.struct.norm.", normal.k, sep=""))
nm <- paste(args[3], ",normal.",normal.k,".clusters", sep="")
funct.dist.norm <- normal_funct_dist[nm][[1]]

dist.str <- as.dist(struc.dist.norm)
dist.fun <- funct.dist.norm$dist
vec.a <- as.vector(dist.str)
vec.b <- as.vector(dist.fun)
ind <- ! is.na(vec.a) & ! is.na(vec.b) # vec.a!=0 & vec.b!=0 

distances <- data.frame(func = vec.b[ind], struc = vec.a[ind])
#distances <- data.frame(func = as.vector(dist.fun), struc = as.vector(dist.str))

png(file="hierarchical.clustering.normal.png", units="in", width=10, height=5, res=200)
par(mfrow=c(1,2)) 
hc.str <- hclust(dist.str)
plot(hc.str, main = paste("Normal, k=", normal.k, sep=""), ylab = "Structural Cluster Distances", xlab = "Cluster", 
     sub = "", cex.lab = 1.5, cex.main=1.5, cex=1.3)

if (sum(is.na(dist.fun))) {
  aa <- as.matrix(dist.fun)
  idx <- apply(aa, 1, sum, na.rm=TRUE)
  dist.fun <- as.dist(aa[idx != 0, idx != 0])
  }

hc.fun <- hclust(dist.fun)
plot(hc.fun, main = paste("Normal, k=", normal.k, sep=""), ylab = "Functional Cluster Distances", xlab = "Cluster", 
     sub = "", cex.lab = 1.5, cex.main=1.5, cex=1.3)
dev.off()

png(file="funct_struct_distance.normal.png", units="in", width=10, height=10, res=200)
tmp <- cor.test(distances$struc, distances$func, alternative="two.sided", method = "spearman")
corVal <- tmp$estimate
pVal <- tmp$p.value
corVal
pVal
ggplot(distances, aes(x = struc, y = func)) + 
  geom_point(size = 4) + xlab("Structural") + ylab("Functional") + ggtitle(paste("Cluster Distances, normal ", normal.k, sep="")) +   
  geom_text(data = NULL, aes(family="serif"), size = 8, x=Inf, y = -Inf, vjust=-1, hjust=1,
            label = paste("rho = ", substr(as.character(corVal), 1, 4), " \n ", 
                          "p-value = ", signif(pVal,3), sep = "")) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), plot.title = element_text(size = 20))
#corVal <- cor(distances$struc, distances$func, method = "spearman")
#corVal
#pVal <- cor.test(distances$struc, distances$func, method = "spearman")$p.value
#pVal
#ggplot(distances, aes(x = struc, y = func)) + 
#  geom_point(size = 4) + xlab("Structural") + ylab("Functional") + ggtitle(paste("Cluster Distances, normal ", normal.k, sep="")) + 
#  geom_text(data = NULL, aes(family="serif"), x = Inf, y = Inf, 
#            label = paste("rho = ", substr(as.character(corVal), 1, 4), " \n ", 
#                          "p-value = ", signif(pVal,3), sep = ""), size = 10) + 
#  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), plot.title = element_text(size = 20))
dev.off()


##########################
## compressed group
##########################
## this is what we have seen before.
compressed.k <- args[4]
struc.dist.comp <- get(paste("dist.struct.comp.", compressed.k, sep=""))
nm <- paste(args[5], ",compressed.",compressed.k,".clusters", sep="")
funct.dist.comp <- compressed_funct_dist[nm][[1]]

dist.str <- as.dist(struc.dist.comp)
dist.fun <- funct.dist.comp$dist
vec.a <- as.vector(dist.str)
vec.b <- as.vector(dist.fun)
ind <- ! is.na(vec.a) & ! is.na(vec.b) # vec.a!=0 & vec.b!=0 

distances <- data.frame(func = vec.b[ind], struc = vec.a[ind])
#distances <- data.frame(func = as.vector(dist.fun), struc = as.vector(dist.str))

png(file="hierarchical.clustering.compressed.png", units="in", width=10, height=5, res=200)
par(mfrow=c(1,2))
hc.str <- hclust(dist.str)
plot(hc.str, main = paste("Compressed, k=", compressed.k, sep=""), ylab = "Structural Cluster Distances", xlab = "Cluster", 
     sub = "", cex.lab = 1.5, cex.main=1.5, cex=1.3)

if (sum(is.na(dist.fun))) {
  aa <- as.matrix(dist.fun)
  idx <- apply(aa, 1, sum, na.rm=TRUE)                                                                                                                                
  dist.fun <- as.dist(aa[idx != 0, idx != 0])
  }

hc.fun <- hclust(dist.fun)
plot(hc.fun, main = paste("Compressed, k=", compressed.k, sep=""), ylab = "Functional Cluster Distances", xlab = "Cluster", 
     sub = "", cex.lab = 1.5, cex.main=1.5, cex=1.3)
dev.off()

png(file="funct_struct_distance.compressed.png", units="in", width=10, height=10, res=200)
tmp <- cor.test(distances$struc, distances$func, alternative="two.sided", method = "spearman")
corVal <- tmp$estimate
pVal <- tmp$p.value
corVal
pVal
ggplot(distances, aes(x = struc, y = func)) + 
  geom_point(size = 4) + xlab("Structural") + ylab("Functional") + ggtitle(paste("Cluster Distances, compressed ", compressed.k, sep="")) +   
  geom_text(data = NULL, aes(family="serif"), size = 8, x=Inf, y = -Inf, vjust=-1, hjust=1,
            label = paste("rho = ", substr(as.character(corVal), 1, 4), " \n ", 
                          "p-value = ", signif(pVal,3), sep = "")) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), plot.title = element_text(size = 20))
#corVal <- cor(distances$struc, distances$func, method = "spearman")
#corVal
#pVal <- cor.test(distances$struc, distances$func, method = "spearman")$p.value
#pVal
#ggplot(distances, aes(x = struc, y = func)) + 
#  geom_point(size = 4) + xlab("Structural") + ylab("Functional") + ggtitle(paste("Cluster Distances, compressed ", compressed.k, sep="")) + 
#  geom_text(data = NULL, aes(family="serif"), x = Inf, y = Inf, 
#            label = paste("rho = ", substr(as.character(corVal), 1, 4), " \n ", 
#                          "p-value = ", signif(pVal,3), sep = ""), size = 10) + 
#  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), plot.title = element_text(size = 20))
dev.off()
