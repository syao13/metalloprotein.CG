#!/usr/bin/Rscript

###########################################################################################
#
##   Written by Sen Yao, 07/09/2015
##   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/09/2015. All rights reserved.
#
###########################################################################################

library(ggplot2)
options(stringsAsFactors=FALSE)
setwd("../output")

load("dists_struct_normal.RData")
load("dists_struct_compressed.RData")
load("normal_funct_dist.RData")
load("compressed_funct_dist.RData")


##########################
## normal group, k=10
##########################
normal.k <- 9
struc.dist.norm <- get(paste("dist.struct.norm.", normal.k, sep=""))
nm <- paste("4,normal.",normal.k,".clusters", sep="")
funct.dist.norm <- normal_funct_dist[nm][[1]]

dist.str <- as.dist(struc.dist.norm)
dist.fun <- funct.dist.norm$dist
distances <- data.frame(func = as.vector(dist.fun), struc = as.vector(dist.str))

png(file="Figure9_1.png", units="in", width=10, height=5, res=300)
par(mfrow=c(1,2)) 
hc.str <- hclust(dist.str)
plot(hc.str, main = "Normal, k=9", ylab = "Structural Cluster Distances", xlab = "Cluster", 
     sub = "", cex.lab = 1.5, cex.main=1.5, cex=1.3)

hc.fun <- hclust(dist.fun)
plot(hc.fun, main = "Normal, k=9", ylab = "Functional Cluster Distances", xlab = "Cluster", 
     sub = "", cex.lab = 1.5, cex.main=1.5, cex=1.3)
dev.off()

png(file="Figure9_2.png", units="in", width=10, height=10, res=300)
corVal <- cor(distances$struc, distances$func, method = "spearman")
corVal
pVal <- cor.test(distances$struc, distances$func, method = "spearman")$p.value
pVal
ggplot(distances, aes(x = struc, y = func)) + 
  geom_point(size = 4) + xlab("Structural") + ylab("Functional") + ggtitle("Cluster Distances, normal 9") + 
  geom_text(data = NULL, aes(family="serif"), x = 60, y = 0.97, 
            label = paste("rho = ", substr(as.character(corVal), 1, 4), " \n ", 
                          "p-value = ", signif(pVal,3), sep = ""), size = 10) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), plot.title = element_text(size = 20))
dev.off()


##########################
## compressed group, k=8
##########################
## this is what we have seen before.
compressed.k <- 10
struc.dist.comp <- get(paste("dist.struct.comp.", compressed.k, sep=""))
nm <- paste("4,compressed.",compressed.k,".clusters", sep="")
funct.dist.comp <- compressed_funct_dist[nm][[1]]

dist.str <- as.dist(struc.dist.comp)
dist.fun <- funct.dist.comp$dist
distances <- data.frame(func = as.vector(dist.fun), struc = as.vector(dist.str))

png(file="Figure10_1.png", units="in", width=10, height=5, res=300)
par(mfrow=c(1,2))
hc.str <- hclust(dist.str)
plot(hc.str, main = "Compressed, k=10", ylab = "Structural Cluster Distances", xlab = "Cluster", 
     sub = "", cex.lab = 1.5, cex.main=1.5, cex=1.3)

hc.fun <- hclust(dist.fun)
plot(hc.fun, main = "Compressed, k=10", ylab = "Functional Cluster Distances", xlab = "Cluster", 
     sub = "", cex.lab = 1.5, cex.main=1.5, cex=1.3)
dev.off()

png(file="Figure10_2.png", units="in", width=10, height=10, res=300)
corVal <- cor(distances$struc, distances$func, method = "spearman")
corVal
pVal <- cor.test(distances$struc, distances$func, method = "spearman")$p.value
pVal
ggplot(distances, aes(x = struc, y = func)) + 
  geom_point(size = 4) + xlab("Structural") + ylab("Functional") + ggtitle("Cluster Distances, compressed 10") + 
  geom_text(data = NULL, aes(family="serif"), x = 70, y = 0.97, 
            label = paste("rho = ", substr(as.character(corVal), 1, 4), " \n ", 
                          "p-value = ", signif(pVal,3), sep = ""), size = 10) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18), plot.title = element_text(size = 20))
dev.off()
