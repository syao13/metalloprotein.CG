#!/mlab/data/software/R-3.2.1-F22/bin/Rscript
## !/usr/bin/Rscript
############################################################################################
##
##   Written by Sen Yao, 07/20/2016
##   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/20/2016. All rights reserved.
##
###########################################################################################

################
## simulation
################

library(IRanges)
library(atom2seq)
library(Biobase)
library(magrittr)
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
load("non_redundant.RData")

interproRes <- read.table("non_redundant.ipr.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE, quote="", comment.char="")
names(interproRes) <- c("id", "md5", "length", "analysis", "sigAccession", "sigDescription", "start", "stop", "score", "status", "date", "iprAccession", "iprDescription", "GO", "Pathways")

iprData <- unique(interproRes[, c("iprAccession", "iprDescription")])
rownames(iprData) <- iprData$iprAccession

# split this into a list based on the ID string, and then transform to ranges with associated meta data. 
splitInter <- split(interproRes, interproRes$id)
splitRanges_all <- lapply(splitInter, interpro2ranges, unique_on = c("sigAccession", "start", "stop"), require_char = "sigAccession")
splitRanges_ipr <- lapply(splitInter, interpro2ranges, unique_on = c("iprAccession", "start", "stop"), require_char = "iprAccession")

# These are what we will work with for all subsequent calculations. `splitInter` has the 
# **InterProScan** results split up by the chain, and `splitRanges` has the **InterProScan** 
# ranges with the other data as *mcols*. Here they are done two different ways, with *all* of 
# the annotations, and with those that have an **InterProScan** database ID with a decent description. 

ligand_2_ipr <- subset_annotation_by_ligand(splitRanges_ipr, ligand_ranges)
ligand_ranges_full <- lapply(ligand_ranges, convert_ligands_to_ranges)
ligand_2_ipr_full <- subset_annotation_by_ligand(splitRanges_ipr, ligand_ranges_full)

chain_to_ipr <- lapply(ligand_2_ipr, get_unique_mcol)
chain_to_ipr_full <- lapply(ligand_2_ipr_full, get_unique_mcol)
site_to_chain <- split(non_red$chain.id, non_red$zinc.id)
site_to_ipr <- collapse_chain_to_annotation(chain_to_ipr, site_to_chain)


min_ligands <- c(1, 2, 3, 4)

chain_min_ligands <- lapply(min_ligands, function(min_count){
  non_red$chain.id[non_red$ligand.count >= min_count]
})
names(chain_min_ligands) <- min_ligands


annotation_subset_ligands <- lapply(chain_min_ligands, function(in_chains){
  tmp_subset <- subset_annotation_by_ligand(splitRanges_ipr, ligand_ranges[in_chains])
  tmp_annotation <- lapply(tmp_subset, get_unique_mcol)
  tmp_to_ipr <- collapse_chain_to_annotation(tmp_annotation, site_to_chain)
  annotation_2_chain <- reverseSplit(tmp_to_ipr)
})

annotation_subset_ligands_full <- lapply(chain_min_ligands, function(in_chains){
  tmp_subset <- subset_annotation_by_ligand(splitRanges_ipr, ligand_ranges_full[in_chains])
  tmp_annotation <- lapply(tmp_subset, get_unique_mcol)
  tmp_to_ipr <- collapse_chain_to_annotation(tmp_annotation, site_to_chain)
  annotation_2_chain <- reverseSplit(tmp_to_ipr)
})

setwd(paste("./simulation/", args[2], "_", args[3],sep=""))

load("normal_cluster_assg.RData")
normal_cluster_names <- paste("normal", seq(2, 30), "clusters", sep = ".")
normal_clusters <- lapply(normal_cluster_names, function(in_cluster){
  tmp <- eval(parse(text = in_cluster))
  split(tmp[,1], tmp[,2])
})
names(normal_clusters) <- normal_cluster_names

#load("compressed_cluster_assg.RData")
#compressed_cluster_names <- paste("compressed", seq(2, 30), "clusters", sep = ".")
#compressed_clusters <- lapply(compressed_cluster_names, function(in_cluster){
#  tmp <- eval(parse(text = in_cluster))
#  split(tmp[,1], tmp[,2])
#})
#names(compressed_clusters) <- compressed_cluster_names

#load("combined_cluster_assg.RData")
#combined_cluster_names <- paste("combined", seq(2, 30), "clusters", sep = ".")
#combined_clusters <- lapply(combined_cluster_names, function(in_cluster){
#  tmp <- eval(parse(text = in_cluster))
#  split(tmp[,1], tmp[,2])
#})
#names(combined_clusters) <- combined_cluster_names



normal_combs <- expand.grid(names(annotation_subset_ligands_full), names(normal_clusters), stringsAsFactors = FALSE)

get_n_clust <- function(in_name){
  as.numeric(regmatches(in_name, regexpr("\\d{1,2}", in_name)))
}

normal_combs_names <- paste(normal_combs[,1], normal_combs[,2], sep = ",")
normal_combs$n_clust <- get_n_clust(normal_combs[,2])

normal_funct_dist <- lapply(seq(1, nrow(normal_combs)), function(in_comb){
  annotation_index <- normal_combs[in_comb, 1]
  cluster_index <- normal_combs[in_comb, 2]
  calculate_functional_distance2(annotation_subset_ligands_full[[annotation_index]],
                                normal_clusters[[cluster_index]])
})

names(normal_funct_dist) <- normal_combs_names

#compressed_combs <- expand.grid(names(annotation_subset_ligands_full), names(compressed_clusters), stringsAsFactors = FALSE)
#compressed_combs_names <- paste(compressed_combs[,1], compressed_combs[,2], sep = ",")
#compressed_funct_dist <- lapply(seq(1, nrow(compressed_combs)), function(in_comb){
#  annotation_index <- compressed_combs[in_comb, 1]
#  cluster_index <- compressed_combs[in_comb, 2]
#  calculate_functional_distance2(annotation_subset_ligands_full[[annotation_index]],
#                                compressed_clusters[[cluster_index]])
#})
#names(compressed_funct_dist) <- compressed_combs_names

#compressed_combs$n_clust <- get_n_clust(compressed_combs[,2])

#combined_combs <- expand.grid(names(annotation_subset_ligands_full), names(combined_clusters), stringsAsFactors = FALSE)
#combined_combs_names <- paste(combined_combs[,1], combined_combs[,2], sep = ",")
#combined_combs$n_clust <- get_n_clust(combined_combs[,2])
#combined_funct_dist <- lapply(seq(1, nrow(combined_combs)), function(in_comb){
#  annotation_index <- combined_combs[in_comb, 1]
#  cluster_index <- combined_combs[in_comb, 2]
#  calculate_functional_distance2(annotation_subset_ligands_full[[annotation_index]],
#                                combined_clusters[[cluster_index]])
#})

normal_funct_dist <- Map(cleanup_dist, normal_funct_dist, normal_combs$n_clust)
names(normal_funct_dist) <- normal_combs_names

#compressed_funct_dist <- Map(cleanup_dist, compressed_funct_dist, compressed_combs$n_clust)
#names(compressed_funct_dist) <- compressed_combs_names

#combined_funct_dist <- Map(cleanup_dist, combined_funct_dist, combined_combs$n_clust)
#names(combined_funct_dist) <- combined_combs_names

save(normal_funct_dist, normal_combs, file = "normal_funct_dist.RData")
#save(compressed_funct_dist, compressed_combs, file = "compressed_funct_dist.RData")
#save(combined_funct_dist, combined_combs, file = "combined_funct_dist.RData")

# save(normal_clusters, compressed_clusters, file = "named_clusters.RData")
#save(normal_clusters, compressed_clusters, combined_clusters, file = "named_clusters.RData")
#save(annotation_subset_ligands_full, file = "annotation_subset_ligands_full.RData")
