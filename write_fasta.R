#!/usr/bin/Rscript
############################################################################################
##
##   Written by Sen Yao, 07/20/2016
##   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/20/2016. All rights reserved.
##
###########################################################################################

library(Biostrings)
library(magrittr)
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
non_red <- read.table("nonredundant.list.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
non_red$chain.id <- paste(non_red$zinc.id, non_red$chain.counts, sep = ":")
single_chains <- !duplicated(non_red$chain.id)
non_red <- non_red[single_chains,]
non_red_seq <- AAStringSet(non_red$sequence)
names(non_red_seq) <- non_red$chain.id

ligand_locs <- strsplit(non_red$ligand.position.combo, ".", fixed = TRUE) %>% lapply(., as.integer)
upper_code <- toupper(AMINO_ACID_CODE)
single_code <- names(upper_code)
ligand_combo <- strsplit(non_red$ligand.combo, ".", fixed = TRUE) %>% lapply(., function(x){paste(single_code[match(x, upper_code)], sep="", collapse="")})
seq_combo <- lapply(seq(1, nrow(non_red)), function(x){as.character(non_red_seq[[x]][ligand_locs[[x]]])})

# sanity check that it all looks right
all.equal(ligand_combo, seq_combo)
non_red$ligand.count <- sapply(ligand_locs, length)

# For each chain, we will record the location of the ligands as positions in an `IRanges` instance 
# that can be queried for overlap with the results from `Interproscan`.
ligand_ranges <- lapply(ligand_locs, function(x){IRanges(x, width = 1)})
names(ligand_ranges) <- non_red$chain.id

save(non_red, ligand_ranges, file = "non_redundant.RData")
writeXStringSet(non_red_seq, file = "non_redundant.fa")
