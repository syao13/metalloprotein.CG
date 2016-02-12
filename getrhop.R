#!/usr/bin/Rscript

## args: metal, numLig, normal/compressed/combined, n_lig, k

args = commandArgs(trailingOnly=TRUE)

if (args[1] == "all") {args[1] <- "allMetal"}

if (args[2] == "all") { setwd(paste("../", args[1], "_results/output_", args[1], "_single6", sep=""))} 
if (args[2] != "all") { setwd(paste("../", args[1], "_results/", args[2], "ligand_6/", sep=""))}

load("four_measures_over_k.RData")

group <- NULL
if (args[3] == "normal") { 
  group <- "norm"
} else if (args[3] == "compressed") {
  group <- "comp"
} else if (args[3] == "combined") {
  group <- "all"
}

rhop <- get(paste("rhop.", args[5], ".", group, sep=""))
round(rhop$rho[as.numeric(args[4])-2],4)
round(rhop$p_value[as.numeric(args[4])-2],4)

