#!/usr/bin/Rscript

## args: metal, numLig, normal/compressed/combined, k, n_lig
## i.e.: Zn, all, normal, 7, 3
args = commandArgs(trailingOnly=TRUE)

#if (args[1] == "all") {args[1] <- "allMetal"}

if (args[2] == "all" || args[2] == "combined") { setwd(paste("../output_", tolower(args[1]), "/allLig", sep=""))
} else if (args[2] == "noheme") { setwd(paste("../output_", tolower(args[1]), "/noheme", sep=""))
} else { setwd(paste("../output_", tolower(args[1]), "/", args[2], "ligand", sep=""))}

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
print(c(round(rhop$rho[as.numeric(args[4])-2],4), round(rhop$p_value[as.numeric(args[4])-2],4)))

