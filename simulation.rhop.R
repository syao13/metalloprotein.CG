#!/usr/bin/Rscript
############################################################################################
##
##   Written by Sen Yao, 07/20/2016
##   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/20/2016. All rights reserved.
##
###########################################################################################

###
## args: simulation sample number, normal/compressed/combined, cluster number k

args = commandArgs(trailingOnly=TRUE)

rhos <- NULL
ps <- NULL

for (i in 1:20) {
  load(paste("./", args[1], "_", i, "/", "four_measures_over_k.RData", sep=""))

  group <- NULL
  if (args[2] == "normal") { 
    group <- "norm"
  } else if (args[2] == "compressed") {
    group <- "comp"
  } else if (args[2] == "combined") {
    group <- "all"
  }

  rhop <- get(paste("rhop.", 4, ".", group, sep=""))
  rhos <- c(rhos, round(rhop$rho[as.numeric(args[3])-2],4))
  ps <- c(ps, round(rhop$p_value[as.numeric(args[3])-2],4))
  }

sum(rhos)/20
rhos
ps

