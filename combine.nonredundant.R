#!/mlab/data/software/R-3.2.1-F22/bin/Rscript
options(stringsAsFactors=FALSE)
nrs <- NULL
nrlists <- NULL
for (metal in c("zn", "mg", "ca", "fe", "na")) {
  setwd(paste("/mlab/data/sen/projects/metal/go_ec_analysis/output_", metal, "/", sep=""))

  load("finalMetalList.RData")
  nrs <- c(nrs, finalZnList)

  non_red <- read.table("nonredundant.list.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  nrlists <- rbind(nrlists, non_red)
}

finalZnList <- nrs
save(finalZnList, file="/mlab/data/sen/projects/metal/go_ec_analysis/output_all/finalMetalList.RData")

non_red <- nrlists
write.table(non_red, file="/mlab/data/sen/projects/metal/go_ec_analysis/output_all/nonredundant.list.txt", sep="\t")

