#!/mlab/data/software/R-3.2.1-F22/bin/Rscript
options(stringsAsFactors=FALSE)
predictions <- NULL
for (metal in c("zn", "ca", "na", "mg", "fe")) {
  setwd(paste("/mlab/data/sen/projects/metal/go_ec_analysis/output_", metal, "/", sep=""))

  load("rf.results.RData")
  predictions <- c(predictions, prediction.all)
}

prediction.all <- predictions
save(prediction.all, file="/mlab/data/sen/projects/metal/go_ec_analysis/output_all/rf.results.RData")

