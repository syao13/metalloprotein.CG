#!/usr/bin/Rscript

# setup -----
library(dplyr)
library(atom2seq)
library(Biobase)
library(GO.db)
ancestor_bimap <- c(as.list(GOBPANCESTOR),
                    as.list(GOMFANCESTOR))
source("normal_compressed_enrichment_functions.R")
# Run Code ----
# get the directories to run on
output_dirs <- dir("..", pattern = "output_", full.names = TRUE)
output_dirs <- grep("_out", output_dirs, invert = TRUE, value = TRUE)

# get the compressed and normal members
compressed_normal_lists <- lapply(output_dirs, function(in_metal){
  ligand_dirs <- c(dir(in_metal, pattern = "ligand"), "allLig")

  ligand_path <- file.path(in_metal, ligand_dirs)

  compress_normal <- lapply(ligand_path, get_compressed_normal_data)
  names(compress_normal) <- ligand_dirs
  compress_normal
})

# check the numbers. We expect that
# * allLigand should equal sum of others in a metal
# * all 4ligand should equal sum of ca 4ligand, fe 4ligand, etc
# * all allLigand should be sum of others, or sum of metal allLigand's
# These should all work **within** the normal and the compressed

## Check each metal
within_metal_sums <- lapply(compressed_normal_lists, function(in_metal){
  n_compressed_allLigand <- length(unique(in_metal$allLig$compressed))

  n_compressed_other <- sum(vapply(in_metal[!(names(in_metal) %in% "allLig")],
                                   function(in_group){
                                     length(unique(in_group$compressed))
                                   }, numeric(1)))

  n_normal_allLigand <- length(unique(in_metal$allLig$normal))

  n_normal_other <- sum(vapply(in_metal[!(names(in_metal) %in% "allLig")],
                                   function(in_group){
                                     length(unique(in_group$normal))
                                   }, numeric(1)))
  data.frame(normal.all = n_normal_allLigand,
             normal.other = n_normal_other,
             compressed.all = n_compressed_allLigand,
             compressed.other = n_compressed_other)

})

within_metal_sums <- do.call(rbind, within_metal_sums)
within_metal_sums$metal <- output_dirs

indiv_metal_all <- within_metal_sums[2:6, 1:4]
sum_indiv <- as.data.frame(t(colSums(indiv_metal_all)))
sum_indiv$metal <- "summing indiv metals"

metal_sums <- rbind(within_metal_sums, sum_indiv)

# ok, some really, really weird things are going on here. Lets start looking
# at intersects compared to sizes and see where we end up.


# basic enrichment -----
# get first round of tables on every metal ligand combination
results_tables <- lapply(output_dirs, function(in_metal){
  use_dir <- getwd()
  metal <- strsplit(in_metal, "_", fixed = TRUE)[[1]][2]
  #metal <- substr(in_metal, nchar(in_metal) - 1, nchar(in_metal))

  ligand_dirs <- c(dir(in_metal, pattern = "ligand"), "allLig")

  n_lig <- substr(ligand_dirs, 1, 1)

  metal_lig <- paste0(metal, "_", n_lig)

  ligand_path <- file.path(in_metal, ligand_dirs)


  ligand_enrich <- lapply(ligand_path, function(in_path){
    print(in_path)
    out_enrich <- try(run_enrichment(in_path))
    setwd(use_dir)
    out_enrich
  })
  names(ligand_enrich) <- metal_lig
  ligand_enrich

})

# which metal ----
# rename them by metal so they are easy to know which is which
metal <- vapply(strsplit(output_dirs, "_", fixed = TRUE), function(x){x[2]}, character(1))
names(results_tables) <- metal

# get compressed and normal from each metal so we can see what is responsible
# for enrichment later
compressed_normal_members <- lapply(results_tables[c("ca", "fe", "mg", "na", "zn")], function(x){
  is_a <- grep("_a", names(x))
  x[[is_a]][["members"]]
})

# run code for getting which is responsible for enrichment
all_loc <- "output_all"
results_all2 <- lapply(all_loc, function(in_metal){
  use_dir <- getwd()
  metal <- strsplit(in_metal, "_", fixed = TRUE)[[1]][2]
  #metal <- substr(in_metal, nchar(in_metal) - 1, nchar(in_metal))

  ligand_dirs <- c(dir(in_metal, pattern = "ligand"), "allLig")

  n_lig <- substr(ligand_dirs, 1, 1)

  metal_lig <- paste0(metal, "_", n_lig)

  ligand_path <- file.path(in_metal, ligand_dirs)

  ligand_enrich <- lapply(ligand_path, function(in_path){
    print(in_path)
    out_enrich <- try(run_enrichment(in_path, compressed_normal_members))
    setwd(use_dir)
    out_enrich
  })
  names(ligand_enrich) <- metal_lig
  ligand_enrich
})

# remove anything that didn't actually have results
for (imetal in names(results_tables)) {
  for (ilig in names(results_tables[[imetal]])) {
    if (!is.list(results_tables[[imetal]][[ilig]])) {
      results_tables[[imetal]][[ilig]] <- NULL
    }
  }
}

# consistency -----

# check consistency on "all_a"
results_all2 <- results_all2[[1]]

all2_all <- results_all2$all_a$sig_only
for (icomp in c("all_4", "all_5", "all_6", "all_7", "all_8")) {
  table_2 <- results_all2[[icomp]]$sig_only
  out_data <- check_consistency(all2_all, table_2)
  all2_all <- out_data$table_1
  results_all2[[icomp]]$sig_only <- out_data$table_2
}

# check against the individual metals
for (imetal in c("ca", "fe", "mg", "na", "zn")) {
  a_metal <- paste0(imetal, "_a")
  table_2 <- results_tables[[imetal]][[a_metal]]$sig_only
  out_data <- check_consistency(all2_all, table_2)
  all2_all <- out_data$table_1
  results_tables[[imetal]][[a_metal]]$sig_only <- out_data$table_2
}

results_all2$all_a$sig_only <- all2_all

# check for each metal, the *all* against the *sub-groups*
for (imetal in c("ca", "fe", "mg", "na", "zn")) {
  a_metal <- paste0(imetal, "_a")
  metal_ligs <- names(results_tables[[imetal]])
  metal_ligs <- metal_ligs[!(metal_ligs %in% a_metal)]
  table_1 <- results_tables[[imetal]][[a_metal]]$sig_only
  for (ilig in metal_ligs) {
    if (is.list(results_tables[[imetal]][[ilig]])) {
      table_2 <- results_tables[[imetal]][[ilig]]$sig_only

      out_data <- check_consistency(table_1, table_2)
      results_tables[[imetal]][[ilig]]$sig_only <- out_data$table_2
      table_1 <- out_data$table_1
    }

  }
  results_tables[[imetal]][[a_metal]]$sig_only <- table_1
}

results_tables$all <- results_all2

# overlap -----
check_overlap_normal_compressed <- lapply(results_tables, function(x){
  lapply(x, function(y){
    z <- y[[1]]
    if (class(z) == "data.frame") {
      z <- group_by(z, IPR.group)
      group_summary <- summarize(z, n_normal = sum(normal), n_compress = sum(compress))
      out_val <- nrow(filter(group_summary, (n_normal != 0) & (n_compress != 0)))
    } else {
      out_val <- 0
    }
    out_val

  })
})

# write files ------

# we write the files to two locations, the actual output directories, and then
# to a folder in the main directory
save_loc <- "all_enrichments"
lapply(results_tables, function(in_res){
  for (ires in names(in_res)) {
    if (is.list(in_res[[ires]])) {
      write.table(in_res[[ires]]$sig_only, file = in_res[[ires]]$out_file, sep = ",", row.names = FALSE,
                  col.names = TRUE)
      other_file <- file.path(save_loc, paste0(ires, "_enrichment.csv"))
      write.table(in_res[[ires]]$sig_only, file = other_file, sep = ",", row.names = FALSE,
                  col.names = TRUE)
    }

  }
})


save(list = ls(), file = "all_results_data.RData")
# render results -----
#rmarkdown::render("normal_compressed_enrichments_pdf.Rmd", output_file = "normal_compressed_enrichment_tables.pdf")
rmarkdown::render("normal_compressed_enrichments_padj05.Rmd", "Table_S149_S168_five_metal_sppl_tables.pdf")
#rmarkdown::render("normal_compressed_enrichments_manuscriptTable.Rmd")
rmarkdown::render("examine_differences.Rmd", "Figure_S29_five_metal_sppl_figures.pdf")
#rmarkdown::render("")


