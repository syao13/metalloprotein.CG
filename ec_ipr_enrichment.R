#!/usr/bin/Rscript
library(rex)
library(categoryCompareANY)
setwd("../output")

#### EC enrichment
test_ec1 <- c("KEGG: 00910+4.2.1.1|MetaCyc: PWY-6142", "KEGG: 00310+3.4.-.-|",
              "KEGG: 00627+3.1.3.1|", "KEGG: 00627+3.-.-.-")

ec_rex <- rex(
  "KEGG:", space, n_times(numbers, 5), "+",
  capture(numbers, ".", numbers %or% "-", ".", numbers %or% "-", ".", numbers %or% "-"),
  maybe(anything))

re_matches(test_ec1, ec_rex)

# function to do the actual parsing of EC numbers
ec_parse <- function(ec_regex, id, pathway_result){
  keep_paths <- nchar(pathway_result) > 0
  id <- id[keep_paths]
  pathway_result <- pathway_result[keep_paths]
  
  split_paths <- strsplit(pathway_result, "|", fixed = TRUE)
  
  ec_result <- lapply(split_paths, function(in_path){
    ec_match <- re_matches(in_path, ec_regex)
    ec_match <- unique(ec_match)
    ec_match <- ec_match[!is.na(ec_match)]
    
    ec_match
  })
  ec_length <- sapply(ec_result, length)
  id <- id[ec_length > 0]
  ec_result <- ec_result[ec_length > 0]
  ec_length <- ec_length[ec_length > 0]
  
  ec_out <- data.frame(id = rep(id, ec_length), ec = unlist(ec_result), stringsAsFactors = FALSE)
  unique(ec_out)
}

collapse_ec <- function(ec_matrix, ec_annotation, ec_loc){
  out_2_loc <- apply(ec_matrix, 1, function(x){paste(x[1:ec_loc],  collapse = ".")})
  split_on_loc <- split(rownames(ec_matrix), out_2_loc)
  collapsed_on_loc <- lapply(split_on_loc, function(x){
    unique(unlist(ec_annotation[x], use.names = FALSE))
  })
  collapsed_on_loc
}


# Load interproscan results
interproRes <- read.table("non_redundant_zinc.ipr.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE, quote="", comment.char="")
names(interproRes) <- c("id", "md5", "length", "analysis", "sigAccession", "sigDescription", "start", "stop", "score", "status", "date", "iprAccession", "iprDescription", "GO", "Pathways")

iprData <- unique(interproRes[, c("iprAccession", "iprDescription")])
rownames(iprData) <- iprData$iprAccession


load("non_redundant.RData")
load("named_clusters.RData")

# get them out
inter_ec <- ec_parse(ec_rex, interproRes$id, interproRes$Pathways)
rownames(non_red) <- non_red$chain.id
inter_ec$zinc.id <- non_red[inter_ec$id, "zinc.id"]

# map between value and the 4 digits
unique_ec <- unique(inter_ec$ec)
split_ec <- strsplit(unique_ec, ".", fixed = TRUE)
num_ec <- lapply(split_ec, function(in_ec){
  tmp_ec <- as.numeric(in_ec)
  tmp_ec[is.na(tmp_ec)] <- 0
  tmp_ec
})
num_ec <- do.call(rbind, num_ec)
rownames(num_ec) <- unique_ec
num_ec <- num_ec[order(num_ec[,1], num_ec[,2], num_ec[,3], num_ec[,4]),]

ec_2_zincid <- split(inter_ec$zinc.id, inter_ec$ec)

has_ec <- unique(inter_ec$zinc.id)
no_ec <- unique(setdiff(non_red$zinc.id, has_ec))
ec_2_zincid$'0.0.0.0' <- no_ec
num_ec <- rbind(num_ec, "0.0.0.0" = c(0, 0, 0, 0))

# do the enrichments
ec_2 <- collapse_ec(num_ec, ec_2_zincid, 2)

ec_2_normal <- new("ANYHyperGParamsCC",
                   geneIds = unlist(normal_clusters$normal.30.clusters),
                   universeGeneIds = unique(non_red$zinc.id),
                   categoryName = "ANY",
                   annotation = ec_2,
                   fdr = 0)

ec_2_normal_res <- hyperGTestCC(ec_2_normal)

ec_2_compress <- new("ANYHyperGParamsCC",
                     geneIds = unlist(compressed_clusters$compressed.30.clusters),
                     universeGeneIds = unique(non_red$zinc.id),
                     categoryName = "ANY",
                     annotation = ec_2,
                     fdr = 0)

ec_2_compress_res <- hyperGTestCC(ec_2_compress)

comp_res <- data.frame(compress = ec_2_compress_res@pvalues)
norm_res <- data.frame(normal = ec_2_normal_res@pvalues)

all_ec <- unique(c(rownames(comp_res), rownames(norm_res)))
all_res <- data.frame(compress = rep(1, length(all_ec)), normal = rep(1, length(all_ec)))
rownames(all_res) <- all_ec

all_res[rownames(comp_res), "compress"] <- comp_res
all_res[rownames(norm_res), "normal"] <- norm_res

all_res$comp_norm <- (-1 * log10(all_res$compress)) - (-1 * log10(all_res$normal))
all_res <- all_res[(order(all_res$comp_norm, decreasing = TRUE)),]
save(all_res, file = "ec2_pvalues.RData")
write.table(all_res, file = "ec2_pvalues.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


#### IPR Enrichment
load("annotation_subset_ligands_full.RData")
ipr_annot <- annotation_subset_ligands_full[[1]]

ipr_normal <- new("ANYHyperGParamsCC",
                  geneIds = unlist(normal_clusters$normal.30.clusters),
                  universeGeneIds = unique(non_red$zinc.id),
                  categoryName = "ANY",
                  annotation = ipr_annot,
                  fdr = 0)

ipr_normal_res <- hyperGTestCC(ipr_normal)
pvalueCutoff(ipr_normal_res) <- 1
ipr_normal_table <- summary(ipr_normal_res)

ipr_compress <- new("ANYHyperGParamsCC",
                    geneIds = unlist(compressed_clusters$compressed.30.clusters),
                    universeGeneIds = unique(non_red$zinc.id),
                    categoryName = "ANY",
                    annotation = ipr_annot,
                    fdr = 0)
ipr_compress_res <- hyperGTestCC(ipr_compress)
pvalueCutoff(ipr_compress_res) <- 1
ipr_compress_table <- summary(ipr_compress_res)

full_ipr <- merge(ipr_normal_table, ipr_compress_table, by.x = "ID", by.y = "ID", all = TRUE, suffixes = c(".normal", ".compressed"))

full_ipr <- full_ipr[order(full_ipr$Pvalue.normal, full_ipr$Pvalue.compressed, decreasing = FALSE), ]
full_ipr$description <- iprData[full_ipr$ID, "iprDescription"]

ipr_tmp_normal <- full_ipr$Pvalue.normal
ipr_tmp_normal[is.na(ipr_tmp_normal)] <- 1
ipr_tmp_compressed <- full_ipr$Pvalue.compressed
ipr_tmp_compressed[is.na(ipr_tmp_compressed)] <- 1

ipr_comp_normal <- (-1 * log10(ipr_tmp_compressed)) - (-1 * log10(ipr_tmp_normal))
full_ipr$compressed_to_normal <- ipr_comp_normal

full_ipr <- full_ipr[, c("ID", "description", "Pvalue.normal", "Pvalue.compressed", "compressed_to_normal", "FDR.normal",
                         "OddsRatio.normal", "ExpCount.normal",
                         "Count.normal", "Size.normal", "FDR.compressed", "OddsRatio.compressed",
                         "ExpCount.compressed", "Count.compressed", "Size.compressed")]
full_ipr <- full_ipr[order(full_ipr$compressed_to_normal),]
write.table(full_ipr, file = "ipr_enrichment.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


pcut <- 0.05

ipr_spec_compress <- subset(full_ipr, (is.na(Pvalue.normal) & (Pvalue.compressed <= pcut)))
ipr_spec_normal <- subset(full_ipr, (is.na(Pvalue.compressed) & (Pvalue.normal <= pcut)))

write.table(ipr_spec_compress, file = "ipr_compressed_only.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(ipr_spec_normal, file = "ipr_normal_only.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
