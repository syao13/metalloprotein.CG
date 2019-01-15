# description -----
# does enrichment for each ligand count in normal and compressed, and compares
# them to each other. Uses categoryCompare2 as the enrichment engine, doing
# simple hypergeometric enrichment.


# Function Definitions -----

get_metal_site_list <- function(in_path){
  list_file <- file.path(in_path, "finalMetalList.RData")
  tmp_env <- new.env()
  load(list_file, envir = tmp_env)
  tmp_env$finalZnList
}

get_prot_id <- function(annotation_frame){
  prot_id <- unique(annotation_frame$X1)
  assertthat::assert_that(length(prot_id) == 1)
  prot_id
}

get_compressed_normal_data <- function(in_path){
  compressed_file <- file.path(in_path, "compressed_cluster_assg.RData")
  normal_file <- file.path(in_path, "normal_cluster_assg.RData")

  if (file.exists(compressed_file)) {
    compressed_sites <- get_cluster_members(compressed_file)
  } else {
    compressed_sites <- NA
  }

  if (file.exists(normal_file)) {
    normal_sites <- get_cluster_members(normal_file)
  } else {
    normal_sites <- NA
  }


  return(list(normal = normal_sites, compressed = compressed_sites))

}

get_go_annotation <- function(annotation_frame, ancestor_bimap){
  go_data <- mcols(annotation_frame)$X14
  go_data <- go_data[!is.na(go_data)]

  go_data <- unique(unlist(lapply(strsplit(go_data, "|", fixed = TRUE), function(x){x})))
  go_has_ancestor <- intersect(go_data, names(ancestor_bimap))
  if (length(go_has_ancestor) > 0) {
    go_ancestor <- unique(unlist(ancestor_bimap[go_has_ancestor]))
  } else {
    go_ancestor <- character(0)
  }

  all_go <- unique(c(go_data, go_ancestor))
  all_go

}

get_ipr_annotation <- function(annotation_frame){
  go_data <- mcols(annotation_frame)$X12
  go_data <- unique(go_data[!is.na(go_data)])

  go_data

}

get_cluster_members <- function(in_file){
  file_env <- new.env()
  load(in_file, envir = file_env)
  file_mem <- ls(file_env)[1]
  cluster_data <- get(file_mem, envir = file_env)
  cluster_data[,1]
}

get_go_ipr_only <- function(in_data, search_str = "GO:"){
  in_data <- grep(search_str, in_data, value = TRUE)
  table_data <- readr::read_tsv(paste(in_data, collapse = "\n"), col_names = FALSE)
  table_data
}



get_ligand_ranges <- function(in_file){
  file_env <- new.env()
  load(in_file, envir = file_env)
  file_env$ligand_ranges
}

create_labels <- function(community_number){
  all_pos_letter <- (c(LETTERS,
                      apply(combn(LETTERS, 2), 2, paste, collapse = "")))

  all_comm_number <- unique(community_number)

  out_comm_letter <- character(length(community_number))

  iletter <- 1
  for (inum in all_comm_number) {
    curr_num <- all_comm_number[inum]
    curr_pos <- which(community_number == curr_num)
    out_comm_letter[curr_pos] <- all_pos_letter[iletter]
    iletter <- iletter + 1
  }
  names(out_comm_letter) <- names(community_number)
  out_comm_letter
}

find_max_contribution <- function(annotation_lists, other_lists, use_list = "normal"){
  other_metals <- lapply(other_lists, function(x){x[[use_list]]})
  all_other <- unique(unlist(other_metals))

  max_contribution <- lapply(annotation_lists, function(in_annotation){
    in_annotation <- base::intersect(in_annotation, all_other)
    if (length(in_annotation) == 0) {
      out_data <- data.frame(metal = "NA", perc = NA)
    } else {
      out_contribution <- vapply(other_lists, function(in_list) {
        query_list <- in_list[[use_list]]
        length(base::intersect(in_annotation, query_list)) / length(in_annotation)
      }, numeric(1))

      if (all(out_contribution == 0)) {
        out_data <- data.frame(metal = "NA", perc = NA)
      } else {
        max_loc <- which.max(out_contribution)
        out_data <- data.frame(metal = names(out_contribution)[max_loc], perc = out_contribution[max_loc])
      }
    }
    out_data
  })

  do.call(rbind, max_contribution)
}

filter_annotation_features <- function(annotation_features, sites){
  sites <- unique(sites)
  annotation_universe <- unique(unlist(annotation_features))
  sites <- base::intersect(sites, annotation_universe)
  annotation_features <- lapply(annotation_features, intersect,
                              sites)
  n_feature <- sapply(annotation_features, length)
  keep_annot <- n_feature > 0
  annotation_features <- annotation_features[keep_annot]
  annotation_features
}

run_enrichment <- function(in_dir, other_lists = NULL){
  setwd(in_dir)
  # check that the necessary files exist, and stop if not
  assertthat::assert_that(file.exists("non_redundant.ipr.tsv"))
  assertthat::assert_that(file.exists("compressed_cluster_assg.RData"))
  assertthat::assert_that(file.exists("normal_cluster_assg.RData"))
  assertthat::assert_that(file.exists("non_redundant.RData"))

  ligand_range_data <- get_ligand_ranges("non_redundant.RData")


  annotation_data <- scan("non_redundant.ipr.tsv", what = character(), sep = "\n")

  # get GO annotations first
  go_functional_data <- as.data.frame(get_go_ipr_only(annotation_data, "GO:"))
  go_split <- split(go_functional_data, go_functional_data$X1)
  go_split_ranges <- lapply(go_split, atom2seq::interpro2ranges,
                            unique_on = c("X14", "X7", "X8"), require_char = "X14")

  go_split_ranges <- atom2seq::subset_annotation_by_ligand(go_split_ranges, ligand_range_data)
  go_annot <- lapply(go_split_ranges, get_go_annotation, ancestor_bimap)

  names(go_annot) <- substr(names(go_annot), start = 1, stop = nchar(names(go_annot)) - 2)
  go_2_gene <- reverseSplit(go_annot)
  go_2_gene <- lapply(go_2_gene, unique)
  go_2_gene <- go_2_gene[!(names(go_2_gene) %in% "all")]

  # remove CC because it is not that informative
  go_class <- AnnotationDbi::Ontology(names(go_2_gene))
  go_2_gene <- go_2_gene[!(go_class %in% "CC")]
  go_desc <- select(GO.db, keys = names(go_2_gene), columns = "TERM", keytype = "GOID")$TERM
  names(go_desc) <- names(go_2_gene)



  # then IPR annotations
  # ipr_functional_data <- as.data.frame(get_go_ipr_only(annotation_data, "IPR"))
  # ipr_split <- split(ipr_functional_data, ipr_functional_data$X1)
  # ipr_split_ranges <- lapply(ipr_split, atom2seq::interpro2ranges,
  #                            unique_on = c("X12", "X7", "X8"), require_char = "X12")
  #
  # ipr_annot <- lapply(ipr_split_ranges, get_ipr_annotation)
  # ipr_2_gene <- reverseSplit(ipr_annot)
  # ipr_2_gene <- lapply(ipr_2_gene, unique)
  # ipr_2_gene <- ipr_2_gene[!(names(ipr_2_gene) %in% "all")]
  #
  # ipr_table <- unique(ipr_functional_data[, c("X12", "X13")])
  # ipr_desc <- ipr_table$X13
  # names(ipr_desc) <- ipr_table$X12
  # ipr_desc <- ipr_desc[names(ipr_2_gene)]

  # functional_annot <- c(go_2_gene, ipr_2_gene)
  # functional_desc <- c(go_desc, ipr_desc)

  functional_annot <- go_2_gene
  functional_desc <- go_desc

  compress_members <- get_cluster_members("compressed_cluster_assg.RData")
  normal_members <- get_cluster_members("normal_cluster_assg.RData")
  g_universe <- unique(c(compress_members, normal_members))

  ipr_annotation <- categoryCompare2::annotation(annotation_features = functional_annot,
                                                 description = functional_desc,
                                                 type = "IPR")
  normal_enrich <- categoryCompare2::hypergeometric_feature_enrichment(
    new("hypergeom_features", significant = normal_members,
        universe = g_universe, annotation = ipr_annotation),
    p_adjust = "BH"
  )
  compress_enrich <- categoryCompare2::hypergeometric_feature_enrichment(
    new("hypergeom_features", significant = compress_members,
        universe = g_universe, annotation = ipr_annotation),
    p_adjust = "BH"
  )

  ipr_combine <- categoryCompare2::combine_enrichments(normal = normal_enrich,
                                                       compress = compress_enrich)
  ipr_sig <- categoryCompare2::get_significant_annotations(ipr_combine, p <= 0.05, counts >= 2)

  ipr_sig2 <- categoryCompare2::get_significant_annotations(ipr_combine, p <= 1)
  ipr_table2 <- as.data.frame(ipr_sig2@statistics@significant@significant) %>% dplyr::mutate(., id = rownames(.))
  ipr_table2 <- cbind(ipr_table2, as.data.frame(ipr_sig2@statistics@statistic_data[ipr_table2$id, ]))

  ipr_graph <- categoryCompare2::generate_annotation_graph(ipr_sig, low_cut = 2, hi_cut = 1000)
  ipr_graph <- categoryCompare2::remove_edges(ipr_graph, 0.8)

  ipr_sig_go <- as.data.frame(ipr_sig@statistics@significant@significant) %>% dplyr::mutate(., id = rownames(.)) %>%
    dplyr::filter(., (normal == TRUE) | (compress == TRUE))
  ipr_sig_table <- cbind(as.data.frame(ipr_sig@statistics@statistic_data[ipr_sig_go$id, ]), ipr_sig_go)
  ipr_sig_table$description <- functional_desc[ipr_sig_table$id]
  ipr_sig_table <- ipr_sig_table[order(ipr_sig_table$normal, decreasing = TRUE), ]
  ipr_sig_table <- dplyr::filter(ipr_sig_table, id %in% graph::nodes(ipr_graph))

  graph_2 <- igraph::graph_from_graphnel(ipr_graph)
  ipr_community <- igraph::membership(igraph::cluster_walktrap(graph_2))[ipr_sig_table$id]
  ipr_sig_table$IPR.group <- create_labels(ipr_community)
  ipr_sig_table <- ipr_sig_table[order(ipr_sig_table$compress, ipr_sig_table$normal, ipr_sig_table$IPR.group), ]


  if (!is.null(other_lists)) {
    compress_annot <- filter_annotation_features(functional_annot, compress_members)
    compress_max_cont <- find_max_contribution(compress_annot, other_lists, "compressed")
    names(compress_max_cont) <- paste0("compress.", names(compress_max_cont))
    compress_max_cont$compress.metal <- as.character(compress_max_cont$compress.metal)

    all_id <- unique(c(ipr_table2$id, ipr_sig_table$id))
    compress_max <- data.frame(compress.metal = rep("NA", length(all_id)),
                               compress.perc = rep(NA, length(all_id)),
                               stringsAsFactors = FALSE)
    rownames(compress_max) <- all_id
    normal_max <- compress_max
    names(normal_max) <- c("normal.metal", "normal.perc")

    normal_annot <- filter_annotation_features(functional_annot, normal_members)
    normal_max_cont <- find_max_contribution(normal_annot, other_lists, "normal")
    names(normal_max_cont) <- paste0("normal.", names(normal_max_cont))
    normal_max_cont$normal.metal <- as.character(normal_max_cont$normal.metal)

    compress_max[rownames(compress_max_cont), ] <- compress_max_cont
    normal_max[rownames(normal_max_cont), ] <- normal_max_cont
    max_contribution <- cbind(normal_max, compress_max)

    ipr_sig_table <- cbind(ipr_sig_table, max_contribution[ipr_sig_table$id, ])

    ipr_table2 <- cbind(ipr_table2, max_contribution[ipr_table2$id, ])
  }
  ipr_sig_table$type <- AnnotationDbi::Ontology(ipr_sig_table$id)
  ipr_table2$type <- AnnotationDbi::Ontology(ipr_table2$id)

  list(sig_only = ipr_sig_table, all = ipr_table2,
       members = list(compressed = compress_members,
                      normal = normal_members),
       out_file = file.path(in_dir, "normal_compress_ipr_enrich.csv"))
  #setwd("../..")
}

check_consistency <- function(table_1, table_2){
  rownames(table_1) <- table_1$id
  rownames(table_2) <- table_2$id

  if (is.null(table_1$consistent)) {
    table_1$consistent <- TRUE
  }
  if (is.null(table_2$consistent)) {
    table_2$consistent <- TRUE
  }

  can_compare <- base::intersect(table_1$id, table_2$id)

  not_consistent <- vapply(can_compare, function(in_loc){
    !identical(table_1[in_loc, c("normal", "compress")],
              table_2[in_loc, c("normal", "compress")])
  }, logical(1))

  table_1[can_compare[not_consistent], "consistent"] <- FALSE
  table_2[can_compare[not_consistent], "consistent"] <- FALSE
  rownames(table_1) <- NULL
  rownames(table_2) <- NULL
  list(table_1 = table_1, table_2 = table_2)
}

