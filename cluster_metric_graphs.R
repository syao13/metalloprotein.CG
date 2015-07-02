#!/usr/bin/Rscript

# creating multiple k-graphs

# There is a graph we need examining cluster distances, jaccard similarity, structural vs functional distances rho 
# and p-values. We want this in a single graphic that makes it easy to compare across num_ligands and metric
library(ggplot2)
library(Cairo)
setwd("../firstManuscriptOutput")
load("four_measures_over_k.RData")

level_order <- c("jaccard", "max - sumdiff", "rho", "-1 * log10(p)")

# normal clusters

normal_jaccard <- data.frame(cluster = seq(1, 30), value = jaccard.norm, type = "jaccard", n_lig = "all")
normal_sumdiff <- data.frame(cluster = seq(1, 30), value = max(sumdiff.norm) - sumdiff.norm, type = "max - sumdiff", n_lig = "all")

normal_rho <- data.frame(cluster = rep(seq(3, 30), 4), 
                         value = c(rhop.1.norm$rho,
                                   rhop.2.norm$rho,
                                   rhop.3.norm$rho,
                                   rhop.4.norm$rho),
                         type = "rho",
                         n_lig = rep(c("1", "2", "3", "4"), each = 28))
normal_p <- data.frame(cluster = rep(seq(3, 30), 4), 
                       value = c(rhop.1.norm$p_value,
                                 rhop.2.norm$p_value,
                                 rhop.3.norm$p_value,
                                 rhop.4.norm$p_value),
                       type = "-1 * log10(p)",
                       n_lig = rep(c("1", "2", "3", "4"), each = 28))

min_p <- min(normal_p$value[normal_p$value != 0])
normal_p$value[normal_p$value <= 1e-10] <- 1e-10
normal_p$value <- -1 * log10(normal_p$value)

normal_metrics <- do.call(rbind, list(normal_jaccard, normal_sumdiff, normal_rho, normal_p))
normal_metrics$type <- ordered(normal_metrics$type, levels = level_order)

Cairo(file = "Figure7.png", type = "png")
ggplot(normal_metrics, aes(x = cluster, y = value, color = n_lig)) + geom_point() + geom_line() + 
  facet_grid(type ~ ., scales = "free_y") + ggtitle("normal")
dev.off()

### compressed clusters
### 
compressed_jaccard <- data.frame(cluster = seq(1, 30), value = jaccard.comp, type = "jaccard", n_lig = "all")
compressed_sumdiff <- data.frame(cluster = seq(1, 30), value = max(sumdiff.comp) - sumdiff.comp, type = "max - sumdiff", n_lig = "all")

compressed_rho <- data.frame(cluster = rep(seq(3, 30), 4),
                             value = c(rhop.1.comp$rho,
                                       rhop.2.comp$rho,
                                       rhop.3.comp$rho,
                                       rhop.4.comp$rho),
                             type = "rho",
                             n_lig = rep(c("1", "2", "3", "4"), each = 28))

compressed_p <- data.frame(cluster = rep(seq(3, 30), 4),
                           value = c(rhop.1.comp$p_value,
                                     rhop.2.comp$p_value,
                                     rhop.3.comp$p_value,
                                     rhop.4.comp$p_value),
                           type = "-1 * log10(p)",
                           n_lig = rep(c("1", "2", "3", "4"), each = 28))
compressed_p$value <- -1 * log10(compressed_p$value)

compressed_metrics <- do.call(rbind, list(compressed_jaccard, compressed_sumdiff, compressed_rho, compressed_p))
compressed_metrics$type <- ordered(compressed_metrics$type, levels = level_order)

Cairo(file = "Figure8.png", type = "png")
ggplot(compressed_metrics, aes(x = cluster, y = value, color = n_lig)) + geom_point() + geom_line() + 
  facet_grid(type ~ ., scales = "free_y") + ggtitle("Compressed")
dev.off()

### combined cluster measures
all_jaccard <- data.frame(cluster = seq(1, 30), value = jaccard.all, type = "jaccard", n_lig = "all")
all_sumdiff <- data.frame(cluster = seq(1, 30), value = max(sumdiff.all) - sumdiff.all, type = "max - sumdiff", n_lig = "all")

all_rho <- data.frame(cluster = rep(seq(3, 30), 4), 
                         value = c(rhop.1.all$rho,
                                   rhop.2.all$rho,
                                   rhop.3.all$rho,
                                   rhop.4.all$rho),
                         type = "rho",
                         n_lig = rep(c("1", "2", "3", "4"), each = 28))
all_p <- data.frame(cluster = rep(seq(3, 30), 4), 
                       value = c(rhop.1.all$p_value,
                                 rhop.2.all$p_value,
                                 rhop.3.all$p_value,
                                 rhop.4.all$p_value),
                       type = "-1 * log10(p)",
                       n_lig = rep(c("1", "2", "3", "4"), each = 28))

min_p <- min(all_p$value[all_p$value != 0])
all_p$value[all_p$value <= 1e-10] <- 1e-10
all_p$value <- -1 * log10(all_p$value)

all_metrics <- do.call(rbind, list(all_jaccard, all_sumdiff, all_rho, all_p))
all_metrics$type <- ordered(all_metrics$type, levels = level_order)

Cairo(file = "FigureS3.png", type = "png")
ggplot(all_metrics, aes(x = cluster, y = value, color = n_lig)) + geom_point() + geom_line() + 
  facet_grid(type ~ ., scales = "free_y") + ggtitle("comined")
dev.off()
