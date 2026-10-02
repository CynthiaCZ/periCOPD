library(data.table)
library(igraph)
library(HMNM)
library(Matrix)
library(WGCNA)
library(dplyr)
library(stringr)

set.seed(1)

counts_raw <- fread("path/to/limma_corrected_allSamples.csv", data.table = F)
weekly_exac <- fread("path/to/periCOPD_exac_weekly.csv", data.table = F)
gene_id_to_name <- read.csv("path/to/gene_id_to_name_map.csv")

weekly_exac <- weekly_exac %>%
  mutate(week_num = as.integer(str_extract(study_week, "(?<=week_)\\d+"))) %>%
  mutate(week_num = if_else(is.na(week_num), 0L, week_num))

exac_time_df = read.table("path/to/ImpulseDE2_exac_time_label.txt", sep="\t", header = T)

top50perc_genes <- readLines("path/to/top_50perc_expressed_genes.txt")
counts_top50perc <- counts_raw %>% dplyr::select(subject_id, week_num, all_of(top50perc_genes))

# removed S1003 week 10 since missing from pheno
counts <- merge(counts_top50perc, weekly_exac %>% dplyr::select(SID, week_num, exacerbation_week_flag), by.x = c("subject_id", "week_num"), by.y = c("SID", "week_num")) %>%
  relocate(exacerbation_week_flag, .after = week_num) %>%
  arrange(subject_id, week_num)

if (!all(counts$subject_id == exac_time_df$subject_id)) {
  stop("subject ids do not match")
}

counts <- counts %>%
  mutate(exac_time=exac_time_df$time) %>%
  relocate(exac_time, .after = exacerbation_week_flag)

sum(counts$exac_time==0)

# filter for only day 0 of exac
counts_exac_day0 <- counts %>% filter(exac_time==0)
table(counts_exac_day0$subject_id)

counts_exac_day0_mat <- as.matrix(counts_exac_day0[ , grepl("^ENSG", names(counts_exac_day0))])

load("path/to/WGCNA_day0_tom_spt.RData")
load("path/to/WGCNA_day0_adjacency_mtx.RData")

TOM.dissimilarity <- 1-TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 40)
ModuleColors <- labels2colors(Modules)

gene2module <- data.frame(
  gene = colnames(counts_exac_day0_mat),
  module_color = ModuleColors
)
# replicated impulse genes
short_list <- read.csv("path/to/impulseDE_top50perc_transient_centered0_genes_norm_overlap_count.csv")
nrow(short_list)

names(all_subj_list)[names(all_subj_list) == "ensembl_id"] <- "gene_id"

cat("short in module (top 5): \n")
impulse_short_in_module <- gene2module[gene2module$gene %in% short_list$gene_id, ]
module_colors <- impulse_short_in_module$module_color
sort(table(module_colors), decreasing = TRUE)[1:5]

cat("short in merged module: \n")
impulse_short_in_merged_module <- gene2module_merged[gene2module_merged$gene %in% short_list$gene_id, ]
module_colors <- impulse_short_in_merged_module$module_color

cat("number of modules: \n")
length(unique(module_colors))
module_count <- sort(table(module_colors), decreasing = TRUE)
top_modules <- head(names(module_count), 5)
cat("Running AMEND on top modules: \n")
print(top_modules)

amend_results <- list()
for (mod in top_modules) {
  
  cat("Running AMEND for module:", mod, "\n")
  
  mod_genes <- gene2module$gene[gene2module$module_color == mod]
  
  adj_mod <- adjacency[mod_genes, mod_genes]
  
  num_nonzero_before <- sum(adj_mod > 0) / 2 # undirected
  
  adj_threshold <- 0.95
  adj_mod[adj_mod < adj_threshold] <- 0
  
  num_nonzero_after <- sum(adj_mod > 0) / 2
  
  cat("Number of edges before:", num_nonzero_before, "\n")
  cat("Number of edges after:", num_nonzero_after, "\n")

  cat("dim of adjacency matrix of module:", dim(adj_mod), "\n")
  
  g_mod <- graph_from_adjacency_matrix(adj_mod, mode = "undirected", weighted = TRUE, diag = FALSE)
  V(g_mod)$name <- mod_genes
  cat("ecount:", ecount(g_mod), "\n") 
  cat("vcount:", vcount(g_mod), "\n")
  cat("avg degree:", mean(degree(g_mod)), "\n")
  
  comp <- components(g_mod)
  cat("size of each connected component:\n", comp$csize, "\n")
  cat("fraction of nodes in the largest component:", max(comp$csize) / vcount(g_mod), "\n")
  
  bin_vec_mod <- ifelse(mod_genes %in% short_list$gene_id, 1, 0) 

  cat("overlap between model genes and impulse list:", sum(bin_vec_mod), "\n")
  names(bin_vec_mod) <- mod_genes
  
  res_mod <- HMNM::AMEND(
    network_layers = g_mod,
    n = 20,
    in_parallel=TRUE,
    data = bin_vec_mod,
    verbose = TRUE
  )
  
  amend_results[[mod]] <- res_mod

  # over write in loop in case loop brakes
  save(amend_results, file = "path/to/impulse_short_in_module_adj0.95_n20.RData")

}

print("AMEND DONE")

