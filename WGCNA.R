library(data.table)
library(ImpulseDE2)
library(WGCNA)
library(dplyr)
library(stringr)

# load data
counts_raw <- fread("path/to/limma_corrected_allSamples.csv", data.table = F)
weekly_exac <- fread("path/to/periCOPD_exac_weekly.csv", data.table = F)
gene_id_to_name <- read.csv("path/to/gene_id_to_name_map.csv")

weekly_exac <- weekly_exac %>%
  mutate(week_num = as.integer(str_extract(study_week, "(?<=week_)\\d+"))) %>%
  mutate(week_num = if_else(is.na(week_num), 0L, week_num))


exac_time_df = read.table("path/to/ImpulseDE2_exac_time_label.txt", sep="\t")

top50perc_genes<- readLines("path/to/top_50perc_expressed_genes.txt")

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

spt <- pickSoftThreshold(counts_exac_day0_mat) 

softPower <- 7
adjacency <- adjacency(counts_exac_day0_mat, power = softPower)

save(adjacency, file = "path/to/WGCNA_day0_adjacency_mtx.RData")
cat("saved adjacency matrix \n")

TOM <- TOMsimilarity(adjacency)

save(TOM, spt, file = "path/to/WGCNA_day0_tom_spt.RData")
cat("saved TOM and spt")
