library(data.table)
library(dplyr)
library(tableone)
library(MASS)
library(ggplot2)
library(glmnet)
library(caret)


call_model <- function(df, outcome_name, top_genes, covars, adjustment=NULL, model_type) {
  results <- data.frame(gene = character(), beta = numeric(), ci_lower = numeric(), ci_upper = numeric(), beta_ci = character(), pval = numeric(), stringsAsFactors = FALSE)
  
  for (gene in top_genes) {
    gene <- ifelse(grepl("\\|", gene), paste0("`", gene, "`"), gene)
    
    formula_components <- c(gene, covars, adjustment)
    formula_str <- paste0(outcome_name, " ~ ", paste0(formula_components, collapse = " + "))
    
    tryCatch({
      if (model_type == "nb") {
        model <- glm.nb(as.formula(formula_str), data = df)
      } else if (model_type == "binomial") {
        model <- glm(as.formula(formula_str), data = df, family = binomial())
      }
      
      beta <- coef(model)[gene]
      ci <- confint(model, gene, level = 0.95)
      beta_ci <- sprintf("%.3f (%.3f, %.3f)", beta, ci[1], ci[2])
      pval <- summary(model)$coefficients[gene, 4]
      results <- rbind(results, data.frame(gene = gene, beta = beta, ci_lower = ci[1], ci_upper = ci[2], beta_ci = beta_ci, pval = pval))
    }, error = function(e) {
      message(sprintf("Skipping gene %s due to error: %s", gene, e$message))
    })
  }
  return(results)
}

rnaseq_all <- fread("path/to/counts_data.tsv", data.table = F)
rnaseq_all$GENEID <- sub("\\..*", "", rnaseq_all$GENEID)

counts_raw <- rnaseq_all[, -1]
rownames(counts_raw) <- rnaseq_all[, 1]
gene_means <- rowMeans(counts_raw, na.rm = TRUE)

threshold <- quantile(gene_means, 0.5)
top_genes <- names(gene_means[gene_means > threshold])
rnaseq_all_top <- rnaseq_all[rnaseq_all[, 1] %in% top_genes, ]


GENEID_all <- rnaseq_all_top[, 1]
rnaseq_mtx_all <- rnaseq_all_top[, -1]
rnaseq_transpose_all <- as.data.frame(t(rnaseq_mtx_all))
colnames(rnaseq_transpose_all) <- GENEID_all

rnaseq_transpose_all$sid <- rownames(rnaseq_transpose_all)
rownames(rnaseq_transpose_all) <- NULL

rnaseq_transpose_all <- rnaseq_transpose_all %>%
  mutate(across(-sid, ~ log2(.x + 1)))


load("path/to/pheno_data.Rdata")

lfu <- fread("path/to/followup_data.txt", data.table = F)
lfu_exac <- lfu %>% 
  filter(!is.na(days_since_P2)) %>%
  group_by(sid) %>%
  summarise(total_days_since_P2 = max(days_since_P2, na.rm = TRUE))

pheno_p2 <- pheno_p2 %>% dplyr::select(-total_days_since_P2)
pheno_p2 <- merge(pheno_p2, lfu_exac, by='sid') 


covars <- c("ATS_PackYears_P2", "Exacerbation_Frequency_P2")

sapply(pheno_p2[c(covars, "total_days_since_P2")], function(x) sum(is.na(x)))
pheno_p2 <- pheno_p2[complete.cases(pheno_p2[, c(covars, "total_days_since_P2")]), ]

pheno_p2$copd <- ifelse(pheno_p2$FEV1_FVC_post_P2 < 0.7 & pheno_p2$fev1_pp_gli_global_P2 < 80, 1,
                        ifelse(pheno_p2$FEV1_FVC_post_P2 >= 0.7 & pheno_p2$fev1_pp_gli_global_P2 >= 80, 0, NA))

pheno_copd <- pheno_p2 %>% filter(copd == 1)
pheno_copd$exacerbation_rate = pheno_copd$total_Net_Exacerbations / pheno_copd$total_days_since_P2 * 365

# fill pack years with 0 if smoking is 0
pheno_copd$ATS_PackYears_P2[is.na(pheno_copd$ATS_PackYears_P2) & pheno_copd$SmokCigNow_P2 == 0] <- 0



df_all_rnaseq <- merge(rnaseq_transpose_all, pheno_copd, on="sid")
# binary exacerbation
df_all_rnaseq <- df_all_rnaseq %>% mutate(has_exacerbation = if_else(total_Net_Exacerbations>=1, 1, 0))

df_all_rnaseq$has_exacerbation <- as.factor(df_all_rnaseq$has_exacerbation)
df_all_rnaseq$gender <- as.factor(df_all_rnaseq$gender)
df_all_rnaseq$race <- as.factor(df_all_rnaseq$race)

table(df_all_rnaseq$has_exacerbation)

set.seed(138)
train_indices <- createDataPartition(df_all_rnaseq$has_exacerbation, p = 0.5, list = FALSE)
df_train_all_rnaseq <- df_all_rnaseq[train_indices, ]
df_test_all_rnaseq  <- df_all_rnaseq[-train_indices, ]

regress_vars_age_sex <- c("Age_P2","gender")

gene_resid_models <- list()
for(gene in top_genes){

  gene_resid_models[[gene]] <- lm(as.formula(paste0(gene,"~",paste0(regress_vars_age_sex,collapse = "+"))),data=df_train_all_rnaseq)
  # apply train model to both train and test
  df_train_all_rnaseq[, gene] <- as.numeric(resid(gene_resid_models[[gene]]))
  df_test_all_rnaseq[, gene]  <- df_test_all_rnaseq[, gene] -
    as.numeric(predict(gene_resid_models[[gene]], newdata = df_test_all_rnaseq))
}

X_train_all_rnaseq <- as.matrix(df_train_all_rnaseq[, top_genes])
y_train_all_rnaseq <- df_train_all_rnaseq$has_exacerbation

pos_wt <- mean(df_train_all_rnaseq$has_exacerbation == 0)
neg_wt <- mean(df_train_all_rnaseq$has_exacerbation == 1)
obs_wt <- ifelse(df_train_all_rnaseq$has_exacerbation == 1, pos_wt, neg_wt)

X_covars <- model.matrix(~ ATS_PackYears_P2 + Exacerbation_Frequency_P2, data = df_train_all_rnaseq)[, -1]

X_train_all_rnaseq_w_covars <- cbind(X_train_all_rnaseq, X_covars)
penalty_vec_ridge <- c(rep(1, ncol(X_train_all_rnaseq)), rep(0, ncol(X_covars)))

ridge_all_rnaseq <- cv.glmnet(
  x = X_train_all_rnaseq_w_covars,
  y = y_train_all_rnaseq,
  type.measure = "deviance",
  family = "binomial",
  nfolds = 5,
  alpha = 0,
  weights = obs_wt,
  trace.it=1,
  penalty.factor = penalty_vec_ridge,
  standardize = TRUE
)

ridge_coef_all_rnaseq <- as.numeric(coef(ridge_all_rnaseq, s = "lambda.min"))[-1]
beta_genes <- ridge_coef_all_rnaseq[1:ncol(X_train_all_rnaseq)]
w_genes <- 1 / abs(beta_genes)
penalty_vec_adalasso <- c(w_genes, rep(0, ncol(X_covars)))

adalasso_all_rnaseq <- cv.glmnet(
  x = X_train_all_rnaseq_w_covars,
  y = y_train_all_rnaseq,
  type.measure = "auc",
  family = "binomial",
  nfolds = 5,
  alpha = 1,
  weights = obs_wt,
  trace.it=1,
  standardize = TRUE,
  penalty.factor = penalty_vec_adalasso
)

coef_ada <- coef(adalasso_all_rnaseq, s = "lambda.min")
sum(coef_ada != 0 & !is.na(coef_ada))

coef_ada_dense <- as.matrix(coef_ada)

selected_genes <- rownames(coef_ada_dense)[which(coef_ada_dense != 0 & !is.na(coef_ada_dense))]
selected_genes <- setdiff(selected_genes, "(Intercept)")
selected_genes <- intersect(selected_genes, colnames(X_train_all_rnaseq))


beta_selected_all_rnaseq_adalasso <- coef_ada_dense[selected_genes, 1]
names(beta_selected_all_rnaseq_adalasso) <- selected_genes

beta_selected_all_rnaseq_adalasso

resid_coefs_all <- t(sapply(selected_genes, function(gene) {
  coef(gene_resid_models[[gene]])
}))

age_weight_all <- sum(-beta_selected_all_rnaseq_adalasso * resid_coefs_all[, "Age_P2"])
sex_weight_all <- sum(-beta_selected_all_rnaseq_adalasso * resid_coefs_all[, "gender2"])
intercept_adj_all <- sum(-beta_selected_all_rnaseq_adalasso * resid_coefs_all[, "(Intercept)"])

formula_weights <- c(
  beta_selected_all_rnaseq_adalasso,
  age = age_weight_all,
  sex = sex_weight_all,
  intercept = intercept_adj_all
)

formula_weights_df <- data.frame(Term = names(formula_weights), Coefficient = formula_weights)

write.table(formula_weights_df, "path/to/cgRNAseq_all_genes_new_scores_regress_age_sex.txt", sep = "\t", quote = FALSE, row.names = FALSE)


X_test_all_rnaseq <- as.matrix(df_test_all_rnaseq[, selected_genes])
all_rnaseq_adalasso_score <- as.vector(X_test_all_rnaseq %*% beta_selected_all_rnaseq_adalasso)

df_test_all_rnaseq$all_rnaseq_adalasso_score <- scale(all_rnaseq_adalasso_score)

formula_components <- c("all_rnaseq_adalasso_score", covars, 'offset(log(total_days_since_P2/365))')
formula_str <- paste0("has_exacerbation ~ ", paste0(formula_components, collapse = " + "))
model_all_genes <- glm(as.formula(formula_str), data = df_test_all_rnaseq, family = binomial())


# no covars
formula_str <- paste0("has_exacerbation ~ ", "all_rnaseq_adalasso_score")
model_all_genes_no_covars <- glm(as.formula(formula_str), data = df_test_all_rnaseq, family = binomial())


save(model_all_genes, model_all_genes_no_covars, df_test_all_rnaseq, file = "path/to/cg_RNAseq_all_genes_new_scores_w_covars_smk_exac_regress_age_sex_seed138_model2.RData")
