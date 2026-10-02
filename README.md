# PERI-COPD exacerbation biomarker analysis

This code derives and validates the **network-based exacerbation risk discrimination (NERD) score**, a 6-gene blood expression score for COPD exacerbation risk.

## Pipeline

| Step | File | Purpose |
|---|---|---|
| 1 | `impulseDE2_individual_top50perc.Rmd` | Find exacerbation-associated genes with ImpulseDE2, fit within each participant |
| 2 | `WGCNA.R` | Build the weighted gene co-expression network |
| 3 | `WGCNA_AMEND.R` | Identify the exacerbation module with AMEND (random walk with restart) |
| 4 | `make_pheno_AMEND_module_gene_plots.ipynb` | Plot module expression against daily EXACT symptom scores |
| 5 | `train_lasso_score.Rmd` | Select the NERD score genes from the module with adaptive LASSO |
| 6 | `train_de_novo_score_cg_RNAseq.R` | Train a de novo comparison score in COPDGene |
| 7 | `train_test_in_cg_RNAseq.Rmd` | Reweight the NERD score and validate it in COPDGene |
| 8 | `test_in_EC_microarray.Rmd` | Validate the NERD score in ECLIPSE |

### 1. Exacerbation-associated genes: `impulseDE2_individual_top50perc.Rmd`

1. **Filter genes.** Keep the top 50% most highly expressed genes.
2. **Align time to each exacerbation.** Each sample is re-indexed to weeks relative to exacerbation onset. Any week in which a participant met exacerbation criteria on one or more days is an exacerbation week (week 0).
3. **Fit models.** ImpulseDE2 is fit within each participant, adjusting for within-participant sequencing batch.
4. **Keep transient genes.** A gene is kept if it has a transient peak centered on an exacerbation (peak within ±1 week, width ≤ 2 weeks) with BH-adjusted P < .05.
5. **Keep replicated genes.** Genes found in at least 2 of 4 participants become the seed genes for the network.

This file also makes the Venn diagram, runs the Reactome pathway analysis (`sigora`), and draws the module heatmaps from the AMEND results.

### 2–3. Co-expression network and exacerbation module: `WGCNA.R`, `WGCNA_AMEND.R`

1. **Build the network.** WGCNA is built from exacerbation-week samples.
2. **Detect modules.** Modules are found by dynamic tree cut.
3. **Run AMEND.** The 5 modules with the most seed genes are pruned. AMEND then connects the seed genes within each module.
4. **Select the module.** The blue module is carried forward as the exacerbation module.

### 4. Module vs. symptoms: `make_pheno_AMEND_module_gene_plots.ipynb`

For each participant, this plots mean module expression (95% CI), scaled to the EXACT score range, against daily EXACT scores.

### 5. Score derivation: `train_lasso_score.Rmd`

- **Model.** Adaptive LASSO is fit to the exacerbation module genes across all PERI-COPD samples, with exacerbation week vs. non-exacerbation week as the outcome. Ridge regression supplies the adaptive penalty weights.
- **Cross-validation.** Participant-stratified, with one fold per participant. Class weights correct for how few exacerbation weeks there are.

### 6. De novo score: `train_de_novo_score_cg_RNAseq.R`

- **Purpose.** Trains a score on top 50% most expressed genes in COPDGene RNA-seq. It is a negative control for the NERD score.
- **Model.**
  1. Use the same 50:50 COPDGene training split as the NERD score.
  2. Fit adaptive LASSO, predicting any exacerbation during follow-up.

### 7–8. Validation: `train_test_in_cg_RNAseq.Rmd` and `test_in_EC_microarray.Rmd`

- **Reweighting.** COPDGene is split 50:50. In the training half, the coefficients are reweighted by logistic regression adjusting for age and sex. The trained score is then tested in the held-out test set of COPDGene and all of ECLIPSE.
- **Analyses.** Each analysis is run in the whole cohort and stratified by prior exacerbation:
  - multivariable logistic regression;
  - AUC (DeLong's test);
  - Cox regression;
  - Kaplan–Meier with Fleming–Harrington tests;
  - NERD quartiles analyses;
  - white blood cell adjustment (COPDGene only);
  - comparison with the de novo score.


## Setup
Paths written as `path/to/...` point to controlled-access individual-level data from PERI-COPD, COPDGene and ECLIPSE, which are not included. 

**R 4.3.3**

| Package | Version |
|---|---|
| ImpulseDE2 | 0.99.10 |
| glmnet | 5.0 |
| survival | 3.8-3 | 
| survminer | 0.5.2 |

**Python 3.8.18** 