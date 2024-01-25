setwd("<<working_dir>>")

library(TwoSampleMR)
library(MRPRESSO)
library(foreach)
library(doParallel)

######################

Exposure_clumped = TwoSampleMR::mv_extract_exposures_local("<Mental_Disorder_GWAS and Covariate_GWAS>",
                                                           phenotype_col = "Phenotype",
                                                           id_col = "ID",
                                                           sep = "\t",
                                                           snp_col = "SNP",
                                                           effect_allele_col = "ALT_ALLELE",
                                                           other_allele_col = "REF_ALLELE",
                                                           beta_col = "EFFECT_SIZE",
                                                           se_col = "EFFECT_SIZE_STANDARD_ERROR",
                                                           eaf_col = "MAF",
                                                           samplesize_col = "N_All",
                                                           pval_col = "PVALUE",
                                                           min_pval = 1e-200,
                                                           log_pval = F,
                                                           pval_threshold = 5e-08,
                                                           clump_r2 = 0.01,
                                                           clump_kb = 250,
                                                           harmonise_strictness = 2)

Outcome = TwoSampleMR::read_outcome_data(filename = "<Cancer_GWAS>",
                                         sep = "\t",
                                         snp_col = "SNP",
                                         beta_col = "EFFECT_SIZE",
                                         se_col = "EFFECT_SIZE_STANDARD_ERROR",
                                         effect_allele_col = "ALT_ALLELE",
                                         other_allele_col = "REF_ALLELE",
                                         eaf_col = "MAF",
                                         samplesize_col = "N_All",
                                         pval_col = "PVALUE")

Outcome$outcome = "<Cancer Type>"

Harmonise_data = TwoSampleMR::mv_harmonise_data(Exposure_clumped, 
                                                Outcome, 
                                                harmonise_strictness = 2)

Results = TwoSampleMR::mv_multiple(Harmonise_data,
                                   intercept = T,
                                   instrument_specific = FALSE,
                                   pval_threshold = 5e-08,
                                   plots = FALSE)

