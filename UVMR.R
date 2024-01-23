setwd("<working_dir>")

library(TwoSampleMR)
library(MRPRESSO)
library(foreach)
library(doParallel)

######################


Exposure = TwoSampleMR::read_exposure_data(filename = "<Mental_Disorder_GWAS>",
                                           sep = "\t",
                                           snp_col = "SNP",
                                           effect_allele_col = "A1",
                                           other_allele_col = "A2",
                                           beta_col = "EFFECT_SIZE",
                                           se_col = "EFFECT_SIZE_STANDARD_ERROR",
                                           eaf_col = "MAF",
                                           samplesize_col = "N_All",
                                           pval_col = "PVALUE")

Exposure$exposure = "<Mental Disorder>"

Exposure_clumped = TwoSampleMR::clump_data(Exposure,
                                           clump_r2 = 0.01,
                                           clump_kb = 250,
                                           clump_p1 = 5e-08,
                                           clump_p2 = 5e-08,
                                           pop = "EUR")

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

Harmonise_Data = TwoSampleMR::harmonise_data(exposure_dat = Exposure_clumped, 
                                             outcome_dat = Outcome,
                                             action = 2)

Harmonise_Data = Harmonise_Data[Harmonise_Data$mr_keep == T & Harmonise_Data$palindromic == F, ]

Outlier = MRPRESSO::mr_presso(BetaExposure = "beta.exposure", 
                              SdExposure = "se.exposure", 
                              BetaOutcome = "beta.outcome", 
                              SdOutcome = "se.outcome", 
                              OUTLIERtest = T, 
                              DISTORTIONtest = T, 
                              data = Harmonise_Data, 
                              NbDistribution = 5000,  
                              SignifThreshold = 0.05)

if(Outlier$`MR-PRESSO results`$`Global Test`$Pvalue < 0.05){
  Harmonise_Data = Harmonise_Data[Outlier$`MR-PRESSO results`$`Outlier Test`$Pvalue >= 0.05, ]
  Outlier = MRPRESSO::mr_presso(BetaExposure = "beta.exposure", 
                                SdExposure = "se.exposure", 
                                BetaOutcome = "beta.outcome", 
                                SdOutcome = "se.outcome", 
                                OUTLIERtest = T, 
                                DISTORTIONtest = T, 
                                data = Harmonise_Data, 
                                NbDistribution = 5000,  
                                SignifThreshold = 0.05)
}

Results = TwoSampleMR::mr(Harmonise_Data)
Heterogeneity = TwoSampleMR::mr_heterogeneity(Harmonise_Data)
Pleiotropy = TwoSampleMR::mr_pleiotropy_test(Harmonise_Data)

PVE = 2 * (1 - Harmonise_Data$eaf.exposure) * Harmonise_Data$eaf.exposure * Harmonise_Data$beta.exposure * Harmonise_Data2$beta.exposure
F_vaule = ((Harmonise_Data$samplesize.exposure[1] - nrow(Harmonise_Data) - 1)/nrow(Harmonise_Data)) * (PVE/(1 - PVE))

Directionality = TwoSampleMR::directionality_test(Harmonise_Data)




