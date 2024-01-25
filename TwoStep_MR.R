setwd("<working_dir>")

library(TwoSampleMR)
library(MRPRESSO)
library(foreach)
library(doParallel)

######################

Step1_Exposure = TwoSampleMR::read_exposure_data(filename = "<Mental_Disorder_GWAS>",
                                           sep = "\t",
                                           snp_col = "SNP",
                                           effect_allele_col = "A1",
                                           other_allele_col = "A2",
                                           beta_col = "EFFECT_SIZE",
                                           se_col = "EFFECT_SIZE_STANDARD_ERROR",
                                           eaf_col = "MAF",
                                           samplesize_col = "N_All",
                                           pval_col = "PVALUE")

Step1_Exposure$exposure = "<Mental Disorder>"

Step1_Exposure_clumped = TwoSampleMR::clump_data(Step1_Exposure,
                                           clump_r2 = 0.01,
                                           clump_kb = 250,
                                           clump_p1 = 5e-08,
                                           clump_p2 = 5e-08,
                                           pop = "EUR")

Step1_Outcome = TwoSampleMR::read_outcome_data(filename = "<Clinical_Trait_GWAS>",
                                         sep = "\t",
                                         snp_col = "SNP",
                                         beta_col = "EFFECT_SIZE",
                                         se_col = "EFFECT_SIZE_STANDARD_ERROR",
                                         effect_allele_col = "ALT_ALLELE",
                                         other_allele_col = "REF_ALLELE",
                                         eaf_col = "MAF",
                                         samplesize_col = "N_All",
                                         pval_col = "PVALUE")

Step1_Harmonise_Data = TwoSampleMR::harmonise_data(exposure_dat = Step1_Exposure_clumped, 
                                             outcome_dat = Step1_Outcome,
                                             action = 2)

Step1_Harmonise_Data = Step1_Harmonise_Data[Step1_Harmonise_Data$mr_keep == T & Step1_Harmonise_Data$palindromic == F, ]

Step1_Outlier = MRPRESSO::mr_presso(BetaExposure = "beta.exposure", 
                              SdExposure = "se.exposure", 
                              BetaOutcome = "beta.outcome", 
                              SdOutcome = "se.outcome", 
                              OUTLIERtest = T, 
                              DISTORTIONtest = T, 
                              data = Step1_Harmonise_Data, 
                              NbDistribution = 5000,  
                              SignifThreshold = 0.05)

if(Step1_Outlier$`MR-PRESSO results`$`Global Test`$Pvalue < 0.05){
  Step1_Harmonise_Data = Harmonise_Data[Outlier$`MR-PRESSO results`$`Outlier Test`$Pvalue >= 0.05, ]
  Step1_Outlier = MRPRESSO::mr_presso(BetaExposure = "beta.exposure", 
                                SdExposure = "se.exposure", 
                                BetaOutcome = "beta.outcome", 
                                SdOutcome = "se.outcome", 
                                OUTLIERtest = T, 
                                DISTORTIONtest = T, 
                                data = Step1_Harmonise_Data, 
                                NbDistribution = 5000,  
                                SignifThreshold = 0.05)
}

Step1_Results = TwoSampleMR::mr(Step1_Harmonise_Data)
Step1_Heterogeneity = TwoSampleMR::mr_heterogeneity(Step1_Harmonise_Data)
Step1_Pleiotropy = TwoSampleMR::mr_pleiotropy_test(Step1_Harmonise_Data)

Step1_PVE = 2 * (1 - Step1_Harmonise_Data$eaf.exposure) * Step1_Harmonise_Data$eaf.exposure * Step1_Harmonise_Data$beta.exposure * Step1_Harmonise_Data2$beta.exposure
Step1_F_vaule = ((Step1_Harmonise_Data$samplesize.exposure[1] - nrow(Step1_Harmonise_Data) - 1)/nrow(Step1_Harmonise_Data)) * (Step1_PVE/(1 - Step1_PVE))

Step1_Directionality = TwoSampleMR::directionality_test(Step1_Harmonise_Data)


######################

Step2_Exposure = TwoSampleMR::read_exposure_data(filename = "<Clinical_Trait_GWAS>",
                                                 sep = "\t",
                                                 snp_col = "SNP",
                                                 effect_allele_col = "A1",
                                                 other_allele_col = "A2",
                                                 beta_col = "EFFECT_SIZE",
                                                 se_col = "EFFECT_SIZE_STANDARD_ERROR",
                                                 eaf_col = "MAF",
                                                 samplesize_col = "N_All",
                                                 pval_col = "PVALUE")

Step2_Exposure$exposure = "<Clinical Trait>"

Step2_Exposure_clumped = TwoSampleMR::clump_data(Step2_Exposure,
                                                 clump_r2 = 0.01,
                                                 clump_kb = 250,
                                                 clump_p1 = 5e-08,
                                                 clump_p2 = 5e-08,
                                                 pop = "EUR")

Step2_Outcome = TwoSampleMR::read_outcome_data(filename = "<Cancer_GWAS>",
                                               sep = "\t",
                                               snp_col = "SNP",
                                               beta_col = "EFFECT_SIZE",
                                               se_col = "EFFECT_SIZE_STANDARD_ERROR",
                                               effect_allele_col = "ALT_ALLELE",
                                               other_allele_col = "REF_ALLELE",
                                               eaf_col = "MAF",
                                               samplesize_col = "N_All",
                                               pval_col = "PVALUE")

Step2_Harmonise_Data = TwoSampleMR::harmonise_data(exposure_dat = Step2_Exposure_clumped, 
                                                   outcome_dat = Step2_Outcome,
                                                   action = 2)

Step2_Harmonise_Data = Step2_Harmonise_Data[Step2_Harmonise_Data$mr_keep == T & Step2_Harmonise_Data$palindromic == F, ]

Step2_Outlier = MRPRESSO::mr_presso(BetaExposure = "beta.exposure", 
                                    SdExposure = "se.exposure", 
                                    BetaOutcome = "beta.outcome", 
                                    SdOutcome = "se.outcome", 
                                    OUTLIERtest = T, 
                                    DISTORTIONtest = T, 
                                    data = Step2_Harmonise_Data, 
                                    NbDistribution = 5000,  
                                    SignifThreshold = 0.05)

if(Step2_Outlier$`MR-PRESSO results`$`Global Test`$Pvalue < 0.05){
  Step2_Harmonise_Data = Harmonise_Data[Outlier$`MR-PRESSO results`$`Outlier Test`$Pvalue >= 0.05, ]
  Step2_Outlier = MRPRESSO::mr_presso(BetaExposure = "beta.exposure", 
                                      SdExposure = "se.exposure", 
                                      BetaOutcome = "beta.outcome", 
                                      SdOutcome = "se.outcome", 
                                      OUTLIERtest = T, 
                                      DISTORTIONtest = T, 
                                      data = Step2_Harmonise_Data, 
                                      NbDistribution = 5000,  
                                      SignifThreshold = 0.05)
}

Step2_Results = TwoSampleMR::mr(Step2_Harmonise_Data)
Step2_Heterogeneity = TwoSampleMR::mr_heterogeneity(Step2_Harmonise_Data)
Step2_Pleiotropy = TwoSampleMR::mr_pleiotropy_test(Step2_Harmonise_Data)

Step2_PVE = 2 * (1 - Step2_Harmonise_Data$eaf.exposure) * Step2_Harmonise_Data$eaf.exposure * Step2_Harmonise_Data$beta.exposure * Step2_Harmonise_Data2$beta.exposure
Step2_F_vaule = ((Step2_Harmonise_Data$samplesize.exposure[1] - nrow(Step2_Harmonise_Data) - 1)/nrow(Step2_Harmonise_Data)) * (Step2_PVE/(1 - Step2_PVE))

Step2_Directionality = TwoSampleMR::directionality_test(Step2_Harmonise_Data)






