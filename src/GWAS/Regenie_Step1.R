library(tidyverse)
library(magrittr)

# Need to update paths
inpath <- '/Phenotypes_covariates/Output'
inpath2 <- '/Regenie_step1/Input'
regenie_path <- '/Software/regenie/regenie'

# Repeat for RV traits with updated file paths
system(str_interp(paste('${regenie_path}/regenie',
                        "--step 1",
                        "--bed ${inpath2}/ukb_directly_genotyped_qc_pass",
                        "--extract ${inpath2}/ukb_directly_genotyped_qc_pass.snplist", 
                        "--phenoFile ${inpath}/CMR_phenotypes_regenie_15112023.txt", 
                        "--phenoCol LVEDV,LVESV,LVSV,LVEF,LVM,LVMVR,LVGFI,LVMCF,LVEDV_BSA,LVESV_BSA,LVSV_BSA,LVM_BSA",
                        "--qt",
                        "--bsize 1000",
                        "--loocv",
                        "--threads 40",
                        "--out /Regenie/Output/UKB_step1_commonvariant_GrCh38")))
