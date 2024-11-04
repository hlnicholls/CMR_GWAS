library(tidyverse)
library(magrittr)

plink_path <- "/Software/plink1.9"
plink2_path <- "/Software/plink2"

file_path <- "/LV_GWAS/Regenie_step1/Input"

#Following recommendation from https://rgcgithub.github.io/regenie/recommendations/

#Pruning:
#window size is 1000kb, step size (50 SNPs), and LD threshold (R^2 < 0.4).
#pruning from https://github.com/pjgreer/ukb-rap-tools/blob/main/GWAS_pipeline/GTfile_prep/03-GTprep-ldprune.sh

# Genotype files need to be GRCh38

system(str_interp(paste("${plink2_path}/plink2_Jan23",
                        "--bfile /Input/all_chr_genotyped_keep_allele_order_GrCh38",
                        "--indep-pairwise 1000 50 0.4 ",
                        "--out ${file_path}/ukb-pruning")))

system(str_interp(paste("${plink2_path}/plink2_Jan23",
                        "--bfile /Input/all_chr_genotyped_keep_allele_order_GrCh38",
                        "--extract ${file_path}/ukb-pruning.prune.in",
                        "--make-bed",
                        "--out ${file_path}/ukb_directly_genotyped_pruned_cohort")))

# Variant QC filters:
#QC of directly genotype variants to extract a set of high-qual variants
system(str_interp(paste("${plink2_path}/plink2",
                        "--bfile ${file_path}/ukb_directly_genotyped_pruned_cohort",
                        "--maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1",
                        "--write-snplist --write-samples --no-id-header --make-bed",
                        "--out ${file_path}}/ukb_directly_genotyped_qc_pass")))
