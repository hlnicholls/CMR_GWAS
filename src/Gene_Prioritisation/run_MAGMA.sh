#!/bin/bash

# Repeat for LV phenotypes
phenotypes=("RVEDV" "RVESV" "RVSV" "RVEF" "RV_LV_ratio" "RVEDV_BSA" "RVESV_BSA" "RVSV_BSA")

for phenotype in "${phenotypes[@]}"
do
    /raid/genetics/Software/magma/magma \
    --bfile /TOPmed_imputation/PostGWAS/9_MAGMA/Input/1000G.EUR \
    --gene-annot /TOPmed_imputation/PostGWAS/9_MAGMA/Input/magma_0kb.genes.annot \
    --pval /TOPmed_imputation/PostGWAS/9_MAGMA/Input/${phenotype}_GWAS_38_rsid_magma.txt ncol=N \
    --gene-model snp-wise=mean \
    --out /TOPmed_imputation/PostGWAS/9_MAGMA/Output/magma_${phenotype}
done
