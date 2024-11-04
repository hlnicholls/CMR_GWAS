#!/bin/bash

phenotypes=("RVEDV" "RVESV" "RVSV" "RVEF" "RV_LV_ratio" "RVEDV_BSA" "RVESV_BSA" "RVSV_BSA")
for phenotype in "${phenotypes[@]}"
    do
    /opt/conda/envs/pops_env/bin/python3.7 /TOPmed_imputation/PostGWAS/10_Pops/pops/pops.py \
    --gene_annot_path /TOPmed_imputation/PostGWAS/10_Pops/pops/example/data/utils/gene_annot_jun10.txt \
    --feature_mat_prefix /TOPmed_imputation/PostGWAS/10_Pops/pops/features_munged_all/pops_features_all \
    --num_feature_chunks 2 \
    --magma_prefix /TOPmed_imputation/PostGWAS/9_MAGMA/Output/magma_${phenotype} \
    --control_features_path /TOPmed_imputation/PostGWAS/10_Pops/pops/example/data/utils/features_jul17_control.txt \
    --out_prefix /TOPmed_imputation/PostGWAS/10_Pops/Output/${phenotype}_pops_all_features
 done
