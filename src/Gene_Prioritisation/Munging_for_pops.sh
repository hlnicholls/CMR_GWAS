#!/bin/bash

/opt/conda/envs/pops_env/bin/python3.7 /TOPmed_imputation/PostGWAS/10_Pops/pops/munge_feature_directory.py \
 --gene_annot_path /TOPmed_imputation/PostGWAS/10_Pops/pops/example/data/utils/gene_annot_jun10.txt \
 --feature_dir /TOPmed_imputation/PostGWAS/10_Pops/pops/features_raw_selected/ \
 --save_prefix /TOPmed_imputation/PostGWAS/10_Pops/pops/features_munged_selected/pops_features_selected
