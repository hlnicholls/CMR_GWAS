
inpath="/TOPmed_imputation/PostGWAS/6_LDSC/Input/"

phenotypes=("RVEDV" "RVESV" "RVSV" "RVEF" "RV_LV_ratio" "RVEDV_BSA" "RVESV_BSA" "RVSV_BSA")

for phenotype in "${phenotypes[@]}"
do
    other_phenotypes=()
    for item in "${phenotypes[@]}"; do
        if [[ "$item" != "$phenotype" ]]; then
            other_phenotypes+=("$item")
        fi
    done

    other_phenotypes_string=$(printf ",${inpath}%s_GWAS_37_corr.txt" "${other_phenotypes[@]}")
    other_phenotypes_string=${other_phenotypes_string:1}

    /opt/conda/envs/python2.7/bin/python2.7 /raid/genetics/Software/ldsc_v1.0.1/ldsc.py \
    --rg ${inpath}${phenotype}_GWAS_37_corr.txt,${other_phenotypes_string} \
    --ref-ld /Software/ldsc_v1.0.1/databases/UKBB.ALL.ldscore/UKBB.EUR \
    --w-ld /Software/ldsc_v1.0.1/databases/UKBB.ALL.ldscore/UKBB.EUR \
    --out /TOPmed_imputation/PostGWAS/6_LDSC/Output/genetic_correlation_ldsc_res_ukb_${phenotype}
done

