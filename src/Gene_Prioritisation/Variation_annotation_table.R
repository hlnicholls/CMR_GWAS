library(data.table)
library(dplyr)
library(stringr)
library(readxl)

var_lv <- fread('/Output/all_LV_loci_annotated.tsv')
var_lv <- dplyr::select(var_lv, -Locus_n, -TEST, -rsid)
setnames(var_lv, c('PP4', 'PP3'), c('GTEx_PP4', 'GTEx_PP3'))

var_rv <- fread('/Output/all_RV_loci_annotated.tsv')
var_rv <- dplyr::select(var_rv, -Locus_n, -TEST, -rsid)
setnames(var_rv, c('PP4', 'PP3'), c('GTEx_PP4', 'GTEx_PP3'))
var_lv_mtag <- fread('/Output/all_LV_mtag_loci_annotated.tsv')
var_lv_mtag <- dplyr::select(var_lv_mtag, -Locus_n, -Z, -ID, -rsid)
var_rv_mtag <- fread('/Output/all_mtag_RV_loci_annotated.tsv')
var_rv_mtag <- dplyr::select(var_rv_mtag, -Locus_n, -Z, -ID, -rsid)

var_lv$Method <- 'Single-trait'
var_rv$Method <- 'Single-trait'
var_lv_mtag$Method <- 'Multi-trait'
var_rv_mtag$Method <- 'Multi-trait'

rv_all_vars <- rbind(var_rv, var_rv_mtag, fill=TRUE)
lv_all_vars <- rbind(var_lv, var_lv_mtag, fill=TRUE)

regulome <- fread('/14_RegulomeDB/Input/ENCFF250UJY_Dec23.tsv')
regulome[, chrom := as.character(chrom)]
gwas <- rv_all_vars
gwas$chrom <- paste0('chr', gwas$CHROM)
gwas$start <- gwas$GENPOS
gwas$end <- gwas$GENPOS
setkey(gwas, chrom, start, end)
setkey(regulome, chrom, start, end)
regulome_gwas <- foverlaps(gwas, regulome, nomatch = 0)
regulome_gwas <- dplyr::select(regulome_gwas, SNP, rsid)
regulome_gwas <- regulome_gwas[!duplicated(regulome_gwas$SNP), ]
gwas_out <- merge(gwas, regulome_gwas, by ='SNP', all.x = TRUE)
gwas_out$chrom <- NULL
gwas_out$start <- NULL
gwas_out$end <- NULL
rv_out <- gwas_out

gwas <- lv_all_vars
gwas$chrom <- paste0('chr', gwas$CHROM)
gwas$start <- gwas$GENPOS
gwas$end <- gwas$GENPOS
setkey(gwas, chrom, start, end)
setkey(regulome, chrom, start, end)
regulome_gwas <- foverlaps(gwas, regulome, nomatch = 0)
regulome_gwas <- dplyr::select(regulome_gwas, SNP, rsid)
regulome_gwas <- regulome_gwas[!duplicated(regulome_gwas$SNP), ]
gwas_out <- merge(gwas, regulome_gwas, by ='SNP', all.x = TRUE)
gwas_out$chrom <- NULL
gwas_out$start <- NULL
gwas_out$end <- NULL
lv_out <- gwas_out

rv_out <- rv_out %>%
  dplyr::select(SNP, rsid, everything())
lv_out <- lv_out %>%
  dplyr::select(SNP, rsid, everything())

lv_out <- lv_out %>%
  filter(credible_set_99 == 'Yes' | (SNP == lead_snp & credible_set_99 == 'No' & credible_set_95 == 'No') | credible_set_95 == 'Yes' | R2 > 0.8)

rv_out <- rv_out %>%
  filter(credible_set_99 == 'Yes' | (SNP == lead_snp & credible_set_99 == 'No' & credible_set_95 == 'No') | credible_set_95 == 'Yes' | R2 > 0.8)

lv_out$locus38 <- sub("chr", "", lv_out$locus38)
rv_out$locus38 <- sub("chr", "", rv_out$locus38)

fwrite(lv_out, '/Output/single_and_multi_lv_variants_annotated.csv')
fwrite(rv_out, '/Output/single_and_multi_rv_variants_annotated.csv')

vars_all <- rbind(lv_out, rv_out)

vars_all$GTEx_PP3 <- vars_all$GTEx_PP3 / 100
vars_all$GTEx_PP4 <- vars_all$GTEx_PP4 / 100

fwrite(vars_all, '/Output/single_and_multi_all_variants_annotated.csv')

high_and_moderate_impact <-"transcript_ablation|splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant|stop_lost|start_lost|transcript_amplification|feature_elongation|feature_truncation|inframe_insertion|inframe_deletion|missense_variant|protein_altering_variant"

conseq <- dplyr::select(vars_all, Consequence, SIFT, PolyPhen)

# Filtering for specific variant types in Consequence
nonsyn_variants <- filter(conseq, grepl(high_and_moderate_impact, Consequence))

# Further filtering for deleterious variants
nonsyn_variants <- nonsyn_variants %>%
  filter(grepl("deleterious\\(", SIFT) | (grepl("damaging", PolyPhen) & !grepl("possibly", PolyPhen, ignore.case = TRUE)))

# Creating a new dataframe with an additional column for damaging variant annotation
vars_all <- vars_all %>%
  mutate(CADD_RawScore_Scaled = scale(CADD_RawScore))

# Use the scaled CADD_RawScore in the filtering
annotated_variants <- vars_all %>%
  mutate(CADD_RawScore_Scaled = scale(CADD_RawScore),
         Damaging_Consequence_SIFT_or_PolyPhen = ifelse(grepl(high_and_moderate_impact, Consequence) & 
                                                          (grepl("deleterious\\(", SIFT) | (grepl("damaging", PolyPhen) & !grepl("possibly", PolyPhen))), 1, 0),
         CADD_and_Regulome_functional_filter = ifelse(CADD_RawScore_Scaled > 20 | grepl("1", regulome_ranking), "Yes", "No"),
         Functional_and_Damaging_Prediction = ifelse(CADD_and_Regulome_functional_filter == "Yes" & Damaging_Consequence_SIFT_or_PolyPhen == 1, "Yes", "No"))

nonsyn_variants_damamging <- filter(annotated_variants, grepl(high_and_moderate_impact, Consequence))
total_yes <- nonsyn_variants_damamging %>% # 27 
  filter(Functional_and_Damaging_Prediction == "Yes") %>%
  distinct(SNP) %>%
  nrow()

setDT(annotated_variants)
annotated_variants[is.na(BETA) & !is.na(mtag_beta), `:=`(
  BETA = mtag_beta,
  SE = mtag_se,
  P = mtag_pval
)]

annotated_variants[, `:=`(
  mtag_beta = NULL,
  mtag_se = NULL,
  mtag_pval = NULL
)]

setnames(annotated_variants, c('locus38', 'locus37', 'lead_snp', 'SNP', 'rsid'), c('BP37', 'BP38','Lead_SNP','CredibleSet_SNP', 'CredibleSet_SNP_rsid'))

all_loci <- fread('/Output/Supplementary/supp_table2_all_loci.csv', select=c('Locus_name', 'Lead_SNP', 'Locus_number', 'Reported_Locus'))

df_annotated <- merge(annotated_variants, all_loci, all.x=TRUE, allow.cartesian=TRUE)
df_annotated <- unique(df_annotated)

lookup <- df_annotated[Gene_Symbol != "", .(Gene_Symbol = first(Gene_Symbol)), by = CredibleSet_SNP]
df_annotated <- merge(df_annotated, lookup, by = "CredibleSet_SNP", all.x = TRUE, suffixes = c("", "_lookup"))
df_annotated[Gene_Symbol == "", Gene_Symbol := Gene_Symbol_lookup]
df_annotated[, Gene_Symbol_lookup := NULL]

variants_output <- dplyr::select(df_annotated, Phenotype, Locus_name, Locus_number, Reported_Locus, Lead_SNP, CredibleSet_SNP, CredibleSet_SNP_rsid, R2, credible_set_99,        
                                 credible_set_95, Gene_Symbol, CHROM,  GENPOS,  BP38, BP37,  MAF,   BETA,  SE,   P, Method,
                                  Chromatin_accessibility, PWM,                     Footprint,              
                                  QTL,                     CADD_PHRED,              CADD_RawScore,  regulome_ranking,        regulome_probability,      CADD_and_Regulome_functional_filter,   Damaging_Consequence_SIFT_or_PolyPhen, Functional_and_Damaging_Prediction,Consequence, 
                                 REVEL, SIFT,                    PolyPhen,                CLIN_SIG,                GTEx_eQTL_COLOC_Tissue,  GTEx_PP3,               
                                 GTEx_PP4,                Levin_HF_PP_H4                                  
                                 )

variants_output <- unique(variants_output)
variants_output[] <- lapply(variants_output, function(x) gsub("\\s*,-|-,\\s*", "", x))
variants_output <- filter(variants_output, !is.na(variants_output$Locus_number))

phenotype_order <- c("LVEDV", "LVESV", "LVSV", "LVEF", "LVGFI", "LVMCF", "LVM", "LVMVR", 
                     "RVEDV", "RVESV", "RVSV", "RVEF", "RV_LV_ratio", 
                     "LVEDV_BSA", "LVESV_BSA", "LVSV_BSA", "LVM_BSA", 
                     "RVEDV_BSA", "RVESV_BSA", "RVSV_BSA")

variants_output[, Phenotype := factor(Phenotype, levels = phenotype_order)]
variants_output <- variants_output[order(Locus_number, Phenotype)]


variants_output <- unique(variants_output)

variants_output[variants_output == ""] <- NA

fwrite(variants_output, '/Output/Supplementary/supp_table7_10MAY24_single_and_multi_lv_rv_variants_annotated.csv')

##########################################################################################################################################

functional <- filter(variants_output, CADD_and_Regulome_functional_filter=="Yes")
functional_unique <- functional[!duplicated(functional$CredibleSet_SNP), ] # 11957

exonic_conseq <-"synonymous_variant|missense_variant|inframe_insertion|inframe_deletion|stop_gained|frameshift_variant|coding_sequence_variant|stop_lost|stop_retained_variant|incomplete_terminal_codon_variant|start_retained_variant|start_lost"
intronic_conseq <-"splice_acceptor_variant|splice_donor_variant|splice_donor_5th_base_variant|splice_region_variant|splice_donor_region_variant|splice_polypyrimidine_tract_variant|intron_variant"
other_conseq <- "transcript_ablation|transcript_amplification|feature_elongation|feature_truncation|protein_altering_variant|mature_miRNA_variant|5_prime_UTR_variant|3_prime_UTR_variant|non_coding_transcript_exon_variant|NMD_transcript_variant|non_coding_transcript_variant|coding_transcript_variant|upstream_gene_variant|downstream_gene_variant|TFBS_ablation|TFBS_amplification|TF_binding_site_variant|regulatory_region_ablation|regulatory_region_amplification|regulatory_region_variant|intergenic_variant|sequence_variant"

# Checking variant consequences in the functional impact (CADD and regulome variants)
conseq <- filter(functional_unique, !is.na(Consequence))
conseq <- dplyr::select(conseq, Consequence, SIFT, PolyPhen)
conseq <- filter(conseq, Consequence != "")


intronic_variants <- conseq  %>%
  filter(grepl(intronic_conseq, Consequence)) #1080
nrow(intronic_variants)

intergenic_variants <- conseq %>%
  filter(grepl("intergenic_variant", Consequence)) #187
nrow(intergenic_variants)

# Count variant consequence types for all variants - 19116 variants with 48235 consequences
out <- rbind(lv_out, rv_out)
out_unique <- out[!duplicated(out$SNP), ]

consequence_counts <- table(unlist(strsplit(out_unique$Consequence, ",")))
total_sum <- sum(consequence_counts) #48235 

#consequence_counts
# 3_prime_UTR_variant                 5_prime_UTR_variant             coding_sequence_variant             downstream_gene_variant 
# 421                                 108                                   1                                4812 
# intergenic_variant                      intron_variant                    missense_variant              NMD_transcript_variant 
# 2838                               21732                                 145                                3851 
# non_coding_transcript_exon_variant       non_coding_transcript_variant             splice_acceptor_variant       splice_donor_5th_base_variant 
# 1086                                8195                                   2                                   4 
# splice_donor_region_variant                splice_donor_variant splice_polypyrimidine_tract_variant               splice_region_variant 
# 17                                   1                                  67                                  50 
# stop_gained                           stop_lost                  synonymous_variant               upstream_gene_variant 
# 3                                   2                                 161                                4739 

consequence_counts <- c( "3_prime_UTR_variant" = 421, "5_prime_UTR_variant" = 108, "coding_sequence_variant" = 1, "downstream_gene_variant" = 4812, "intergenic_variant" = 2838, "intron_variant" = 21732, "missense_variant" = 145, "NMD_transcript_variant" = 3851, "non_coding_transcript_exon_variant" = 1086, "non_coding_transcript_variant" = 8195, "splice_acceptor_variant" = 2, "splice_donor_5th_base_variant" = 4, "splice_donor_region_variant" = 17, "splice_donor_variant" = 1, "splice_polypyrimidine_tract_variant" = 67, "splice_region_variant" = 50, "stop_gained" = 3, "stop_lost" = 2, "synonymous_variant" = 161, "upstream_gene_variant" = 4739)

exonic_conseq <- c("synonymous_variant", "missense_variant", "inframe_insertion", "inframe_deletion", "stop_gained", "frameshift_variant", "coding_sequence_variant", "stop_lost", "stop_retained_variant", "incomplete_terminal_codon_variant", "start_retained_variant", "start_lost")
intronic_conseq <- c("splice_acceptor_variant", "splice_donor_variant", "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant", "intron_variant")
other_conseq <- c("transcript_ablation", "transcript_amplification", "feature_elongation", "feature_truncation", "protein_altering_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification", "regulatory_region_variant",  "sequence_variant")


exonic_conseq_no_synon <- c("synonymous_variant","missense_variant", "inframe_insertion", "inframe_deletion", "stop_gained", "frameshift_variant", "coding_sequence_variant", "stop_lost", "stop_retained_variant", "incomplete_terminal_codon_variant", "start_retained_variant", "start_lost")
filtered_variants <- variants_output %>%
  filter(str_detect(Consequence, str_c(exonic_conseq_no_synon, collapse = "|")))
exonic_functional_variants  <- filter(filtered_variants, Functional_and_Damaging_Prediction=="Yes")
exonic_functional_variants <- exonic_functional_variants [!duplicated(exonic_functional_variants $CredibleSet_SNP), ]
exonic_loci <- unique(exonic_functional_variants$Locus_number) #19 loci with exonic variants with functional consequences
exonic_genes <- unique(exonic_functional_variants$Gene_Symbol) # 21 genes in the loci

total_count <- sum(consequence_counts)

exonic_count <- sum(consequence_counts[names(consequence_counts) %in% exonic_conseq])
intronic_count <- sum(consequence_counts[names(consequence_counts) %in% intronic_conseq])
other_count <- sum(consequence_counts[names(consequence_counts) %in% other_conseq])

exonic_percentage <- (exonic_count / total_count) * 100
intronic_percentage <- (intronic_count / total_count) * 100
other_percentage <- (other_count / total_count) * 100 # excluded intergenic
intergenic <- (2838/total_count) * 100

result <- data.frame(
  Group = c("Exonic", "Intronic", "Other"),
  Count = c(exonic_count, intronic_count, other_count),
  Percentage = c(exonic_percentage, intronic_percentage, other_percentage)
)

result

# Group Count Percentage
# 1   Exonic   312  0.6468332
# 2 Intronic 21873 45.3467399
# 3    Other 23212 48.1227325
# > exonic_count
# [1] 312
