library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(magrittr)

lv_df <- fread('Output/LV_Gene_annotation.csv')
rv_df <- fread('/Output/RV_Gene_annotation.csv')
lv_df_mtag <- fread('/Output/LV_mtag_Gene_annotation.csv')
rv_df_mtag <- fread('/Output/RV_mtag_Gene_annotation.csv')

lv_genes <- fread('/Output/all_LV_traits_loci_38_with_ld_genes.txt')
lvmtag_genes <- fread('/Output/all_LV_traits_loci_38_with_ld_genes.txt')
rv_genes <- fread('/Output/all_RV_traits_loci_38_with_ld_genes.txt')
rvmtag_genes <- fread('/Output/all_RV_traits_loci_38_with_ld_genes.txt')


lv_genes <-rbind(lv_df, lv_df_mtag, fill=TRUE)
rv_genes <- rbind(rv_df, rv_df_mtag, fill=TRUE)

# Merge LV and RV genes and condense to unique genes per row
merged_df <- merge(lv_genes, rv_genes, all=TRUE)

# Collapsing columns by gene
unique_loci_at_start <- unique(merged_df$Locus_name)
collapsed_df <- merged_df[, c(list(Phenotype = paste(unique(Phenotype), collapse = ", ")),
                              lapply(.SD, function(x) paste(unique(x), collapse = ", "))), 
                          by = Nearest_Gene_10kb, .SDcols = setdiff(names(merged_df), c("Nearest_Gene_10kb", "Phenotype"))]

collapsed_df <- filter(collapsed_df, !is.na(collapsed_df$Nearest_Gene_10kb))

# Looking at unique locus names
dt_split <- collapsed_df %>%
  separate_rows(Locus_name, sep = ",\\s*")
unique_loci_at_pops <- unique(dt_split$Locus_name) 

# Merging in annotations
pops_data <- data.table::fread( '/Output/Supplementary/supp_table10_all_pops_genes.csv', select=c('Nearest_Gene_10kb', 'Top_PoPS_Score_per_Locus'))
collapsed_df <- merge(collapsed_df, pops_data, all.x=TRUE)
pops_yes <- filter(collapsed_df , Top_PoPS_Score_per_Locus=='Yes') 

coloc_data <- data.table::fread('/Output/Supplementary/supp_table9_total_coloc.csv', select=c('Nearest_Gene_10kb', 'Locus_number'))

coloc_data$HF_Coloc_H4_08 <- "Yes"
coloc_data <- unique(coloc_data)
coloc_test <- unique(coloc_data$Nearest_Gene_10kb)

collapsed_df <- merge(collapsed_df, coloc_data, all.x=TRUE)
collapsed_df <- filter(collapsed_df, !is.na(collapsed_df$Nearest_Gene_10kb))
any(duplicated(collapsed_df$Nearest_Gene_10kb))

merged_df <- unique(collapsed_df)
merged_df <- filter(merged_df, !is.na(merged_df$Nearest_Gene_10kb))

drugnome <- fread('/Input/drugnomeai_LV_RV_results/supervised-learning/ranked-by-proba_predictions/GradientBoostingClassifier.All_genes.Ranked_by_prediction_proba.csv', 
                  header = FALSE, 
                  col.names = c('Nearest_Gene_10kb', 'DrugnomeAI_probability'))

merged_df <- merge(merged_df, drugnome, by='Nearest_Gene_10kb', all.x=T)
merged_df <- filter(merged_df, !is.na(merged_df$Nearest_Gene_10kb))

combine_values <- function(x) {
  unique_values <- unique(na.omit(x))
  if (length(unique_values) > 0) {
    return(paste(unique_values, collapse = ","))
  } else {
    return(NA)
  }
}

condensed_df <- merged_df %>% 
  group_by(Nearest_Gene_10kb) %>% 
  summarise(across(everything(), combine_values), .groups = "drop")
condensed_df <- filter(condensed_df, !is.na(condensed_df$Nearest_Gene_10kb))
any(duplicated(condensed_df$Nearest_Gene_10kb))

has_non_na_gtex <- function(row) {
  gtex_cols <- grep("GTEx", names(row), value = TRUE)
  any(!is.na(row[gtex_cols]))
}


check_for_specific_strings <- function(row, cols, target_strings) {
  any(sapply(row[cols], function(x) any(sapply(target_strings, function(y) grepl(y, x, ignore.case = TRUE)))))
}

condensed_df <- condensed_df %>%
  rowwise() %>%
  mutate(`has_GTEx_coloc_H4_>0.8` = {
    # Concatenate values from specified columns and split by comma
    values <- unlist(cur_data()[, c(32:55, 57:72)])
    numeric_values <- unlist(lapply(values, function(x) as.numeric(str_split(x, pattern = ",", simplify = TRUE))))
    # Check if any value is greater than 80
    result <- any(numeric_values >= 80, na.rm = TRUE)
    # Return 1 if true, else 0
    if_else(result, 1, 0)
  }) %>%
  ungroup()

na_row <- filter(condensed_df, is.na(condensed_df$Nearest_Gene_10kb))
condensed_df <- filter(condensed_df, !is.na(condensed_df$Nearest_Gene_10kb))
any(duplicated(condensed_df$Nearest_Gene_10kb))

mendel_count <- condensed_df %>%
  filter(grepl("cardio|cardiac|artial|myocard|arrhyt|vascular|heart|hypertroph|dilated|left ventric|right ventric", ClinVar_2019, ignore.case = TRUE) | 
           grepl("cardio|cardiac|artial|myocard|arrhyt|vascular|heart|hypertroph|dilated|left ventric|right ventric", OMIM_Expanded, ignore.case = TRUE))

clingen <- fread('/Input/clingen_gene_disease_validities_31JAN2024.csv')
setnames(clingen, c("Gene", "Disease", "Classification"), c('Nearest_Gene_10kb', "ClinGen_Disease", 'ClinGen_Classification'))
clingen <- dplyr::select(clingen,Nearest_Gene_10kb, ClinGen_Disease, ClinGen_Classification)
clingen <- unique(clingen)
clingen <- clingen %>% 
  group_by(Nearest_Gene_10kb) %>% 
  summarise(across(everything(), combine_values), .groups = "drop")

condensed_df <- merge(condensed_df, clingen, by ='Nearest_Gene_10kb', all.x=TRUE)
na_row <- filter(condensed_df, is.na(condensed_df$Nearest_Gene_10kb))
condensed_df <- filter(condensed_df, !is.na(condensed_df$Nearest_Gene_10kb))

clingen_count <- condensed_df %>%
  filter(
    (grepl("cardio|cardiac|artial|myocard|arrhyt|heart|hypertroph|dilated|left ventric|right ventric", ClinGen_Disease, ignore.case = TRUE) &
       grepl("definitive", ClinGen_Classification, ignore.case = TRUE)) |
      grepl("cardio|cardiac|artial|myocard|arrhyt|heart|hypertroph|dilated|left ventric|right ventric", ClinVar_2019, ignore.case = TRUE) | 
      grepl("cardio|cardiac|artial|myocard|arrhyt|heart|hypertroph|dilated|left ventric|right ventric", OMIM_Expanded, ignore.case = TRUE)
  )

dgidb_interact <- fread('/Input/DGidb/interactions.tsv')
setnames(dgidb_interact, c('gene_claim_name', 'drug_name', 'approved', 'interaction_score'), 
         c('Nearest_Gene_10kb','dgidb_drug_name', 'dgidb_approved_drug', 'dgidb_interaction_score' ))
dgidb_interact <- dplyr::select(dgidb_interact, Nearest_Gene_10kb, dgidb_drug_name, dgidb_approved_drug, dgidb_interaction_score)

dgidb_interact$dgidb_drug_name <- ifelse(is.na(dgidb_interact$dgidb_interaction_score),
                                         dgidb_interact$dgidb_drug_name,
                                         paste(dgidb_interact$dgidb_drug_name, 
                                               sprintf("(%0.4f)", as.numeric(dgidb_interact$dgidb_interaction_score))))

dgidb_interact <- dgidb_interact[dgidb_interact$dgidb_drug_name != 'NULL (NA)', ]

dgidb_interact <- dgidb_interact %>%
  filter(grepl("TRUE", dgidb_approved_drug))

dgidb_interact <- dgidb_interact %>%
  group_by(Nearest_Gene_10kb) %>%
  summarise(
    `DGIdb_drug(s) (interaction score)` = paste(dgidb_drug_name, collapse = ", "),
    dgidb_approved_drug = paste(dgidb_approved_drug, collapse = ", "),
    .groups = 'drop' 
  )


filtered_dgidb_interact <- dplyr::select(dgidb_interact, Nearest_Gene_10kb, `DGIdb_drug(s) (interaction score)`)

condensed_df <- merge(condensed_df, filtered_dgidb_interact, by='Nearest_Gene_10kb', all.x=TRUE)
na_row <- filter(condensed_df, is.na(condensed_df$Nearest_Gene_10kb))
condensed_df <- filter(condensed_df, !is.na(condensed_df$Nearest_Gene_10kb))

any(duplicated(condensed_df$Nearest_Gene_10kb))


condensed_df <- condensed_df %>%
  mutate(across(everything(), ~na_if(trimws(.), "")))

na_row <- filter(condensed_df, is.na(condensed_df$Nearest_Gene_10kb))
condensed_df <- filter(condensed_df, !is.na(condensed_df$Nearest_Gene_10kb))

impc <- fread('/Databases/IMPC/phenotypeHitsPerGene-v20_1.csv')
setnames(impc, c('Gene Symbol', 'Phenotype Hits'), c('Nearest_Gene_10kb', 'MGI_IMPC_Phenotypes'))
impc$Nearest_Gene_10kb <- toupper(impc$Nearest_Gene_10kb)
impc <- dplyr::select(impc, Nearest_Gene_10kb, MGI_IMPC_Phenotypes)

condensed_df <- merge(condensed_df, impc, all.x=TRUE)
condensed_df <- filter(condensed_df, !is.na(condensed_df$Nearest_Gene_10kb))


condensed_df$HF_Coloc_H4_08 <- ifelse(grepl("Yes", condensed_df$HF_Coloc_H4_08), "Yes", "No")

fwrite(condensed_df, '/Output/input_data_for_ranking.csv')

#################################
# Applying gene prioritisation criteria

output <- condensed_df %>%
  mutate(Pops = case_when(
    Top_PoPS_Score_per_Locus == 'Yes' ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(HiC_or_GTEx = case_when(
    (!is.na(HiC_tissue) | grepl('heart|artery', GTEx_Tissues_V8_2023, ignore.case = TRUE)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(MGI_or_OMIM_or_ClinGen = case_when(
    (grepl('heart|cardio|cardia|arrhyt', MGI_IMPC_Phenotypes, ignore.case = TRUE) |
       grepl('cardio|cardiac|artial|myocard|arrhyt|heart|hypertroph|dilated|left ventric|right ventric', ClinGen_Disease, ignore.case = TRUE) |
       grepl('cardio|cardiac|artial|myocard|arrhyt|heart|hypertroph|dilated|left ventric|right ventric', OMIM_Expanded, ignore.case = TRUE)) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(HF_coloc = case_when(
    HF_Coloc_H4_08 == "Yes" ~ 1,
    TRUE ~ 0
  )) 



cols_of_interest <- c('Pops', 'HiC_or_GTEx', 'MGI_or_OMIM_or_ClinGen', 'HF_coloc')

for (c in cols_of_interest){
  
  print(str_interp("${c} ${sum(output[[c]])}"))
  
}

output %<>% mutate(Gene_Prioritisation_Score=rowSums(output %>% dplyr::select(all_of(cols_of_interest))))

table(output$Gene_Prioritisation_Score)


coloc_genes <- filter(condensed_df, HF_Coloc_H4_08 == "Yes")

output <- output %>%
  mutate(across(where(is.character), ~sub("^,", "", .)))

output <- output %>%
  dplyr::select(-c(Locus_n))

output <- output %>%
  dplyr::select(Nearest_Gene_10kb, Gene_Prioritisation_Score, Top_PoPS_Score_per_Locus, everything())

clean_duplicates <- function(locus_name) {
  unique_genes <- unique(unlist(strsplit(locus_name, ",\\s*")))
  paste(unique_genes, collapse = ", ")
}

output$Locus_name <- sapply(output$Locus_name, clean_duplicates)

fwrite(output, '/Output/LV_RV_all_cols_15MAY24.csv')

any(duplicated(output$Nearest_Gene_10kb))

#################################
# Curating table columns

dt <- dplyr::select(output, Phenotype, CHROM, Nearest_Gene_10kb, Locus_name, Gene_Prioritisation_Score, `DGIdb_drug(s) (interaction score)`, Top_PoPS_Score_per_Locus, GO_Biological_Process_2023, GO_Cellular_Component_2023,         
                    GO_Molecular_Function_2023,  KEGG_2021_Human, GTEx_Tissues_V8_2023, GWAS_Catalog_2023, 
                    Exomiser_RV_Score, Exomiser_LV_Score, DrugnomeAI_probability, MGI_IMPC_Phenotypes, ClinGen_Disease, ClinVar_2019,
                    OMIM_Expanded, HiC_tissue, `has_GTEx_coloc_H4_>0.8`, HF_Coloc_H4_08)


lv_vep <- fread('/Output/vep_cleaned_results.txt')
rv_vep <- fread('/Output/vep_cleaned_results.txt')
lv_mtag_vep <- fread('/Output/vep_cleaned_results.txt')
rv_mtag_vep <- fread('/Output/vep_cleaned_results.txt')

vep <- rbind(rv_vep, lv_vep, lv_mtag_vep, rv_mtag_vep, fill=TRUE)
vep <- dplyr::select(vep, SYMBOL, Consequence)
vep <- unique(vep)
vep <- filter(vep, !grepl("intron|non_coding|upstream|downstream|UTR|intergenic", Consequence))
vep <- unique(vep)
colnames(vep)[c(1,2)] <- c('Nearest_Gene_10kb', 'Coding_variant_conequences')
condensed_vep <- vep %>% 
  group_by(Nearest_Gene_10kb) %>% 
  summarise(across(everything(), combine_values), .groups = "drop")
condensed_vep <- unique(condensed_vep)

any(duplicated(output$Nearest_Gene_10kb))
dt <- merge(dt, condensed_vep, by='Nearest_Gene_10kb', all.x=TRUE)
any(duplicated(output$Nearest_Gene_10kb))


combine_values_var <- function(values) {
  if (any(values == "Yes")) {
    return("Yes")
  } else {
    return(unique(values))
  }
}

combine_damaging_consequence <- function(values) {
  if (any(grepl("damaging|deleterious", values, ignore.case = TRUE))) {
    return(1)
  } else {
    return(0)
  }
}

var_annot <- fread('/Output/Supplementary/supp_table7_single_and_multi_lv_rv_variants_annotated.csv', 
                   select=c('Gene_Symbol', 'Consequence', 'CADD_and_Regulome_functional_filter', 'Damaging_Consequence_SIFT_or_PolyPhen'))
var_annot <- unique(var_annot)
setnames(var_annot, 'Gene_Symbol', 'Nearest_Gene_10kb')

condensed_var_annot <- var_annot %>% 
  group_by(Nearest_Gene_10kb) %>% 
  summarise(
    CADD_and_Regulome_functional_filter = combine_values_var(CADD_and_Regulome_functional_filter),
    Damaging_Consequence_SIFT_or_PolyPhen = combine_damaging_consequence(Damaging_Consequence_SIFT_or_PolyPhen),
    .groups = "drop"
  )
condensed_var_annot$CADD_and_Regulome_functional_filter[condensed_var_annot$CADD_and_Regulome_functional_filter == "" | is.na(condensed_var_annot$CADD_and_Regulome_functional_filter)] <- "No"
condensed_var_annot$Damaging_Consequence_SIFT_or_PolyPhen[condensed_var_annot$Damaging_Consequence_SIFT_or_PolyPhen == "" | is.na(condensed_var_annot$Damaging_Consequence_SIFT_or_PolyPhen)] <- 0

condensed_var_annot <- unique(condensed_var_annot)

dt <- merge(dt, condensed_var_annot, by='Nearest_Gene_10kb', all.x=TRUE)

any(duplicated(dt$Nearest_Gene_10kb))

####################################
# Adding druggability information

ot_drugs <- fread('/Input/OT-MAR-24-FTP/OT_drug_interactions.csv')
ot_warnings <- fread('/Input/OT-MAR-24-FTP/drugwarnings.csv')

ot_warnings$chemblIds <- gsub("\\['|'\\]|'", "", ot_warnings$chemblIds)
ot_warnings<- ot_warnings %>%
  separate_rows(chemblIds, sep = ",\\s*")
ot_warnings <- dplyr::select(ot_warnings, chemblIds, toxicityClass, country, description, efo_term)

ot_dt <- merge(ot_drugs, ot_warnings, by='chemblIds', all.x=TRUE, allow.cartesian=TRUE)
setnames(ot_dt, 'name_mechanisms', 'Nearest_Gene_10kb')
ot_dt <- unique(ot_dt)
ot_dt <- filter(ot_dt, !is.na(chemblIds))
ot_dt <- ot_dt[chemblIds != ""]
ot_pharmgkb <- fread('/Input/OT-MAR-24-FTP/pharmacogenomics.csv')
setnames(ot_pharmgkb, 'drugId', 'chemblIds')
ot_pharmgkb <- dplyr::select(ot_pharmgkb, chemblIds, datasourceId, pgxCategory, phenotypeText, directionality, variantRsId, isDirectTarget)
ot_dt <- merge(ot_dt, ot_pharmgkb, by='chemblIds', all.x=TRUE, allow.cartesian=TRUE)

ot_dt <- ot_dt %>% 
  group_by(Nearest_Gene_10kb) %>% 
  summarise(across(everything(), combine_values), .groups = "drop")


ot_dt_filt <- dplyr::select(ot_dt, Nearest_Gene_10kb, chemblIds, entity_drug, name_drug, terms_mechanisms, toxicityClass, pgxCategory, directionality, variantRsId, isDirectTarget)
ot_dt_filt <- dplyr::select(ot_dt_filt, Nearest_Gene_10kb,  name_drug)

ot_dt_filt <- as.data.frame(ot_dt_filt)
dt <- merge(dt, ot_dt_filt, by='Nearest_Gene_10kb', all.x=TRUE)

dgidb <- fread('/Input/DGidb/categories.tsv')
setnames(dgidb, 'entrez_gene_symbol', 'Nearest_Gene_10kb')
dgidb$Druggable_dgidb <- ifelse(dgidb$category=='DRUGGABLE GENOME',1,0)
dgidb <- dplyr::select(dgidb, Nearest_Gene_10kb, Druggable_dgidb)
dgidb <- unique(dgidb)
dgidb <- filter(dgidb, Druggable_dgidb==1)
dt <- merge(dt, dgidb, by='Nearest_Gene_10kb', all.x=TRUE)


setnames(dt,'name_drug', 'OT_drug(s)')

loci <- fread("/Output/Supplementary/supp_table2_all_loci.csv", select=c('Locus_name', 'Locus_number', 'Reported_Locus'))

# Split Locus_name in dt and create a long format for merging
dt_long <- dt %>%
  separate_rows(Locus_name, sep = ",\\s*")

# Merge with loci data.frame
merged_dt <- left_join(dt_long, loci, by = "Locus_name")

# Collapse back to the original structure
final_dt <- merged_dt %>%
  group_by(across(-c(Locus_name, Locus_number))) %>%
  summarise(Locus_number = toString(unique(na.omit(Locus_number))), .groups = 'drop')

# Reassign the original Locus_name values to final_dt
final_dt <- final_dt %>%
  left_join(dt %>% select(Nearest_Gene_10kb, Locus_name), by = "Nearest_Gene_10kb")


dt_output <- dplyr::select(final_dt, Phenotype, Nearest_Gene_10kb, Locus_number, Locus_name, Reported_Locus, CHROM, Top_PoPS_Score_per_Locus,
                           Gene_Prioritisation_Score, HF_Coloc_H4_08, MGI_IMPC_Phenotypes, ClinGen_Disease, OMIM_Expanded, 
                           HiC_tissue,GTEx_Tissues_V8_2023, `OT_drug(s)`, `DGIdb_drug(s) (interaction score)`, 
                           Druggable_dgidb,DrugnomeAI_probability, Exomiser_RV_Score, Exomiser_LV_Score, 
                           CADD_and_Regulome_functional_filter, Damaging_Consequence_SIFT_or_PolyPhen,
                           Coding_variant_conequences,`has_GTEx_coloc_H4_>0.8`, GO_Biological_Process_2023, GO_Cellular_Component_2023, GO_Molecular_Function_2023, KEGG_2021_Human,  GWAS_Catalog_2023,  ClinVar_2019) 

dt_split <- dt_output %>%
  separate_rows(Locus_name, sep = ",\\s*")
unique_loci_test <- unique(dt_split$Locus_name) 

test <-filter(dt_output, HF_Coloc_H4_08=='Yes')

setnames(dt_output, 'HF_Coloc_H4_08', 'has_HF_Coloc_H4_>0.8')
dt_output[, 8:ncol(dt_output)] <- lapply(dt_output[, 8:ncol(dt_output)], function(x) gsub("NA,\\s*", "", x))
dt_output[, 8:ncol(dt_output)] <- lapply(dt_output[, 8:ncol(dt_output)], function(x) gsub(",\\s*NA", "", x))

dt_output <- filter(dt_output, Nearest_Gene_10kb!="")
any(duplicated(dt_output$Nearest_Gene_10kb))
fwrite(dt_output, '/Output/Supplementary/supp_table11_LV_RV_gene_annotation.csv')
