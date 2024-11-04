library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(stringr)

rv <- read_excel('/PostGWAS/18_Pleiotropy/Output/RV_all_loci_in_LD_phenoscanner_output_R20.8.xlsx', sheet = 1)
rv_mtag <- read_excel('/PostGWAS/21_MTAG_loci/18_Pleiotropy/Output/RV_mtag_all_loci_in_LD_phenoscanner_output_R20.8.xlsx')
lv <- read_excel('/PostGWAS/18_Pleiotropy/Output/RV_all_loci_in_LD_phenoscanner_output_R20.8.xlsx', sheet = 1)
lv_mtag <- read_excel('/18_Pleiotropy/Output/LV_mtag_all_loci_in_LD_phenoscanner_output_R20.8.xlsx')

df <- rbind(rv, lv, rv_mtag, lv_mtag)
df <- as.data.frame(df)
df <- unique(df)

vars <- fread('/Output/Supplementary/supp_table7_10MAY24_single_and_multi_lv_rv_variants_annotated.csv')
setnames(vars, 'CredibleSet_SNP_rsid', 'rsid')
vars <- dplyr::select(vars, Phenotype, Locus_name, Locus_number, Lead_SNP, CredibleSet_SNP, rsid, Gene_Symbol, CHROM,  GENPOS,  BP38,    
                      BP37, MAF,    BETA,                      
                      SE,   P)
vars <- unique(vars)

vars_with_rsids <- filter(vars, rsid %in% df$rsid)

dt <- merge(vars_with_rsids, df, all.x=TRUE,  allow.cartesian=TRUE)


dt <- dplyr::select(dt, Phenotype, Locus_name, Locus_number, Lead_SNP, CredibleSet_SNP, rsid, Gene_Symbol, CHROM,  GENPOS,  BP38,    
                           BP37, MAF,    BETA,                      
                           SE,   P, Phenoscanner_Trait, Study, PMID,
                           Phenoscanner_Ancestry, Year,   Phenoscanner_A1,                 Phenoscanner_A2,                
                           Phenoscanner_BETA,               Phenoscanner_SE,                 Phenoscanner_P,                  Phenoscanner_Direction,         
                           Phenoscanner_N,                  Phenoscanner_Unit,               Phenoscanner_Dataset)


dt <- unique(dt)

lv_pheno <- c("LVEDV",       "LVESV",       "LVSV",        "LVEF",        "LVGFI",       "LVMCF",
              "LVM",         "LVMVR",  "LVEDV_BSA",   "LVESV_BSA",   "LVSV_BSA",    "LVM_BSA")      
rv_pheno <- c("RVEDV",       "RVESV",      "RVSV",        "RVEF",        "RV_LV_ratio",    "RVEDV_BSA",   "RVESV_BSA",   "RVSV_BSA")


rv_df <- filter(dt, Phenotype %in% rv_pheno)
lv_df <- filter(dt, Phenotype %in% lv_pheno)

rv_loci_count <- unique(rv_df$lead_snp) 
lv_loci_count <- unique(lv_df$lead_snp) 

rv_traits_count <- unique(rv_df$Phenoscanner_Trait)
lv_traits_count <- unique(lv_df$Phenoscanner_Trait) 
rv_traits <- rv_df$Phenoscanner_Trait
lv_traits <- lv_df$Phenoscanner_Trait
total_unique_traits <- unique(c(rv_traits, lv_traits))
total_unique_count <- length(total_unique_traits)
overlap_traits <- intersect(rv_traits, lv_traits)
overlap_count <- length(overlap_traits)
cat("Total unique Phenoscanner_Traits across both dataframes:", total_unique_count, "\n")
cat("Number of overlapping Phenoscanner_Traits:", overlap_count)


key_phrases <- list(
  Smoking = c(
    "Smoking status: previous", 
    "Past tobacco smoking"
  ),
  
  `Respiratory Traits` = c(
    "lung function", 
    "asthma", 
    "respiratory", 
    "expiratory volume", 
    "peak expiratory flow",
    "Forced vital capacity", 
    "Wheeze or whistling in the chest in last year", 
    "Snoring", 
    "Treatment with ventolin 100micrograms inhaler"
  ),
  
  `Renal Function` = c(
    "Unspecified haematuria", 
    "Chronic kidney disease", 
    "Sodium in urine", 
    "Neoplasm of uncertain or unknown behaviour of endocrine glands", 
    "Kidney diseases", 
    "Serum creatinine", 
    "Creatinine levels", 
    "Glomerular filtration rate", 
    "Glomerular filtration rate creatinine", 
    "Glomerular filtration rate in non diabetics creatinine", 
    "Creatinine", 
    "Renal function related traits sCR", 
    "Renal function related traits eGRFcrea",
    "Idiopathic membranous nephropathy",
    "Glomerular filtration rate in non diabetics creatinine",
    "Cystatin C in serum",
    "log eGFR cystatin C"
  ),
  
  `Metabolic Traits` = c(
    "metabolic rate", 
    "glycated hemoglobin", 
    "insulin", 
    "C reactive protein", 
    "Fibrinogen levels", 
    "Triglyceride levels", 
    "HbA1c", 
    "Metabolic syndrome domains Multivariate analysis", 
    "Adiponectin levels", 
    "Adiponectin",
    "Serum urate",
    "Urate levels",
    "diabetes", 
    "Treatment with metformin"
  ),
  
  Lipids = c(
    "cholesterol", 
    "Triglycerides", 
    "Low density lipoprotein", 
    "High density lipoprotein"
  ),
  
  Hematology = c(
    "Blood protein levels",
    "Monocyte percentage of white cells", 
    "Mean corpuscular volume",
    "Beta 2 microglubulin plasma levels",
    "blood cell count", 
    "hemoglobin", 
    "platelet", 
    "eosinophil", 
    "neutrophil", 
    "lymphocyte", 
    "myeloid", 
    "red cell distribution width", 
    "Hematocrit", 
    "Reticulocyte", 
    "Granulocyte count", 
    "Basophil count", 
    "High light scatter", 
    "Plasminogen activator inhibitor 1 PAI 1 concentration", 
    "Red blood cell traits", 
    "Hematological parameters", 
    "Basophil percentage of granulocytes",
    "Blood metabolite levels",
    "Monocyte count",
    "Erythrocyte indices",
    "Plasma Beta 2 microglobulin levels",
    "Plasminogen activator inhibitor type 1 levels PAI 1",
    "Gamma-linolenic acid 18:3n6"  
  ),
  
  `Dilated Cardiomyopathy` = c(
    "Sporadic dilated cardiomyopathy", 
    "Dilated cardiomyopathy", 
    "Cardiomyopathy dilated"
  ),
  
  `Cognitive Traits` = c(
    "neuroticism", 
    "cognitive", 
    "IQ", 
    "intelligence", 
    "mood swings", 
    "anxiety", 
    "depression", 
    "psychological", 
    "Irritability", 
    "Nervous feelings", 
    "Worrier or anxious feelings", 
    "Worry too long after embarrassment", 
    "Tense or highly strung", 
    "Fed-up feelings", 
    "Miserableness", 
    "Sensitivity or hurt feelings", 
    "Suffer from nerves", 
    "Mean time to correctly identify matches", 
    "Relative age of first facial hair", 
    "Overall health rating", 
    "Schizophrenia", 
    "Extraversion", 
    "Sleep duration", 
    "Bipolar disorder or attention deficit hyperactivity disorder",
    "Time spent watching television",
    "Information processing speed"
  ),
  
  `Cardiovascular Traits` = c(
    "Aorta",
    "Left ventricular internal diastolic dimensions LV internal diastolic dimensions",
    "Left ventricular diastolic dimensions LV diastolic dimensions",
    "Arrhythmias cardiac", 
    "heart rate", 
    "P wave duration",
    "RR interval", 
    "QT interval", 
    "QT interval females", 
    "QT interval old age 50 yrs", 
    "PR duration", 
    "PR segment", 
    "QRS interval",
    "P wave terminal force",
    "Electrocardiography", 
    "Electrocardiographic conduction measures", 
    "Electrocardiographic traits",
    "Pulse rate",
    "Pulse wave reflection index",
    "Left ventricle diastolic internal dimension",
    "QRS complex 12 leadsum",
    "Heart function tests",
    "PR interval",
    "QRS complex Cornell",
    "QRS complex Sokolow Lyon",
    "QRS duration",
    "Aortic root size",
    "Cardiac structure and function"
  ),
  
  `Cardiovascular Diseases` = c(
    "cardiovascular", 
    "atrial fibrillation", 
    "coronary artery", 
    "myocardial infarction", 
    "Ischemic stroke", 
    "Ischemic stroke small artery occlusion", 
    "Vascular or heart problems diagnosed by doctor",
    "Treatment with bendroflumethiazide", 
    "Treatment with atenolol", 
    "Treatment with amlodipine", 
    "Treatment with doxazosin",
    "Illnesses of father: heart disease",
    "Tetralogy of Fallot", 
    "Tetrology of fallot"
  ),
  
  `Body Size` = c(
    "height", 
    "sitting height"
  ),
  
  `Body Composition` = c(
    "weight", 
    "body fat", 
    "BMI", 
    "hip circumference", 
    "bone mineral density", 
    "impedance", 
    "fat-free mass", 
    "predicted mass", 
    "Whole body water mass", 
    "Body mass index", 
    "Arm fat mass", 
    "Leg fat mass", 
    "Trunk fat mass", 
    "Arm fat percentage", 
    "Comparative body size",
    "Waist circumference", 
    "Waist hip ratio", 
    "Comparative body size at age 10", 
    "Arm fat percentage left", 
    "Arm fat percentage right",
    "Waist circumference adjusted for smoking", 
    "Waist hip ratio adjusted for smoking",
    "Waist circumference adjusted for smoking in females", 
    "Waist hip ratio in physically active females",
    "Waist circumference in non-smokers", 
    "Waist hip ratio in females",
    "Waist hip ratio in physically active individuals", 
    "Waist circumference in female non-smokers",
    "Waist circumference in physically active individuals", 
    "Sexual dimorphism in anthropometric traits"
  ),
  
  `Blood Pressure` = c(
    "blood pressure", 
    "Hypertension", 
    "Pulse pressure", 
    "Mean arterial pressure"
  ),
  
  `Autoimmune Diseases` = c(
    "IgA deficiency", 
    "lupus", 
    "Arthritis rheumatoid",
    "Allopurinol related Stevens Johnson syndrome or toxic epidermal necrolysis", 
    "Behcet syndrome",
    "Autoimmune", 
    "Ankylosing spondylitis", 
    "Self-reported ankylosing spondylitis", 
    "Primary sclerosing cholangitis",
    "Self-reported sarcoidosis", 
    "Multiple sclerosis", 
    "Self-reported multiple sclerosis",
    "Autism spectrum disorder or schizophrenia", 
    "Rheumatoid arthritis", 
    "Seropositive rheumatoid arthritis",
    "Juvenile idiopathic arthritis including oligoarticular and rheumatoid factor negative polyarticular JIA",
    "Celiac disease",
    "Other rheumatoid arthritis", 
    "Self-reported rheumatoid arthritis", 
    "Rheumatoid arthritis ACPA positive", 
    "Behcet syndrome",
    "Primary biliary cirrhosis",
    "Thyroid peroxidase antibody positivity"
  ),
  
  Miscellaneous_Conditions = c(
    "skin",  
    "psoriasis", 
    "eczema", 
    "allergic", 
    "immunoglobulin", 
    "Vitiligo",
    "Atopic dermatitis", 
    "Allergic disease", 
    "Itch intensity from mosquito bite adjusted by bite size",
    "Taking other prescription medications",
    "Number of self-reported non-cancer illnesses",
    "Number of treatments or medications taken",
    "Illnesses of siblings: none of the above, group 1",
    "Treatment with thyroxine product",
    "Central retinal venular equivalent",
    "Retinal vascular caliber",
    "Retinal vein",
    "Crohns disease",
    "Qualifications: none",
    "Biomedical quantitative traits",
    "Treatment with dovobet ointment",
    "Self-reported psoriatic arthropathy",
    "Medication for pain relief, constipation, heartburn: none of the above",
    "Age-related macular degeneration",
    "Kawasaki disease excluding controls with chronic hepatitis B drug eruption and tuberculosis",
    "Kawasaki disease",
    "BLK gene expression in transformed B cell lines",
    "C8orf13 gene expression in transformed B cell lines",
    "Gene expression of GOSR2 in liver",
    "Breast cancer",
    "Gene expression",
    "Barretts esophagus",
    "HIV 1 control",
    "Mucocutaneous lymph node syndrome",
    "Mouth or teeth dental problems: mouth ulcers",
    "C8orf5FAM167ABLK expression in lymphoblastoid cell lines",
    "Proximal gene expression in monocytes II UFSP1 gene",
    "Proximal gene expression in monocytes UFSP1 gene",
    "Proximal gene expression in monocytes TRIP6 gene",
    "Proximal gene expression in monocytes II TRIP6 gene",
    "Gene expression of CUX2 in visual cortex",
    "Gene expression of CUX2 in prefrontal cortex",
    "Colorectal cancer",
    "Colorectal or endometrial cancer",
    "Inflammatory bowel disease",
    "Treatment with folic acid product",
    "Treatment with sulfasalazine",
    "Treatment with methotrexate",
    "Gene expression of BRAP in lymphoblastoid cell lines",
    "Gene expression of ATXN2 in lymphoblastoid cell lines",
    "Gene expression of ALDH2 in lymphoblastoid cell lines",
    "Gene expression of SH2B3 in lymphoblastoid cell lines",
    "Gene expression of ACAD10 in lymphoblastoid cell lines",
    "Self-reported iritis",
    "mDC:%32+; mDC subset (CD32+)",
    "Photic sneeze reflex",
    "Sporadic neuroblastoma",
    "Cystatin c",
    "SH2B3 expression",
    "Self-reported hypothyroidism or myxoedema",
    "Barretts esophagus  or Esophageal adenocarcinoma",
    "Treatment with levothyroxine sodium",
    "Self-reported hyperthyroidism or thyrotoxicosis",
    "Qualifications: college or university degree",
    "Long-standing illness, disability or infirmity",
    "Tonsillectomy",
    "Intestinal malabsorption",
    "Self-reported malabsorption or coeliac disease",
    "Time spent watching television",
    "Axial length",
    "Hand grip strength right",
    "Shingles",
    "Doctor diagnosed sarcoidosis",
    "Hand grip strength left",
    "Ulcerative colitis",
    "Other serious medical condition or disability diagnosed by doctor",
    "Self-reported enlarged prostate",
    "Chickenpox",
    "Proximal gene expression in the liver MUC3A gene",
    "Glaucoma primary open angle",
   "Age at menarche",                                       
    "Leg fat percentage left",                               
    "Retinal venular caliber",                               
   "Hypothyroidism",                                        
    "Soluble intercellular adhesion molecule 1 ICAM 1",      
   "Hair or balding pattern: pattern 3",                    
    "Hair or balding pattern: pattern 4",                    
    "Leg fat percentage right",                              
     "Trunk fat percentage",                                  
     "Menarche age at onset",                                 
    "Menarche",                                              
    "Advanced glycation end product levels",                 
    "Serum magnesium",                                       
     "Magnesium levels",                                      
    "Magnesium",                                             
    "mDC:%32+; mDC subset (CD32+)",                          
     "Infant length",                                         
   "Intracranial volume",                                   
    "Male pattern baldness",                                 
  "Illnesses of father: lung cancer",                      
     "Morning or evening person",                             
   "Mouth or teeth dental problems: dentures",              
   "Self-reported sjogrens syndrome or sicca syndrome",     
     "Population differentiation markers among Great Britain"
    
  )
)


# Function to categorize phenotypes based on key phrases
categorise_phenotypes <- function(phenotypes, key_phrases) {
  categories <- names(key_phrases)
  categorised_phenotypes <- setNames(vector("list", length(categories)), categories)
  
  for (category in categories) {
    phrases <- key_phrases[[category]]
    for (phrase in phrases) {
      matched <- grep(phrase, phenotypes, ignore.case = TRUE, value = TRUE)
      categorised_phenotypes[[category]] <- unique(c(categorised_phenotypes[[category]], matched))
    }
  }
  
  return(categorised_phenotypes)
}

categorised_phenotypes <- categorise_phenotypes(total_unique_traits, key_phrases)
category_counts <- sapply(categorised_phenotypes, length)
total_items <- sum(category_counts)
flattened_categorised <- unlist(categorised_phenotypes)
not_categorised <- setdiff(total_unique_traits, flattened_categorised)

print(not_categorised)


find_category <- function(trait, categorised_phenotypes) {
  for (category in names(categorised_phenotypes)) {
    if (trait %in% categorised_phenotypes[[category]]) {
      return(category)
    }
  }
  return(NA) # Return NA if the trait is not found in any category
}

# Apply the function to each row in the dataframe to annotate it with the corresponding category
lv_df$Trait_category_manually_curated <- sapply(lv_df$Phenoscanner_Trait, find_category, categorised_phenotypes)
#Nas are misc conditions included already but still given NA somehow
lv_df$Trait_category_manually_curated <- ifelse(is.na(lv_df$Trait_category_manually_curated), 
                                                'Miscellaneous_Conditions', 
                                                lv_df$Trait_category_manually_curated)

lv_df_unique <- lv_df %>%
  distinct(across(-Phenotype), .keep_all = TRUE)



rv_df$Trait_category_manually_curated <- sapply(rv_df$Phenoscanner_Trait, find_category, categorised_phenotypes)

# Cleaning NAs that are misc conditions included already but still given NA
rv_df$Trait_category_manually_curated <- ifelse(is.na(rv_df$Trait_category_manually_curated), 
                                                'Miscellaneous_Conditions', 
                                                rv_df$Trait_category_manually_curated)

rv_df_unique <- rv_df %>%
  distinct(across(-Phenotype), .keep_all = TRUE)

table(rv_df_unique$Trait_category_manually_curated)

table(lv_df_unique$Trait_category_manually_curated)

phenotype_order <- c("LVEDV", "LVESV", "LVSV", "LVEF", "LVGFI", "LVMCF", "LVM", "LVMVR", 
                     "RVEDV", "RVESV", "RVSV", "RVEF", "RV_LV_ratio", 
                     "LVEDV_BSA", "LVESV_BSA", "LVSV_BSA", "LVM_BSA", 
                     "RVEDV_BSA", "RVESV_BSA", "RVSV_BSA")

lv_df_unique <- filter(lv_df_unique, !is.na(lv_df_unique$Phenotype))
lv_df_unique[, Phenotype := factor(Phenotype, levels = phenotype_order)]
lv_df_unique <- lv_df_unique[order(Phenotype, CHROM, Lead_SNP)]
lv_df_unique <- filter(lv_df_unique, !is.na(lv_df_unique$rsid))
setnames(lv_df_unique, 'rsid', 'CredibleSet_SNP_rsid')
lv_df_unique <- dplyr::select(lv_df_unique, Phenotype, Locus_name, Locus_number, Lead_SNP, CredibleSet_SNP, CredibleSet_SNP_rsid, Gene_Symbol, 
                              CHROM, GENPOS, BP38, BP37, MAF, BETA, SE, P, Trait_category_manually_curated, Phenoscanner_Trait, 
                              Study, PMID, Phenoscanner_Ancestry, Year, Phenoscanner_A1, Phenoscanner_A2, Phenoscanner_BETA, Phenoscanner_SE,
                              Phenoscanner_P, Phenoscanner_Direction, Phenoscanner_N, Phenoscanner_Unit, Phenoscanner_Dataset)

rv_df_unique <- filter(rv_df_unique, !is.na(rv_df_unique$Phenotype))
rv_df_unique[, Phenotype := factor(Phenotype, levels = phenotype_order)]
rv_df_unique <- rv_df_unique[order(Phenotype, CHROM, Lead_SNP)]
rv_df_unique <- filter(rv_df_unique, !is.na(rv_df_unique$rsid))
setnames(rv_df_unique, 'rsid', 'CredibleSet_SNP_rsid')
rv_df_unique <- dplyr::select(rv_df_unique, Phenotype, Locus_name, Locus_number, Lead_SNP, CredibleSet_SNP, CredibleSet_SNP_rsid, Gene_Symbol, 
                              CHROM, GENPOS, BP38, BP37, MAF, BETA, SE, P, Trait_category_manually_curated, Phenoscanner_Trait, 
                              Study, PMID, Phenoscanner_Ancestry, Year, Phenoscanner_A1, Phenoscanner_A2, Phenoscanner_BETA, Phenoscanner_SE,
                              Phenoscanner_P, Phenoscanner_Direction, Phenoscanner_N, Phenoscanner_Unit, Phenoscanner_Dataset)


fwrite(lv_df_unique, '/Output/Supplementary/supp_table4_lv_all_phenoscanner_cv.csv')
fwrite(rv_df_unique, '/Output/Supplementary/supp_table5_rv_all_phenoscanner_cv.csv')

############################################################################################################
# Counting number of loci and overlaps 

locus_count_lv <- unique(lv_df_unique$Locus_number)
locus_count_rv <- unique(rv_df_unique$Locus_number) 
overlapping_locus_numbers <- intersect(locus_count_lv, locus_count_rv)

# Count the number of overlapping Locus_number values
num_overlapping_locus_numbers <- length(overlapping_locus_numbers) 
trait_count_lv <- unique(lv_df_unique$Phenoscanner_Trait) 
trait_count_rv <- unique(rv_df_unique$Phenoscanner_Trait) 
overlapping_traits <- intersect(trait_count_lv, trait_count_rv)

# Count the number of overlapping Phenoscanner_Trait values
num_overlapping_traits <- length(overlapping_traits) 

all_traits <- rbind(lv_df_unique, rv_df_unique)
all_traits_loci_count <- unique(all_traits$Locus_number)
all_traits_count <- unique(all_traits$Phenoscanner_Trait) 

rv_loci_traits <- rv_df_unique %>%
  group_by(Locus_number) %>%
  summarise(Trait_Categories_Count = n_distinct(Trait_category_manually_curated), 
            Traits = toString(unique(Trait_category_manually_curated))) %>%
  arrange(desc(Trait_Categories_Count))

lv_loci_traits <- lv_df_unique %>%
  group_by(Locus_number) %>%
  summarise(Trait_Categories_Count = n_distinct(Trait_category_manually_curated), 
            Traits = toString(unique(Trait_category_manually_curated))) %>%
  arrange(desc(Trait_Categories_Count)) 

rv_genes_traits <- rv_df_unique %>%
  group_by(Gene_Symbol) %>%
  summarise(Trait_Categories_Count = n_distinct(Trait_category_manually_curated), 
            Traits = toString(unique(Trait_category_manually_curated))) %>%
  arrange(desc(Trait_Categories_Count)) 

lv_genes_traits <- lv_df_unique %>%
  group_by(Gene_Symbol) %>%
  summarise(Trait_Categories_Count = n_distinct(Trait_category_manually_curated), 
            Traits = toString(unique(Trait_category_manually_curated))) %>%
  arrange(desc(Trait_Categories_Count))

############################################################################################################
# Seeing the most pleiotropic genes/SNPs
all_traits_with_counts_snps <- all_traits %>%
  group_by(CredibleSet_SNP_rsid) %>%
  mutate(Traits_Count = n_distinct(Phenoscanner_Trait)) %>%
  ungroup()

all_traits_with_counts_snps <- select(all_traits_with_counts_snps, CredibleSet_SNP_rsid, Gene_Symbol, Traits_Count)
all_traits_with_counts_snps <- unique(all_traits_with_counts_snps) #4796

##########
# Identifying SNPs with CVD traits
table(all_traits$Trait_category_manually_curated)
filtered_traits_cvd <- all_traits %>%
  filter(grepl("Cardiovascular|Blood Pressure|Cardiomyopathy", Trait_category_manually_curated, ignore.case = TRUE))

filtered_traits_cvd_snps <- unique(filtered_traits_cvd$CredibleSet_SNP_rsid) #2069
filtered_traits_cvd_loci <- unique(filtered_traits_cvd$Locus_number) #48


cardiovascular_traits <- c("Atrial fibrillation", "Coronary artery disease", "Dilated cardiomyopathy", "Hypertrophic cardiomyopathy", "Heart failure")

filtered_traits <- all_traits %>%
  filter(Phenoscanner_Trait %in% cardiovascular_traits)
filtered_traits_loci <- unique(filtered_traits$Locus_number)
filtered_traits_snps <- unique(filtered_traits$CredibleSet_SNP_rsid)


filtered_traits_dcm <- all_traits %>%
  filter(grepl("cardiomyopathy", Trait_category_manually_curated, ignore.case = TRUE))

filtered_traits_dcm_loci <- unique(filtered_traits_dcm$Locus_number)
filtered_traits_dcm_names <- unique(filtered_traits_dcm$Locus_name)
filtered_traits_dcm_snps <- unique(filtered_traits_dcm$CredibleSet_SNP_rsid)

filtered_traits_atrial <- all_traits %>%
  filter(grepl("Atrial fibrillation", Phenoscanner_Trait, ignore.case = TRUE))
filtered_traits_atrial_loci <- unique(filtered_traits_atrial$Locus_number)
filtered_traits_atrial_names <- unique(filtered_traits_atrial$Locus_name)
filtered_traits_atrial_snps <- unique(filtered_traits_atrial$CredibleSet_SNP_rsid)
# 7 loci, 10 genes, 14 SNPs
filtered_traits_atrial_list <- dplyr::select(filtered_traits_atrial, Locus_name, Locus_number, CredibleSet_SNP_rsid)
filtered_traits_atrial_list <- unique(filtered_traits_atrial_list)

filtered_traits_cad <- all_traits %>%
  filter(grepl("coronary artery disease", Phenoscanner_Trait, ignore.case = TRUE))
filtered_traits_cad_loci <- unique(filtered_traits_cad$Locus_number)
filtered_traits_cad_names <- unique(filtered_traits_cad$Locus_name)
filtered_traits_cad_snps <- unique(filtered_traits_cad$CredibleSet_SNP_rsid)
# 334 SNPs in 9 loci with 12 genes

atx <- filter(filtered_traits_cad, Gene_Symbol=="ATXN2")
atx <- unique(atx$CredibleSet_SNP_rsid) # 52
atx <- filter(filtered_traits_cad, Gene_Symbol=="SH2B3")
atx <- unique(atx$CredibleSet_SNP_rsid) # 20

filtered_traits_cad_genes <- filter(filtered_traits_cad, Gene_Symbol !="ATXN2"| Gene_Symbol !="SH2B3")


############################################################################################################
# Seeing the genes that overlap with the cardio group and anything else
rv_genes_cardiovascular <- rv_genes_traits %>%
  filter(grepl("Cardiovascular|Blood Pressure|Dilated Cardiomyopathy", Traits) & Trait_Categories_Count > 1) # 104

lv_genes_cardiovascular <- lv_genes_traits %>%
  filter(grepl("Cardiovascular|Blood Pressure|Dilated Cardiomyopathy", Traits) & Trait_Categories_Count > 1) # 64

# Seeing the genes that overlap with the cardio group and anything else except blood and body
rv_genes_cardiovascular_except <- rv_genes_traits %>%
  filter(grepl("Cardiovascular|Blood Pressure|Dilated Cardiomyopathy", Traits) & 
           !grepl("Hematology", Traits) & 
           !grepl("Body Size|Body Composition", Traits) & 
           Trait_Categories_Count > 1) # 21

lv_genes_cardiovascular_except <- lv_genes_traits %>%
  filter(grepl("Cardiovascular|Blood Pressure|Dilated Cardiomyopathy", Traits) & 
           !grepl("Hematology", Traits) & 
           !grepl("Body Size|Body Composition", Traits) & 
           Trait_Categories_Count > 1) # 14

############################################################################################################
# Plotting pleiotropy:
df <- rbind(lv_df_unique, rv_df_unique)

# Count associations per Gene_Symbol and Phenoscanner_Trait
association_counts <- df %>%
  group_by(Gene_Symbol, Phenoscanner_Trait) %>%
  summarise(Count = n_distinct(rsid)) %>%
  ungroup()

association_counts_manual <- df %>%
  group_by(Gene_Symbol, Trait_category_manually_curated) %>%
  summarise(Count = n_distinct(rsid)) %>%
  ungroup()

# Identify top 20 phenotypes and top 20 genes with the most associations
top_phenotypes <- association_counts_manual %>%
  group_by(Trait_category_manually_curated) %>%
  summarise(Total_Associations = sum(Count)) %>%
  top_n(20, Total_Associations) %>%
  pull(Trait_category_manually_curated)

top_genes <- association_counts_manual %>%
  group_by(Gene_Symbol) %>%
  summarise(Total_Associations = sum(Count)) %>%
  top_n(50, Total_Associations) %>%
  pull(Gene_Symbol)

# Filter the dataset for only top genes and phenotypes
filtered_data <- association_counts_manual %>%
  filter(Gene_Symbol %in% top_genes, Trait_category_manually_curated %in% top_phenotypes)

filtered_data <- filter(filtered_data, !is.na(Gene_Symbol))
filtered_data <- filtered_data[!grepl("Miscellaneous", filtered_data$Trait_category_manually_curated), ]


p <- ggplot(filtered_data, aes(x = Gene_Symbol, y = Trait_category_manually_curated, size = Count)) +
  geom_point(aes(size = Count, color = Count), alpha = 0.6) + # Adjust point transparency with alpha
  theme_light(base_size = 14) + # Use theme_light for a clean look with white background
  labs(
    x = "Locus name",
    y = "Phenotypic Category",
    size = "Number of Associations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(breaks = c(1, 10, 50, 100), range = c(1, 10)) +
  scale_color_gradient(low = "darkblue", high = "darkblue") +
  guides(color = FALSE)

# Rotate x-axis labels for better readability

ggsave("/Output/Figure_phenoscanner_cv_curated_traits.png",
       width=20, height=10, plot = p, dpi = 300)
