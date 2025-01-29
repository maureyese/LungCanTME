# ===============================
# Authors: Ximena Vázquez-Cadena, Mauricio Reyes-Elizondo
#
# Project: An mRNA–miRNA Regulatory Network to Identify Key microRNAs 
#          Driving Immune Infiltration and Prognosis in Lung Cancer
#
# Related peers:
#   Dr. Vianey González Villasana (FCB, UANL)
#   Dr. Diana Reséndez Pérez (FCB, UANL)
# ==============================

# The purpose of this script is to analyze number of available
# samples for our analysis. This step is crucial to know if
# we have sufficient samples to start our research.
# In this case, we will be analyzing TCGA-LUAD and TCGA-LUSC
# projects. We will also download RNAseq, miRNAseq and clinical
# data.

# __        __         _    _                
# \ \      / /__  _ __| | _(_)_ __   __ _    
#  \ \ /\ / / _ \| '__| |/ / | '_ \ / _` |   
#   \ V  V / (_) | |  |   <| | | | | (_| |   
#  __\_/\_/ \___/|_|  |_|\_\_|_| |_|\__, |   
# |  _ \(_)_ __ ___  ___| |_ ___  _ |___/  _ 
# | | | | | '__/ _ \/ __| __/ _ \| '__| | | |
# | |_| | | | |  __/ (__| || (_) | |  | |_| |
# |____/|_|_|  \___|\___|\__\___/|_|   \__, |
#                                      |___/ 

# This is the folder we will be storing our files
setwd("C:/Users/xavie/Documents/LungCanTME")

# Do not uncomment this code
#setwd("~/Research/LungCanTME")

#  _     _ _                    _           
# | |   (_) |__  _ __ __ _ _ __(_) ___  ___ 
# | |   | | '_ \| '__/ _` | '__| |/ _ \/ __|
# | |___| | |_) | | | (_| | |  | |  __/\__ \
# |_____|_|_.__/|_|  \__,_|_|  |_|\___||___/

# Load libraries
suppressPackageStartupMessages(lapply(c("tidyverse",
                                        "TCGAutils", 
                                        "curatedTCGAData", 
                                        "TCGAbiolinks",
                                        "HDF5Array"), # Add more libraries as needed
                                      library, 
                                      character.only = TRUE))

# If a library is not install on your computer
# Check if the library is from CRAN or Bioconductor
# repository by searching on Google.

# Once you identified the repository you can use:
# For CRAN:
# install.packages()
# For Bioconductor:
# BiocManager::install()

#   ____      _                 _   _            _       
#  / ___| ___| |_   _ __   __ _| |_(_) ___ _ __ | |_ ___ 
# | |  _ / _ \ __| | '_ \ / _` | __| |/ _ \ '_ \| __/ __|
# | |_| |  __/ |_  | |_) | (_| | |_| |  __/ | | | |_\__ \
#  \____|\___|\__| | .__/ \__,_|\__|_|\___|_| |_|\__|___/
#                  |_|                                   

# Lung Adenocarcinoma (LUAD) =======================================

# Get general information from the TCGA-LUAD project
luad <- curatedTCGAData(
  diseaseCode = "LUAD", # Name of project
  assays = c("CNASeq", "Mutation", "miRNA*",
             "RNASeq2*", "mRNAArray", "Methyl*"), # Name of molecular data 
  version = "1.1.38", dry.run = FALSE
) # Retrieve a S4 object (Table of tables)

# Review available molecular data per patient
sampleTables(luad)

# Retrieve important table from luad (the table of tables)
luad_samples <- as.data.frame(luad@sampleMap)
# Divide samples in tumor (01), healthy (11) or any other
luad_samples$sample_types <- substring(luad_samples$colname, 14, 15)
# Check which type of samples we have
unique(luad_samples$sample_types)
# Since we have more samples rather than just "01" of "11"
# We will remove them
luad_samples <- luad_samples[!(luad_samples$sample_types == "02"),] # Remove unnecessary samples
luad_samples <- luad_samples[!(luad_samples$sample_types == "10"),] # Remove unnecessary samples
# Check again type of samples
unique(luad_samples$sample_types)
# We only have now "01" and "11" samples

# You can also mimic "sampleTables" function 
# to count available samples per patient
table(luad_samples$assay,luad_samples$sample_types)
# You get the same information as sampleTables(luad)
# But using a cleaned table

# We will now get the patient IDs and classify them as:
#     1. "Cancer" or "Normal"
#     2. "RNASeq" or "miRNAseq" (miRNAs)

# We will use %>% function to manipulate the luad_samples easily readable
# %>% function comes from tidyverse library (dplyr, in particular)

# Retrieve cancer samples with available transcriptomic data
luad_cancer_rna_codes <- luad_samples %>% # Get the table
  filter(assay == "LUAD_RNASeq2GeneNorm-20160128", sample_types == "01") %>% # Filter columns
  pull(primary) # Get the IDs
# Retrieve cancer samples with available miRNA data
luad_cancer_mirna_codes <- luad_samples %>% 
  filter(assay == "LUAD_miRNASeqGene-20160128", sample_types == "01") %>%
  pull(primary)
# Retrieve normal samples with available transcriptomic data
luad_normal_rna_codes <- luad_samples %>% 
  filter(assay == "LUAD_RNASeq2GeneNorm-20160128", sample_types == "11") %>%
  pull(primary)
# Retrieve normal samples with available miRNA data
luad_normal_mirna_codes <- luad_samples %>% 
  filter(assay == "LUAD_miRNASeqGene-20160128", sample_types == "11") %>%
  pull(primary)

# Remove patients that have miRNA data but no RNAseq data
luad_cancer_mirna_codes <- luad_cancer_mirna_codes[luad_cancer_mirna_codes %in% luad_cancer_rna_codes]
length(luad_cancer_mirna_codes) # 342 patients with both RNAseq and miRNA
# Remove normal samples that are not related with IDs from RNAseq
luad_normal_rna_codes <- luad_normal_rna_codes[luad_normal_rna_codes %in% luad_cancer_rna_codes]
length(luad_normal_rna_codes) # 58 samples out of 59
# Remove normal samples that are not related with IDs from miRNA
luad_normal_mirna_codes <- luad_normal_mirna_codes[luad_normal_mirna_codes %in% luad_cancer_mirna_codes]
length(luad_normal_mirna_codes) # 38 samples out of 46

# Lung Squamuous Cell Carcinoma (LUSC) =======================================

# Get general information from the TCGA-LUSC project
lusc <- curatedTCGAData(
  diseaseCode = "LUSC", # Name of project
  assays = c("CNASeq", "Mutation", "miRNA*",
             "RNASeq2*", "mRNAArray", "Methyl*"), # Name of molecular data
  version = "1.1.38", dry.run = FALSE
) # Retrieve a S4 object (Table of tables)

# Review available molecular data per patient
sampleTables(lusc)

# Retrieve important table from lusc (the table of tables)
lusc_samples <- as.data.frame(lusc@sampleMap)
# Divide samples in tumor (01), healthy (11) or any other
lusc_samples$sample_types <- substring(lusc_samples$colname, 14, 15)
# Check which type of samples we have
unique(lusc_samples$sample_types)
# Since we only have "01" and "11"
# We skip filtering our table

# We will now proceed to classify patients as:
#     1. "Cancer" or "Normal"
#     2. "RNASeq" or "miRNAseq" (miRNAs)

# Retrieve cancer samples with available transcriptomic data
lusc_cancer_rna_codes <- lusc_samples %>% 
  filter(assay == "LUSC_RNASeq2GeneNorm-20160128", sample_types == "01") %>%
  pull(primary)
# Retrieve cancer samples with available miRNA data
lusc_cancer_mirna_codes <- lusc_samples %>% 
  filter(assay == "LUSC_miRNASeqGene-20160128", sample_types == "01") %>%
  pull(primary)
# Retrieve normal samples with available transcriptomic data
lusc_normal_rna_codes <- lusc_samples %>% 
  filter(assay == "LUSC_RNASeq2GeneNorm-20160128", sample_types == "11") %>%
  pull(primary)
# Retrieve normal samples with available miRNA data
lusc_normal_mirna_codes <- lusc_samples %>% 
  filter(assay == "LUSC_miRNASeqGene-20160128", sample_types == "11") %>%
  pull(primary)

# Remove patients that have miRNA data but no RNAseq data
lusc_cancer_mirna_codes <- lusc_cancer_mirna_codes[lusc_cancer_mirna_codes %in% lusc_cancer_rna_codes]
length(lusc_cancer_mirna_codes) # 342 patients with both RNAseq and miRNA
# Remove normal samples that are not related with IDs from RNAseq
lusc_normal_rna_codes <- lusc_normal_rna_codes[lusc_normal_rna_codes %in% lusc_cancer_rna_codes]
length(lusc_normal_rna_codes) # 51 samples out of 51
lusc_normal_mirna_codes <- lusc_normal_mirna_codes[lusc_normal_mirna_codes %in% lusc_cancer_mirna_codes]
length(lusc_normal_mirna_codes) # 45 samples out of 45

#   ____      _                     _                 _            
#  / ___| ___| |_   _ __ ___   ___ | | ___  ___ _   _| | __ _ _ __ 
# | |  _ / _ \ __| | '_ ` _ \ / _ \| |/ _ \/ __| | | | |/ _` | '__|
# | |_| |  __/ |_  | | | | | | (_) | |  __/ (__| |_| | | (_| | |   
#  \____|\___|\__| |_| |_| |_|\___/|_|\___|\___|\__,_|_|\__,_|_|   
#   __| | __ _| |_ __ _                                            
#  / _` |/ _` | __/ _` |                                           
# | (_| | (_| | || (_| |                                           
#  \__,_|\__,_|\__\__,_|                                           

# Lung Adenocarcinoma (LUAD) =======================================

# First we will retrieve RNA-Seq data from LUAD

# Start a query to access to GDC database
query <- GDCquery(project = "TCGA-LUAD",
                         data.category = "Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification", 
                         experimental.strategy = "RNA-Seq",
                         workflow.type = "STAR - Counts",
                         barcode = c(luad_cancer_rna_codes,luad_normal_rna_codes)) # Patients

# Download once, make sure to work on the correct directory
# Uncomment if you do not have the files on your folder
#GDCdownload(query)

# Prepare the RNA-Seq from the downloaded files
raw.counts_luad <- GDCprepare(query = query)
# Get the expression table
rnacounts_luad <- raw.counts_luad@assays@data@listData$unstranded
# Label columns using patient IDs
colnames(rnacounts_luad) <- colnames(raw.counts_luad)
# Label rows using ENSEMBL IDs
# Split human genome version from ENSEMBL IDs
ensg <- unlist(strsplit(rownames(raw.counts_luad), split = "[.]"))
ensg <- ensg[c(TRUE, FALSE)] # Remove human genome version
rownames(rnacounts_luad) <- ensg # Add ENSEMBL to row names

# The expression table from RNAseq is now ready to use

# Save non-normalized expression table
saveRDS(rnacounts_luad, "~/Research/LungCanTME/ExpressionData/rnacounts_luad.RDS") # Store as R object
write.csv(rnacounts_luad, "~/Research/LungCanTME/ExpressionData/counts_luad.csv", row.names = TRUE) # Store as table
# csv files can be opened on Excel

# We will now retrieve miRNASeq data from LUAD

query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling", 
                  data.type = "miRNA Expression Quantification", 
                  experimental.strategy = "miRNA-Seq",
                  barcode = c(luad_cancer_mirna_codes,luad_normal_mirna_codes)) # Patients

# Download once, make sure to work on the correct directory
# Uncomment if you do not have the files on your folder
#GDCdownload(query)

# Prepare the miRNA-Seq from the downloaded files
raw.mircounts_luad <- GDCprepare(query = query)
# A table was downloaded
# First column of table is miRNA names
# Remaining columns are miRNA expression counts and
# Other columns not needed

# Create a new table and extract only
# the miRNA expression count columns that start with "read_count"
mircounts_luad <- raw.mircounts_luad[, grep("^read_count", colnames(raw.mircounts_luad))]
# Remove "read_count" from the column names
colnames(mircounts_luad) <- gsub("^read_count_", "", colnames(mircounts_luad))
# Add miRNA names to row names
rownames(mircounts_luad) <- raw.mircounts_luad[,1]

# The expression table is ready to use

# Save non-normalized expression table
saveRDS(mircounts_luad, "~/Research/LungCanTME/ExpressionData/mircounts_luad.RDS") # Store as R object
write.csv(mircounts_luad, "~/Research/LungCanTME/ExpressionData/mircounts_luad.csv", row.names = TRUE) # Store as table
# csv files can be opened on Excel

# Lung Squamous Cell Carcinoma (LUSC) =======================================

# First we will retrieve RNA-Seq data from LUAD

# Start a query to access to GDC database
query <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts",
                  barcode = c(lusc_cancer_rna_codes,lusc_normal_rna_codes)) # Patients

# Download once, make sure to work on the correct directory
# Uncomment if you do not have the files on your folder
#GDCdownload(query)

# Prepare the RNA-Seq from the downloaded files
raw.counts_lusc <- GDCprepare(query = query)
# Get the expression table
rnacounts_lusc <- raw.counts_lusc@assays@data@listData$unstranded
# Label columns using patient IDs
colnames(rnacounts_lusc) <- colnames(raw.counts_lusc)
# Label rows using ENSEMBL IDs
# Split human genome version from ENSEMBL IDs
ensg <- unlist(strsplit(rownames(raw.counts_lusc), split = "[.]"))
ensg <- ensg[c(TRUE, FALSE)] # Remove human genome version
rownames(rnacounts_lusc) <- ensg # Add ENSEMBL to row names

# The expression table from RNAseq is now ready to use

# Save non-normalized expression table
saveRDS(rnacounts_lusc, "~/Research/LungCanTME/ExpressionData/rnacounts_lusc.RDS") # Store as R object
write.csv(rnacounts_lusc, "~/Research/LungCanTME/ExpressionData/counts_lusc.csv", row.names = TRUE) # Store as table
# csv files can be opened on Excel

# We will now retrieve miRNASeq data from LUAD

query <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Transcriptome Profiling", 
                  data.type = "miRNA Expression Quantification", 
                  experimental.strategy = "miRNA-Seq",
                  barcode = c(lusc_cancer_mirna_codes,lusc_normal_mirna_codes)) # Patients

# Download once, make sure to work on the correct directory
# Uncomment if you do not have the files on your folder
#GDCdownload(query)

# Prepare the miRNA-Seq from the downloaded files
raw.mircounts_lusc <- GDCprepare(query = query)
# A table was downloaded
# First column of table is miRNA names
# Remaining columns are miRNA expression counts and
# Other columns not needed

# Create a new table and extract only
# the miRNA expression count columns that start with "read_count"
mircounts_lusc <- raw.mircounts_lusc[, grep("^read_count", colnames(raw.mircounts_lusc))]
# Remove "read_count" from the column names
colnames(mircounts_lusc) <- gsub("^read_count_", "", colnames(mircounts_lusc))
# Add miRNA names to row names
rownames(mircounts_lusc) <- raw.mircounts_lusc[,1]

# The expression table is ready to use

# Save non-normalized expression table
saveRDS(mircounts_lusc, "~/Research/LungCanTME/ExpressionData/mircounts_lusc.RDS") # Store as R object
write.csv(mircounts_lusc, "~/Research/LungCanTME/ExpressionData/mircounts_lusc.csv", row.names = TRUE) # Store as table
# csv files can be opened on Excel

#   ____      _          _ _       _           _       _       _        
#  / ___| ___| |_    ___| (_)_ __ (_) ___ __ _| |   __| | __ _| |_ __ _ 
# | |  _ / _ \ __|  / __| | | '_ \| |/ __/ _` | |  / _` |/ _` | __/ _` |
# | |_| |  __/ |_  | (__| | | | | | | (_| (_| | | | (_| | (_| | || (_| |
#  \____|\___|\__|  \___|_|_|_| |_|_|\___\__,_|_|  \__,_|\__,_|\__\__,_|

# Lung Adenocarcinoma (LUAD) =======================================

# We will now download clinical data

# Start a query to access to GDC database
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Clinical",
  data.format = "bcr xml",
  barcode = luad_cancer_rna_codes
)

# Download once, make sure to work on the correct directory
# Uncomment if you do not have the files on your folder
GDCdownload(query)

# We will prepare the data
# Prepare general data
luad_clinical <- GDCprepare_clinic(query, clinical.info = "patient")
# Export data
saveRDS(luad_clinical, "~/Research/LungCanTME/ClinicalData/luad_clinical.RDS")
write.csv(luad_clinical, "~/Research/LungCanTME/ClinicalData/luad_clinical.csv", row.names = FALSE)

# Prepare follow-up data
luad_followup <- GDCprepare_clinic(query, clinical.info = "follow_up")
# Export data
saveRDS(luad_followup, "~/Research/LungCanTME/ClinicalData/luad_followup.RDS")
write.csv(luad_followup, "~/Research/LungCanTME/ClinicalData/luad_followup.csv", row.names = FALSE)

# Prepare drug data
luad_drug <- GDCprepare_clinic(query, clinical.info = "drug")
# Export data
saveRDS(luad_drug, "~/Research/LungCanTME/ClinicalData/luad_drug.RDS")
write.csv(luad_drug, "~/Research/LungCanTME/ClinicalData/luad_drug.csv", row.names = FALSE)

# Prepare admin

# Prepare radiation

# Prepare stage_event

# Prepare new_tumor_event

# Lung Squamous Cell Carcinoma (LUSC) =======================================

#  ____                                       _                             
# / ___|  __ ___   _____  __      _____  _ __| | _____ _ __   __ _  ___ ___ 
# \___ \ / _` \ \ / / _ \ \ \ /\ / / _ \| '__| |/ / __| '_ \ / _` |/ __/ _ \
#  ___) | (_| |\ V /  __/  \ V  V / (_) | |  |   <\__ \ |_) | (_| | (_|  __/
# |____/ \__,_| \_/ \___|   \_/\_/ \___/|_|  |_|\_\___/ .__/ \__,_|\___\___|
#                                                     |_|                   

# We will be working now on the R script "parte2_agrupacionLC.R"
# to analyze clinical data and get our design groups based on
# ESTIMATE and CIBERSORTx algorithms

save.image("parte1_datosclinicosLC.RData")

# END