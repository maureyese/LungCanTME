x = 2

#x es igual a 2
setwd("C:/Users/xavie/Documents/LungCanTME")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("curatedTCGAData")

#abrir la libreria de TCGA utils
library("TCGAutils")
library("curatedTCGAData")
luad <- curatedTCGAData(
  diseaseCode = "LUAD", assays = c("CNASeq", "Mutation", "miRNA*",
                                   "RNASeq2*", "mRNAArray", "Methyl*"), version = "1.1.38", dry.run = FALSE
)
#revisar la cant de pacientes por datos moleculares
sampleTables(luad)
install.packages(tidyverse)
BiocManager::install("TCGAbiolinks")
BiocManager::install("TCGAbiolinks")



