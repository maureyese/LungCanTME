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

# The purpose of this script is to compare clinical data from our
# patients previously retrieved using 'parte1_datosclinicosLC.R" script.
# We will then get the immune-stromal composition of patients
# using ESTIMATE algorithm and cell composition of TME using
# CIBERSORTx algorithm. We will finally define our design groups.

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
                                        "estimate",
                                        "maftools"), # Add more libraries as needed
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

# To install ESTIMATE algorithm uncomment:
#install.packages("estimate", 
#                 repos="http://r-forge.r-project.org", 
#                 dependencies=TRUE)

#   ____                                            _ _       _           _ 
#  / ___|___  _ __ ___  _ __   __ _ _ __ ___    ___| (_)_ __ (_) ___ __ _| |
# | |   / _ \| '_ ` _ \| '_ \ / _` | '__/ _ \  / __| | | '_ \| |/ __/ _` | |
# | |__| (_) | | | | | | |_) | (_| | | |  __/ | (__| | | | | | | (_| (_| | |
#  \____\___/|_| |_| |_| .__/ \__,_|_|  \___|  \___|_|_|_| |_|_|\___\__,_|_|
#   __| | __ _| |_ __ _|_|                                                  
#  / _` |/ _` | __/ _` |                                                    
# | (_| | (_| | || (_| |                                                    
#  \__,_|\__,_|\__\__,_|                                                    



#  _____ ____ _____ ___ __  __    _  _____ _____ 
# | ____/ ___|_   _|_ _|  \/  |  / \|_   _| ____|
# |  _| \___ \ | |  | || |\/| | / _ \ | | |  _|  
# | |___ ___) || |  | || |  | |/ ___ \| | | |___ 
# |_____|____/ |_| |___|_|  |_/_/   \_\_| |_____|



#   ____ ___ ____  _____ ____  ____   ___  ____ _____      
#  / ___|_ _| __ )| ____|  _ \/ ___| / _ \|  _ \_   _|_  __
# | |    | ||  _ \|  _| | |_) \___ \| | | | |_) || | \ \/ /
# | |___ | || |_) | |___|  _ < ___) | |_| |  _ < | |  >  < 
#  \____|___|____/|_____|_| \_\____/ \___/|_| \_\|_| /_/\_\

# CIBERSORTx is only available via:
# https://cibersortx.stanford.edu/
# We can also install it on our computer, but a lot of
# complex steps are required. We will use their website.

# We will import the results

#  ____            _                                               
# |  _ \  ___  ___(_) __ _ _ __     __ _ _ __ ___  _   _ _ __  ___ 
# | | | |/ _ \/ __| |/ _` | '_ \   / _` | '__/ _ \| | | | '_ \/ __|
# | |_| |  __/\__ \ | (_| | | | | | (_| | | | (_) | |_| | |_) \__ \
# |____/ \___||___/_|\__, |_| |_|  \__, |_|  \___/ \__,_| .__/|___/
#                    |___/         |___/                |_|        



#  ____                                       _                             
# / ___|  __ ___   _____  __      _____  _ __| | _____ _ __   __ _  ___ ___ 
# \___ \ / _` \ \ / / _ \ \ \ /\ / / _ \| '__| |/ / __| '_ \ / _` |/ __/ _ \
#  ___) | (_| |\ V /  __/  \ V  V / (_) | |  |   <\__ \ |_) | (_| | (_|  __/
# |____/ \__,_| \_/ \___|   \_/\_/ \___/|_|  |_|\_\___/ .__/ \__,_|\___\___|
#                                                     |_|                   

#save.image("parte2_agrupacionLC.RData")

# END