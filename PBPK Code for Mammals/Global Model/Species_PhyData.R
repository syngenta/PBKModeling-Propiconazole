# load library
library(stationaRy)
library(isdparser)
library(tidyr)
library(lubridate)
library(WriteXLS)
library(readxl)
library(writexl)
library("tidyr")
library(dplyr)
library(Metrics)
library(cowplot)
library(knitr)
library(reshape2)
library(httk)
library(deSolve)
library(pksensi)
library(ggplot2)


#=============================================================================================================
#                                   Generic Maternal PBPK Model Parameters                                   #
#                              Tissue and physiological data for mammal species                              #
#=============================================================================================================
##############################################################################################################
#species                 <- 'Mouse'              #  "Rat", "Rabbit", "Dog", "Mouse", or "Human"
pH                      <- 7.45
R_rat                   <- 0.11                  # 0.2cm # radius of rat jejunum; https://github.com/wfsrqivive/rat_PBK/blob/main/model/generate_input_data.r
R_human                 <- 1.25                  # cm; The small intestine, which is 670 to 760 cm (22 to 25 feet) in length and 3 to 4 cm (about 2 inches) in diameter, https://www.britannica.com/science/human-digestive-system/Anatomy
##############################################################################################################
# Species-specific physiology parameters
df_physio  <- data.frame(physiology.data$Parameter,
                         physiology.data$Units,
                         physiology.data[,grepl(species, names(physiology.data))]) 

# Tissue composition and species-specific physiology parameters
df_tissue  <- tissue.data[which(tissue.data$Species == species),]

# existing database: chem <- chem.physical_and_invitro.data; kg
BW         <- subset(df_physio, physiology.data.Parameter == 'Average BW')$"physiology.data...grepl.species..names.physiology.data..."   # kg

############################################################
###               tissue volumes  (L/kg BW)              ###  
############################################################
# 12 for rat, 14 for human; subset(df_tissue , variable == "Vol (L/kg)"); placenta not included for rat
Vadipose   <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'adipose')$value
Vbone      <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'bone')$value
Vbrain     <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'brain')$value
Vgut       <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'gut')$value
Vheart     <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'heart')$value
Vkidney    <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'kidney')$value
Vliver     <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'liver')$value
Vlung      <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'lung')$value
Vmuscle    <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'muscle')$value
Vskin      <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'skin')$value
Vspleen    <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'spleen')$value
#thyroid   <- subset(df_tissue , variable == "Vol (L/kg)" & 
#                      tolower(Tissue) == 'thyroid')$value      # human only
#placenta                                                       # human only
Vrest_httk <- subset(df_tissue , variable == "Vol (L/kg)" & 
                       tolower(Tissue) == 'rest')$value

Vrest      <- Vbone  + Vheart + Vskin + Vspleen + Vmuscle + Vkidney + Vrest_httk    # L/kg BW
# quality check: 
Vol        <- subset(df_tissue , variable == "Vol (L/kg)")
Vrest+Vlung+Vgut+Vliver+Vkidney + Vbrain + Vadipose  == sum(Vol$value)

Hematocrit_fraction <- subset(df_physio, physiology.data.Parameter == 'Hematocrit')$"physiology.data...grepl.species..names.physiology.data..."
Vplasma    <- subset(df_physio, physiology.data.Parameter == 'Plasma Volume')$"physiology.data...grepl.species..names.physiology.data..."/1000   # L/kg

Vart       <- Vplasma / (1 - Hematocrit_fraction) * 0.25             # L/kg 72/91 Sitovition report
Vven       <- Vplasma / (1 - Hematocrit_fraction) * 0.75


##########################################################
###              tissue flows (L/h/kg BW)              ###  
##########################################################
# 12 for rat, 14 for human
Qadipose   <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'adipose')$value /(BW^0.25)/1000*60 
Qbone      <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'bone')$value /(BW^0.25)/1000*60 
Qbrain     <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'brain')$value /(BW^0.25)/1000*60 
Qgut       <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'gut')$value /(BW^0.25)/1000*60 
Qheart     <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'heart')$value /(BW^0.25)/1000*60 
Qkidney    <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'kidney')$value /(BW^0.25)/1000*60 
Qliver     <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'liver')$value /(BW^0.25)/1000*60 
Qlung      <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'lung')$value /(BW^0.25)/1000*60 
Qmuscle    <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'muscle')$value /(BW^0.25)/1000*60 
Qskin      <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'skin')$value /(BW^0.25)/1000*60 
Qspleen    <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'spleen')$value /(BW^0.25)/1000*60 
Qrest_httk <- subset(df_tissue , variable == "Flow (mL/min/kg^(3/4))" & 
                       tolower(Tissue) == 'rest')$value /(BW^0.25)/1000*60 

Qcardiac   <- subset(df_physio, physiology.data.Parameter == 'Cardiac Output')$"physiology.data...grepl.species..names.physiology.data..." /(BW^0.25)/1000*60                  # L/h/kg BW
Qrest      <- Qcardiac - Qliver - Qgut - Qadipose - Qbrain 


#----------------------------                   End of general species parameters         -----------------------------------

