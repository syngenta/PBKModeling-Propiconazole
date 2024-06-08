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
library(PKNCA)
library(cowplot)
library(knitr)
library(reshape2)
library(httk)
library(deSolve)
library(pksensi)
library(ggplot2)
# make sure HHTK is the latest version
#update.packages('httk') 
packageVersion("httk")
# rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)

#==========================================================================================
#                                         Start                                           #
#==========================================================================================
# Reference
# PKNCA : https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk  : https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

# - Target compound : PPZ (parent)
# - Species         : Human (Maternal only)
# - Author          : Yaoxing Wu
# - Date            : Jun, 2022
# - Structure       : GI tract, Liver, Lungs, Adipose, Brain, Rest, Vein, Arterial

# Clean the environment
rm(list = ls())
#dev.off()

species                 <- 'Human'  

source(paste("C:/Users/s1036120/OneDrive - Syngenta/HTTK/Mammals/General_code/Maternal_8compt_InductionModel_Parent.R"))
source(paste("C:/Users/s1036120/OneDrive - Syngenta/HTTK/Mammals/General_code/Species_PhyData.R"))

MPPGL      <- MPPGL_human
MPPGGI     <- MPPGGI_human
kdeg       <- kdeg_human

BW_mouse    <- 0.02
BW_rat      <- 0.25
BW_human    <- 70

Qrest       ==  Qcardiac - Qliver - Qgut - Qadipose - Qbrain 


#======================================================================================================
#                          Compound physico-chemical & in vitro properties                            #
#======================================================================================================
####################### Parent (Table 1)    Acibenzolar ######################
# physico-chemical properties   
Compound_name_parent    <- 'Propiconazole (PPZ)'
Compound_type_parent    <- 'base'
CAS_parent              <- "60207-90-1-a"                                # add a to CAS since PPZ has been included in HTTK database
LogP_parent             <-  3.72                                         # log Kow; from raav
MW_parent               <-  342.2
pKa_a_parent            <-  'None'                                       # very weak base without any numerical value; chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://hal-unilim.archives-ouvertes.fr/hal-00929974/document
pKa_b_parent            <-  1.09                                         # pka accept (base); Teb?H+ -> Teb + H+ 

# in vitro properties   -
fub_parent              <-  0.05 #0.0955                                         # Unbound fraction in plasma; measured in sed
Rblood2plasma_parent    <-  0.597                                          #  if not available, place "None" here.
Papp_parent             <-  40.8                                           # unit: 10-6 cm/s                 
#a <- chem.physical_and_invitro.data[chem.physical_and_invitro.data$CAS == '60207-90-1',]

BW_scaledfrom                    <- BW_rat
# in vitro metabolic-related properties   
Vmax_unit                        <- 'umol/h/kg bw'                         # 'umol/h/kg bw' or 'none'
# liver
type_hep_clearance               <- 'liver'
incubation_hep_parent            <- 'microsome'                            # for rat metabolic stability test (hepatic)
# 13 L/h/kg Bw (216 mL/min/kg)
Km_hep_parent                    <- 4.6                                    # unit: uM (umol/L)
Vmax_hep_parent                  <- 15.5 * Km_hep_parent                   # umol/h/kg bw (59.8)


df_tissue_human                     <- tissue.data[which(tissue.data$Species == 'Human'),]
Vliver_human                        <- subset(df_tissue_human , variable == "Vol (L/kg)" &
                                             tolower(Tissue) == 'liver')$value

parent_hep_conc_metabolic        <- 0.5                                   # 0.5 mg/ml protein; Shen et al. 2013
fuinc_hep_parent                 <- 9999                                  # binding was not measured in in vitro metabolsim assays
parent_hep_conc_binding          <- 9999
Clint_ori_hep_parent             <- 9999
# intestine
type_int_clearance               <- 'intestine'
incubation_int_parent            <- 'microsome'                           # for rat intestinal stability test
Vmax_int_parent                  <- Vmax_hep_parent  / MPPGL_human  * MPPGGI_human          # unit: nmol/min/mg protein; Table 1 -->   umol/h/mg protein; Intestinal metabolism determined at 50 ug/mL Frost 2017
Km_int_parent                    <- Km_hep_parent                         # unit: uM (umol/L)
parent_int_mic_conc_metabolic    <- 9999                                  # Intestinal metabolism determined at 50 ug/mL Frost 2017
fuinc_int_parent                 <- 9999                                  # Measured at 1 ug/mL microsome binding experiment; Frost 2017 (16/90)
parent_int_mic_conc_binding      <- 9999  
# plasma
type_plasma_clearance            <- 'plasma'
incubation_plasma_parent         <- 'half-life'   
Clint_plasma_parent              <-  1e6                                  # half life (h)
fu_plasma_parent                 <- 'None'

# Induction related parapeters
EC50         <- 6.26                   # uM mouse
Emax         <- 5.62                   # fold
E0           <- 1 
fmCYP        <- 0.65                   # CYP3A4

# absorption extent
fa           <- 1                                                         # absorption fraction

#calc_fup_correction(fup = 0.0721, chem.cas = '107534-96-3-a',species = "Rabbit")
###################################################################################
#                     Input to HTTK for parameterization                          #
###################################################################################
# Reference: add_chemtable (Httk manual dated Sep. 22, 2022, page 7)
# Add acibenzolar and acibenzolar acid data to HTTK chemical table for further compound-specific parameterization using HTTK
my.new.data           <- as.data.frame(cbind(Compound_name_parent, CAS_parent, "DTX_parent", MW_parent, LogP_parent, fub_parent, 0, pKa_a_parent, pKa_b_parent))
colnames(my.new.data) <- c("Name","CASRN","DTXSID","MW","LogP","Fup","CLint", "pKa.a", 'pKa.b')
chem.physical_and_invitro.data <- add_chemtable(my.new.data,
                                                current.table=
                                                  chem.physical_and_invitro.data,
                                                data.list=list(
                                                  Compound  ="Name",
                                                  CAS       ="CASRN",
                                                  DTXSID    ="DTXSID",
                                                  MW        ="MW",
                                                  logP      ="LogP",
                                                  Funbound.plasma="Fup",
                                                  Clint     ="CLint",
                                                  pKa_Donor = 'pKa.a',
                                                  pKa_Accept= 'pKa.b'),
                                                species     = species,
                                                reference   ="None")


#=========================================================================================
#                    Chemical-specific Maternal PBPK Model Parameters                    #
#=========================================================================================
####################################################################### 
###                    partition coefficients                       ###  to unbound plasma; 12 for rat, 14 for human
#######################################################################
source(paste("C:/Users/s1036120/OneDrive - Syngenta/HTTK/Mammals/General_code/Maternal_DM_Parent.R"))
# Validation only
tissuelist <- list(liver=c("liver"),lung=c("lung"),gut=c("gut"), brain=c("brain"), adipose=c("adipose"),
                    muscle.bone=c('bone', 'heart', 'skin', 'spleen', 'muscle', 'rest', "kidney"))
round(lump_tissues(PC_parent,tissuelist=tissuelist,species = species)$Kmuscle.bone2pu,1)  == round(Krest2pu_parent, 1) 
 
##################################################################
###                     absorption (1/h)                       ###  
##################################################################
###                 Acibenzolar                   ###
# calculation of Ka using caco2 permeability was not used as experimentally fitting values are available
# Ka was obtained by fitting to the observed blood concentration-time data for acibenzolar-acid after oral dosing (1, 10 and 100 mg/kg acibenzolar)
# as acibenzolar was undetectable in blood at these doses.
# Reference: Acibenzolar manuscript equation 1
calc_ka     <- function(dose.mg.kg, Papp){
  if(is.numeric(Papp)){
    Peff_human  <- 10^(0.4926*log10(Papp) - 0.1454)                  # 1e-4 cm/s; Papp in the unit of 1e-6 cm/s
    ka          <- 2 * Peff_human / R_human / 1e4 * 3600
  }else{
    ka_scaledfrom   <- 2
    ka              <- ka_scaledfrom  * (BW/BW_scaledfrom)^(-0.25)                       # eqaution A11 in appendix
  }
  
  return(ka)                                                # unit: 1/h
}

##########################################################################
###                        excretion L/h/kg                            ###  
##########################################################################
#  renal clearance 
CLR_rat              <- 'None'     # L/h/kg no urinary excretion for teb parent compound.

calc_GFR         <- function(CLR){
  if(CLR == "None"){
    Qgfr       <- subset(df_physio, physiology.data.Parameter == 'GFR')$"physiology.data...grepl.species..names.physiology.data..."/(BW^0.25)/1000*60     # L/h/kg
  }else{
    Qgfr       <- (CLR * BW_rat) * (BW / BW_rat)^0.75 / BW   
  }
  
  return(Qgfr)
}
  
CLR  <- calc_GFR(CLR = CLR_rat)
#----------------------------                   End of compound specific parameters         -----------------------------------  



#================================================================================================================
###                                      Define modeling variables                                            ###  
#================================================================================================================
# in vivo experiment
days          <- 40
Time_max      <- days * 24                             # h
StartTime     <- 0                                     # Time_min 
StopTime      <- Time_max
dt            <- 0.1
Times         <- seq(StartTime, StopTime, dt)

#           Dose regimen            #
# oral dose 
oral_mg.kg        <- 1#400#150#29#400#150#24.7#400#150#24.7                  # mg/kg
oral_input        <- oral_mg.kg  * 1000 / MW_parent   # mg/kg -> umol/kg bw

# iv dose 
iv_mg.kg          <- 1                               # 
iv_input          <- iv_mg.kg  * 1000  / MW_parent   # mg/kg -> umol/kg bw

# infusion dose
Infusion_mg.kg    <- 0                                            # 1, 10 or 100; 2017 in vivo rat study
Infusion_input    <- Infusion_mg.kg  * 1000 / MW_parent           # mg/kg -> ug/kg bw
CTinf             <- 1                                            # 1 h infusion time
CRATE             <- Infusion_input/CTinf                         # ug/kg bw /h
CTIMEinf          <- c(0,CTinf,Time_max)
CRATEinf          <- c(CRATE,0,0)
# Define an interpolation function that returns rate when given time - "const"
Cstep.doseinf     <- approxfun(CTIMEinf, CRATEinf, method = "const")

# absorption calculation
ka                <- 1                                             # h-1
kt                <- kt_human                                                 

ndays         <- 2    #13#0#2#13
ndays_single  <- 13   #2#13
Dose_events_oral   <- data.frame(var = "Agutlumen_parent",           # Agutlumen_parent; Aven_parent
                                 time = seq(from = 0, to = ndays * 24, by = 24),
                                 value = rep(oral_input,ndays+1),
                                 method = "add")

Dose_events_iv      <- data.frame(var = "Aven_parent",              # Agutlumen_parent; Aven_parent
                                  time = 0,
                                  value = iv_input,
                                  method = "add")


Dose_events_oral_Single      <- data.frame(var = "Agutlumen_parent",           # Agutlumen_parent; Aven_parent
                                  time = seq(from = 0, to = ndays_single * 24, by = 24),
                                  value = rep(oral_input,ndays_single+1),
                                  method = "add")

#########################    parameters for exposure scenario    #########################
parms <- c( Oral_mg.kg                = oral_mg.kg, 
            Oral_input                = oral_input,   
            
            pH                        = pH,
            fa                        = fa,             
            
            Rblood2plasma_parent      = Rblood2plasma_parent,
            BW                        = BW,              
            
            fub_parent                = fub_parent,   
            ka                        = ka,
            kt                        = kt,
            
            EC50                      = EC50,              # uM
            Emax                      = Emax,                  # fold
            E0                        = E0, 
            kdeg                      = kdeg,
            fmCYP                     = fmCYP,
            # partiton coefficient
            Kgut2pu_parent            = Kgut2pu_parent,      
            Kliver2pu_parent          = Kliver2pu_parent,     
            Klung2pu_parent           = Klung2pu_parent,  
            Kbrain2pu_parent          = Kbrain2pu_parent,
            Kadipose2pu_parent        = Kadipose2pu_parent,
            Krest2pu_parent           = Krest2pu_parent ,     
 
            
            #incubation_hep_parent     = incubation_hep_parent,                                    
            Vmax_hep_parent            = Vmax_hep_parent,                                      
            Km_hep_parent              = Km_hep_parent,                                         
            Vmax_int_parent            = Vmax_int_parent ,                                           
            Km_int_parent              = Km_int_parent,                                             # unit: uM (umol/L)
            parent_hep_conc_metabolic  = parent_hep_conc_metabolic,                            
            fuinc_hep_parent           = fuinc_hep_parent,                              
            parent_hep_conc_binding    = parent_hep_conc_binding,   
            parent_int_mic_conc_metabolic   = parent_int_mic_conc_metabolic,                        
            fuinc_int_parent           = fuinc_int_parent,                                      
            parent_int_mic_conc_binding     = parent_int_mic_conc_binding,      
            #incubation_plasma_parent  = incubation_plasma_parent,       
            #Clint_plasma_parent        = Clint_plasma_parent,  
            
            CLR                       = CLR,

                        # Parametes for flow 
            Qcardiac                  = Qcardiac, 
            Qgut                      = Qgut,
            Qkidney                   = Qkidney,
            Qliver                    = Qliver,
            Qadipose                  = Qadipose,
            Qbrain                    = Qbrain,
            Qrest                     = Qrest,   
            
            # volume 
            Vart                      = Vart,
            Vgut                      = Vgut,
            Vkidney                   = Vkidney,
            Vliver                    = Vliver,
            Vlung                     = Vlung,
            Vbrain                    = Vbrain,
            Vadipose                  = Vadipose,
            Vrest                     = Vrest,
            Vven                      = Vven)  

initState <- c(Agutlumen_parent    = 0, 
               Agut_parent         = 0,
               Egut                = 1,
               
               Aliver_parent       = 0,
               Eliver              = 1,
               
               Abrain_parent       = 0,
               Aadipose_parent     = 0,
               Arest_parent        = 0,
               Aven_parent         = 0,
               Alung_parent        = 0,
               Aart_parent         = 0,
               AUC_Cplasma_parent  = 0,
               AUC_Cblood_parent   = 0,
               Aurine_parent       = 0,
               
               Agut_parent_in      = 0, 
               Aven_parent_in      = 0,
               Agut_parent_out     = 0,
               Aliver_parent_out   = 0,
               Aven_parent_out     = 0,
               Aint_parent_out     = 0,
               Aart_parent_out     = 0)  



df <- ode(y = initState, 
          times = Times, 
          func = pbpk8cpt, 
          parms = parms,
          atol = 1e-8, rtol = 1e-10,
          events  = list(data = Dose_events_oral),
          method = 'lsoda')
  
df  <- as.data.frame(df)
colnames(df)[1] <- "Time"

#write.csv(df, paste('C:/Users/s1036120/OneDrive - Syngenta/HTTK/Mammals/PPZ/Plots/', species, '_oral_', oral_mg.kg, 'mgkg.csv', sep = ''), row.names = FALSE)

####### Plot (oral)
# Theme Reference: https://ggplot2.tidyverse.org/reference/ggtheme.html
# 2017 in vivo data' 1, 10, 100 mg/kg oral; 
path     <- "C:/Users/s1036120/OneDrive - Syngenta/AI/Propiconazole/PBK_Modeling/Reference compounds/Difenoconazole/Figures/DIFEN.xlsx"
obs      <- read_xlsx(path,  sheet = 'Sheet1')
obs$Time_h    <- as.numeric(obs$Time_h)
obs$Conc_ugl  <- as.numeric(obs$Conc_ugl)
oral   <- obs[obs$Matrix == 'blood'& obs$Method == 'oral' & obs$Dose_mg_kg == oral_mg.kg,]
oral   <- na.omit(oral)
oral   <- oral[oral$Species == 'CD1Mice'| oral$Species == 'C57BL/6Mice',]

iv   <- obs[obs$Matrix == 'blood'& obs$Method == 'iv',]
iv   <- na.omit(iv)
iv   <- iv[iv$Species != 'hCar/hPxrMice' & iv$Species != 'Car/PxrKOMice' & iv$Species != 'WTMice',]

# Mass balance
ggplot() +
  geom_line (data = df, aes(Time, Mass_parent_bal), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14)

ggplot() +
  geom_line (data = df, aes(Time, Mass_parent_out/oral_input*100), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14)+ xlim(0, 10)


ggplot() +
  geom_line (data = df, aes(Time,C_blood_parent), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20)) + theme_bw(base_size = 14) + xlim(0, 48)


##################################        Toxicity endpoint derivation       ################################
# AUC of blood (acute, single oral dose at 22.4 mg/kg)
AUC24   <- df[round(df$Time,1) == 48,]$AUC_Cblood_parent - df[round(df$Time,1) == 0,]$AUC_Cblood_parent
AUC24    # 15.8 umol*L/h

#s
dvalue  <- oral_mg.kg
dose    <- paste(oral_mg.kg, 'mgkg', sep = '')

coeff <- 1/250#1/80#1/4#1/120#1/4#1/40#1/40

ggplot() + 
  geom_line (data =df, aes(Time/24, C_blood_parent, color = 'Predicted blood conc'),  lwd=0.7) + 
  # add data for the second y-axis
  geom_line(data = df, aes(x = Time/24, y = (Eliver /coeff), color = "Predicted CYP3A4 conc"), lwd=0.8) +
   ylab(expression("Blood Concentration ("*mu*"mol/L)")) + 
  scale_y_continuous(sec.axis = sec_axis(~ . * coeff, name = "Relative CYP3A4 Concentration [fold]")) + 
  xlab("Time (d)") + theme(text = element_text(size = 20)) + xlim(0, 15) + theme_bw(base_size = 12)+
  # legend
  scale_color_manual(name = NULL,
                     values = c('#276DC2', "#228C22"),
                     #values = c('#00A1D5FF', "#ffc425"),
                     breaks=c("Predicted blood conc", "Predicted CYP3A4 conc")) +
  scale_shape_manual(name=NULL, values=c(16)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = c(0.8, 0.9),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm')) +
  # Titles and Labels
  ggtitle(
    label = paste("Human [oral dose of ", dvalue, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Propiconazole (PPZ)")) +
  theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=10))

ggsave(paste('C:/XXX/Mammals/PPZ/Plots/HumanwCYP_', dose, '_oral_0509.tiff', sep = ''), width = 6.3, height = 4.5, dpi = 600, compression = 'lzw')



#==========================================================================================================
###                                   Across-species extrapolation                                      ###
#==========================================================================================================
# Rat AUC24 is 14.75891 umol/L/h based on 30 mg/kg/d repeated dose
# 14 days of exposure; 3 meals taken at 7 am, 12 pm and 6pm (0, 5, 11 hr)
days          <- 20#10#14
Time_max      <- days * 24                             # h
StartTime     <- 0                                     # Time_min
StopTime      <- 40*24
dt            <- 0.1
Times         <- seq(StartTime, StopTime, dt)

# dose
Oral_mg_kg         <- 37.9#19.8 #66.8                          # mg/kg BW/d
Oral_input         <- Oral_mg_kg * 1000 / MW_parent            # mg/kg -> ug/kg bw

# # dose pattern 47/91 Sitovition report
var                <- rep("Agutlumen_parent", 3 * days)
time               <- c(0, 5, 11)
zz                 <- time

for (i in 1:(days-1)){
  temp  <- time + 24 * i
  print(temp)
  zz <- cbind(zz, temp)
}

zz             <- data.frame(zz)
time           <- unlist(zz)

value          <- c(Oral_input/3, Oral_input/3, Oral_input/3)            # P47/91 sitovition report; 3 meals taken at 7 am, 12 pm and 6pm (0, 5, 11 hr)
Dose_events    <- data.frame(var,
                             time,
                             value,
                             method = "add")


# #########################    parameters for exposure scenario    #########################
dfC <- ode(y = initState,
           times = Times,
           func  = pbpk8cpt,
           parms = parms,
           atol  = 1e-8, rtol = 1e-10,
           events  = list(data = Dose_events),
           method  = 'lsoda')

dfC  <- as.data.frame(dfC)
colnames(dfC)[1] <- "Time"

ggplot() +
  geom_line (data = dfC, aes(Time,C_blood_parent), col="#00AFBB", lwd=1) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20)) + theme_bw(base_size = 14) #+ xlim(216,240)

# AUC of blood (chronic, single oral dose at  mg/kg)
AUC24   <- dfC[round(dfC$Time,1) == 240,]$AUC_Cblood_parent - dfC[round(dfC$Time,1) == 216,]$AUC_Cblood_parent
AUC24    # 11.6761 umol*L/h (45.6 mg/kg)   # 5.010958

#s
dvalue  <- Oral_mg_kg
dose    <- paste(oral_mg.kg, 'mgkg', sep = '')

coeff <- 1#1/3#1/40

ggplot() + 
  geom_line (data =dfC, aes(Time/24, C_blood_parent, color = 'Predicted blood conc'),  lwd=0.4) + 
  # add data for the second y-axis
  geom_line(data = dfC, aes(x = Time/24, y = (Eliver /coeff), color = "Predicted CYP3A4 conc"), lwd=0.6) +
  ylab(expression("Blood Concentration ("*mu*"mol/L)")) + 
  scale_y_continuous(sec.axis = sec_axis(~ . * coeff, name = "Relative CYP3A4 Concentration [fold]")) + 
  xlab("Time (d)") + theme(text = element_text(size = 20)) + xlim(0, 15) + theme_bw(base_size = 12)+
  # legend
  scale_color_manual(name = NULL,
                     values = c('#276DC2', "#228C22"),
                     #values = c('#00A1D5FF', "#ffc425"),#7CAE00
                     breaks=c("Predicted blood conc", "Predicted CYP3A4 conc")) +
  scale_shape_manual(name=NULL, values=c(16)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(#legend.position = c(0.8, 0.9),
        legend.position = c(0.8, 0.78),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm')) +
  # Titles and Labels
  ggtitle(
    label = paste("Human [oral dose of ", dvalue, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Propiconazole (PPZ)")) +
  theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=10))


ggsave(paste('C:/XXX/Mammals/PPZ/Plots/HumanwCYP_', Oral_mg_kg, 'mgkg_oral_0405.tiff', sep = ''), width = 6.3, height = 4.5, dpi = 600, compression = 'lzw')

