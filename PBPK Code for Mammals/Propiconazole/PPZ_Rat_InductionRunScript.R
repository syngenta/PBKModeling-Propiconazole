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


#==========================================================================================
#                                         Start                                           #
#==========================================================================================
# Reference
# PKNCA : https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk  : https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

# - Target compound : PPZ (parent)
# - Species         : Rat (Maternal only)
# - Author          : Yaoxing Wu
# - Date            : Jun, 2022
# - Structure       : GI tract, Liver, Lungs, Adipose, Brain, Rest, Vein, Arterial

# Clean the environment
rm(list = ls())
#dev.off()

species                 <- 'Rat'  

source(paste("C:/XXX/General_code/Maternal_8compt_InductionModel_Parent.R"))
source(paste("C:/XXX/General_code/Species_PhyData.R"))

MPPGL      <- HPGL_rat
MPPGGI     <- MPPGGI_rat
kdeg       <- kdeg_rat

BW_mouse    <- 0.02
BW_rat      <- 0.25
BW_human    <- 70



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
fub_parent              <-  0.05 #0.0955                                 # Unbound fraction in plasma; measured in sed
Rblood2plasma_parent    <-  0.758                                        #  if not available, place "None" here.
Papp_parent             <-  40.8                                         # unit: 10-6 cm/s                 
#a <- chem.physical_and_invitro.data[chem.physical_and_invitro.data$CAS == '119446-68-3',]

BW_scaledfrom                    <- BW_mouse
# in vitro metabolic-related properties   
Vmax_unit                        <- 'umol/h/million cells'                   # 'umol/h/kg bw' or 'none'
# liver
type_hep_clearance               <- 'liver'
incubation_hep_parent            <- 'hepatocytes'                           # for rat metabolic stability test (hepatic)
Km_hep_parent                    <- 4.6                                       # unit: uM (umol/L)
Vmax_hep_parent                  <- 164 / 1E6 * 60 * Km_hep_parent           # L/h/mg * umol/L  -> umol/h/million cells (rat) 

parent_hep_conc_metabolic        <- 0.5                                   # 0.5 mg/ml protein; Shen et al. 2013
fuinc_hep_parent                 <- 9999                                  # binding was not measured in in vitro metabolsim assays
parent_hep_conc_binding          <- 9999
Clint_ori_hep_parent             <- 9999
# intestine
type_int_clearance               <- 'intestine'
incubation_int_parent            <- 'microsome'                           # for rat intestinal stability test
Vmax_int_parent                  <- Vmax_hep_parent * 2.4 / MPPGL_rat  * MPPGGI_rat          # unit: nmol/min/mg protein; Table 1 -->   umol/h/mg protein; Intestinal metabolism determined at 50 ug/mL Frost 2017
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
EC50         <- 3.57                 # uM mouse
Emax         <- 7.5                  # fold
E0           <- 1 
fmCYP        <- 0.65                 # CYP3A4

# absorption extent
fa           <- 1                                                         # absorption fraction

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
calc_ka     <- function(dose.mg.kg, Papp){
  if(is.numeric(Papp)){
    Peff_human  <- 0.4926*log10(Papp) - 0.1454                  # 1e-4 cm/s; Papp in the unit of 1e-6 cm/s
    Peff_rat    <- (Peff_human - 0.03)/3.6                      # 1e-4 cm/s
    ka_rat      <-  2 * Peff_rat / R_rat / 1e4 * 3600           # h-1
    ka          <-  ka_rat * (BW/BW_rat)^(-0.25)                # eqaution A11 in appendix
  }else{
    #ka_scaledfrom   <- exp(1.2121) * dose.mg.kg^(- 0.4476)                   # dose unit: mg/kg; chemical specific and subject to change_mouse
    #ka_scaledfrom   <- exp(0.9930) * dose.mg.kg^(-0.3667)
    ka              <- 2#exp(0.9930) * dose.mg.kg^(-0.3667)
    #ka              <- ka_scaledfrom  * (BW/BW_scaledfrom)^(-0.25)                       # eqaution A11 in appendix
  }
  
  return(ka)                                                # unit: 1/h
}
# test: ka         <- calc_ka(dose.mg.kg = 400, Papp = 'None')


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
oral_mg.kg        <- 1000#400#150#30#400#30#150#30                  # mg/kg
oral_input        <- oral_mg.kg  * 1000 / MW_parent   # mg/kg -> umol/kg bw

# iv dose 
iv_mg.kg          <- 1                               # 
iv_input          <- iv_mg.kg  * 1000  / MW_parent   # mg/kg -> umol/kg bw

# infusion dose
Infusion_mg.kg    <- 0                                            
Infusion_input    <- Infusion_mg.kg  * 1000 / MW_parent           # mg/kg -> ug/kg bw
CTinf             <- 1                                            # 1 h infusion time
CRATE             <- Infusion_input/CTinf                         # ug/kg bw /h
CTIMEinf          <- c(0,CTinf,Time_max)
CRATEinf          <- c(CRATE,0,0)
# Define an interpolation function that returns rate when given time - "const"
Cstep.doseinf     <- approxfun(CTIMEinf, CRATEinf, method = "const")

# absorption calculation
ka                <- calc_ka(dose.mg.kg = oral_mg.kg, Papp = Papp_parent)    # h-1; 'None'
kt                <-  kt_rat                                                 # 0.11 # 0.1 #0.15 # 0.08  # 0.15  h-1

ndays        <-  2#13 #9#0#2#13 #0#2#13#13#9
Dose_events_oral   <- data.frame(var = "Agutlumen_parent",                   # Agutlumen_parent; Aven_parent
                                 time = seq(from = 0, to = ndays * 24, by = 24),
                                 value = rep(oral_input,ndays+1),
                                 method = "add")

Dose_events_SingleOral   <- data.frame(var    = "Agutlumen_parent",          # Agutlumen_parent; Aven_parent
                                       time   = 0,
                                       value  = oral_input,
                                       method = "add")

Dose_events_iv      <- data.frame(var = "Aven_parent",           # Agutlumen_parent; Aven_parent
                                  time = 0,
                                  value = iv_input,
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
            
            EC50                      = EC50,                  # uM
            Emax                      = Emax,                  # fold
            E0                        = E0, 
            kdeg                      = kdeg,

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
            Km_int_parent              = Km_int_parent,                                            
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
          events  = list(data = Dose_events_SingleOral), #events  = list(data = Dose_events_oral),
          method = 'lsoda')
  
df  <- as.data.frame(df)
colnames(df)[1] <- "Time"

write.csv(df, paste('C:/XXX/PPZ/Plots/', species, '_oral_', oral_mg.kg, 'mgkg.csv', sep = ''), row.names = FALSE)

####### Plot (oral)
# Theme Reference: https://ggplot2.tidyverse.org/reference/ggtheme.html
path     <- "C:/XXX/Difenoconazole/Figures/DIFEN.xlsx"
obs      <- read_xlsx(path,  sheet = 'Sheet1')
obs$Time_h    <- as.numeric(obs$Time_h)
obs$Conc_ugl  <- as.numeric(obs$Conc_ugl)
oral   <- obs[obs$Matrix == 'blood'& obs$Method == 'oral' & obs$Dose_mg_kg == oral_mg.kg,]
oral   <- na.omit(oral)
oral   <- oral[oral$Species == 'CD1Mice'| oral$Species == 'C57BL/6Mice',]

iv   <- obs[obs$Matrix == 'blood'& obs$Method == 'iv',]
iv   <- na.omit(iv)
iv   <- iv[iv$Species != 'hCar/hPxrMice' & iv$Species != 'Car/PxrKOMice' & iv$Species != 'WTMice',]

#obs_plasma$Conc_uM    <- obs_plasma$Conc_ugl / MW_parent
#obs_brain$Conc_uM     <- obs_brain$Conc_ugl  / MW_parent
#obs_liver$Conc_uM     <- obs_liver$Conc_ugl  / MW_parent
#obs_lung$Conc_uM      <- obs_lung$Conc_ugl   / MW_parent

#obs_df           <- obs[obs$Matrix!='lung' & obs$Matrix!= 'plasma',]
#obs_df$Matrix    <- sub("^(.)", "\\U\\1", obs_df$Matrix , perl = TRUE)   # uppercase the first letter in the column
#obs_df$Matrix    <- paste(obs_df$Matrix, ' (exp.)', sep='')
#obs_df$Conc_uM   <- obs_df$Conc_ugl  / MW_parent

# Mass balance
ggplot() +
  geom_line (data = df, aes(Time, Mass_parent_bal), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14)

ggplot() +
  geom_line (data = df, aes(Time,C_blood_parent), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20)) + theme_bw(base_size = 14) #+ xlim(0,24)



#==========================================================================================================
###                                   Across-species extrapolation                                      ###  
#==========================================================================================================
##################################        Toxicity endpoint derivation       ################################
# Acute
# AUC of blood (acute, single oral dose at 30 mg/kg)
AUC24   <- df[round(df$Time,1) == 48,]$AUC_Cblood_parent - df[round(df$Time,1) == 0,]$AUC_Cblood_parent
AUC24    # 15.80492 umol*L/h


# oral gavage
# Chronic
# AUC of blood (chronic, single oral dose at 30 mg/kg)
AUC24   <- df[round(df$Time,1) == 240,]$AUC_Cblood_parent - df[round(df$Time,1) == 216,]$AUC_Cblood_parent
AUC24    # 11.69799 umol*L/h 

# 
ggplot() +
  geom_line (data = df, aes(Time/24, Eliver), col="#00AFBB", lwd=2) +
  #geom_line (data = df, aes(Time/24, Egut ), col="red", lwd=2) + ylab("Net fold increase") +
  xlab("Time (d)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14) + geom_vline(xintercept=4)#+ ylim(0,1)



#######################   Rat Experimental data   #########################
CYP3A1_dose        <- c(150,     400)
CYP3A1_time        <- c(14,      4)
CYP3A1_response    <- c(2.213, 2.106)
CYP3A1_responseSD  <- c(0.8,   0.2766)
class               <- c('150 mg/kg BW', '400 mg/kg BW')
CYP3A1             <- data.frame(CYP3A1_dose, CYP3A1_time, CYP3A1_response, CYP3A1_responseSD, class)
CYP3A1$class       <- factor(CYP3A1$class, levels =  c('150 mg/kg BW', '400 mg/kg BW'))

dvalue  <- oral_mg.kg
dose    <- paste(oral_mg.kg, 'mgkg', sep = '')

coeff <- 1/100#1/100#1/40#1/80#1/20 #1/100#1/30

ggplot() + 
  geom_point(data = CYP3A1[CYP3A1$CYP3A1_dose == dvalue,],aes(CYP3A1_time, CYP3A1_response/coeff, shape = "Cyp3a1 observation"),  color = "#8c8c8c",  size=3.5, alpha = 1)+   #  colour = "#FC4E07",
  geom_errorbar(data = CYP3A1[CYP3A1$CYP3A1_dose == dvalue,], aes(CYP3A1_time, 
                                                                  ymin= CYP3A1_response/coeff - CYP3A1_responseSD/coeff, 
                                                                  ymax = CYP3A1_response/coeff + CYP3A1_responseSD/coeff),  color = "#8c8c8c", width=0.5, linewidth = 0.8) +
  
  geom_line (data =df, aes(Time/24, C_blood_parent, color = 'Predicted blood conc'),  lwd=0.5) + 
  # add data for the second y-axis
  geom_line(data = df, aes(x = Time/24, y = (Eliver /coeff), color = "Predicted Cyp3a1 conc"), lwd=0.7) +
  ylab(expression("Blood Concentration ("*mu*"mol/L)")) + 
  scale_y_continuous(sec.axis = sec_axis(~ . * coeff, name = "Relative Cyp3a1 Concentration [fold]")) + 
  xlab("Time (d)") + theme(text = element_text(size = 20)) + xlim(0, 15) + theme_bw(base_size = 12)+
  # legend
  scale_color_manual(name = NULL,
                     values = c("#00A1D5FF", '#ffc425'), ##7CAE00
                     breaks=c("Predicted blood conc", "Predicted Cyp3a1 conc")) +
  scale_shape_manual(name=NULL, values=c(16)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = c(0.8, 0.84),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm')) +
  # Titles and Labels
  ggtitle(
    label = paste("Rat [oral dose of ", dvalue, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Propiconazole (PPZ)")) +
  theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=10))

ggsave(paste('xxx/PPZ/Plots/RatwCYP_', dose, '_oral_1007.tiff', sep = ''), width = 6.3, height = 4.5, dpi = 600, compression = 'lzw')


