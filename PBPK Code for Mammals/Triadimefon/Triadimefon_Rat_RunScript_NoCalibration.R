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


#==========================================================================================
#                                         Start                                           #
#==========================================================================================
# Reference
# PKNCA : https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk  : https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

# - Target compound : Triadimefon (parent)
# - Species         : Rat (Maternal only)
# - Author          : Yaoxing Wu
# - Date            : Jun, 2022
# - Structure       : GI tract, Liver, Lungs, Adipose, Brain, Rest, Vein, Arterial

# Clean the environment
rm(list = ls())
#dev.off()

species                 <- 'Rat'  

source(paste("C:/XXX/Mammals/General_code/Maternal_8compt_Model_Parent.R"))
source(paste("C:/XXX/Mammals/General_code/Species_PhyData.R"))


MPPGL      <- MPPGL_rat
MPPGGI     <- MPPGGI_rat

Qrest       ==  Qcardiac - Qliver - Qgut - Qadipose - Qbrain 
BW_rat      <- 0.25
Vliver_rat  <- 0.03486              
#======================================================================================================
#                          Compound physico-chemical & in vitro properties                            #
#======================================================================================================
#######################    Parent   ##########################
# physico-chemical properties   
Compound_name_parent    <- 'Triadimefon(TDF)'
Compound_type_parent    <- 'base'
CAS_parent              <- "43121-43-3-a"                                # add a to CAS since acibenzolar has been included in HTTK database
LogP_parent             <-  2.77                                          # log Kow; 3.59 from RAAV at 25c; 3.68 at 20c at pubchem.
MW_parent               <-  293.8 
pKa_a_parent            <-  7.22                                         # very weak base without any numerical value; chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://hal-unilim.archives-ouvertes.fr/hal-00929974/document
pKa_b_parent            <-  1.93                                          # pka accept (base); Teb?H+ -> Teb + H+ 4.89 from httk

# in vitro properties   -
fub_parent              <-  0.1                                          # Unbound fraction in plasma; measured in sed
Rblood2plasma_parent    <-  'None'                                        # Table 1; assumed; if not available, place "None" here.
Papp_parent             <-  38                                            # unit: 10-6 cm/s                 
#a <- chem.physical_and_invitro.data[chem.physical_and_invitro.data$CAS == '43121-43-3',]
# in vitro metabolic-related properties   
Vmax_unit                        <- 'none'                        # 'umol/h/kg bw' or 'none'
# liver
type_hep_clearance               <- 'liver'
incubation_hep_parent            <- 'microsome'                              # for rat metabolic stability test (hepatic)
Vmax_hep_parent                  <- 0.478                             #  unit: umol/h/mg microsome
Km_hep_parent                    <- 47.3                                   # unit: uM (umol/L)
parent_hep_conc_metabolic        <- 9999                                  # hepatocyte metabolic experiment determined at 1000 cells/mL; page 13/90 of frost 2017
fuinc_hep_parent                 <- 9999                                  # measured at 1000 cells/mL hepatocyte binding experiment; Frost 2017 (17/90)
parent_hep_conc_binding          <- 9999
Clint_ori_hep_parent             <- 'None'

# intestine
type_int_clearance               <- 'intestine'
incubation_int_parent            <- 'microsome'                                # for rat intestinal stability test
Vmax_int_parent                  <- Vmax_hep_parent / MPPGL_rat * MPPGGI_rat                  # unit: nmol/min/mg protein; Table 1 -->   umol/h/mg protein; Intestinal metabolism determined at 50 ug/mL Frost 2017
Km_int_parent                    <- 47.3                                    # unit: uM (umol/L)
parent_int_mic_conc_metabolic    <- 9999                                    # Intestinal metabolism determined at 50 ug/mL Frost 2017
fuinc_int_parent                 <- 9999                                    # Measured at 1 ug/mL microsome binding experiment; Frost 2017 (16/90)
parent_int_mic_conc_binding      <- 9999   

incubation_plasma_parent         <- 'plasma'   
incubation_plasma_parent         <- 'half-life'   
Clint_plasma_parent              <-  1e6                                  # half life (h)
fu_plasma_parent                 <- 'None'

a <- chem.physical_and_invitro.data[chem.physical_and_invitro.data$CAS == '43121-43-3',]
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
## predicting the tissue to unbound plasma partition coefficients for the tissues contained in the tissue.data table 
source(paste("C:/Users/s1036120/OneDrive - Syngenta/HTTK/Mammals/General_code/Maternal_DM_Parent.R"))

# Validation only
tissuelist <- list(liver=c("liver"),lung=c("lung"),gut=c("gut"), brain=c("brain"), adipose=c("adipose"),
                    muscle.bone=c('bone', 'heart', 'skin', 'spleen', 'muscle', 'rest', "kidney"))
round(lump_tissues(PC_parent,tissuelist=tissuelist,species = species)$Kmuscle.bone2pu,2)  == round(Krest2pu_parent, 2) 
 


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
    Peff_human  <- 0.4926*log10(Papp) - 0.1454                  # 1e-4 cm/s; Papp in the unit of 1e-6 cm/s
    Peff_rat    <- (Peff_human - 0.03)/3.6                      # 1e-4 cm/s
    ka_rat      <-  2 * Peff_rat / R_rat / 1e4 * 3600           # h-1
    ka          <-  ka_rat * (BW/BW_rat)^(-0.25) 
    #kt          <-  0.68                                        # Transit rate constant (1/h) 0.68'
    
  }else{
    #ka_scaledfrom   <- exp(1.2121) * dose.mg.kg^(- 0.4476)                   # dose unit: mg/kg; chemical specific and subject to change_mouse
    ka_scaledfrom   <- 2
    ka              <- ka_scaledfrom  * (BW/BW_scaledfrom)^(-0.25)                       # eqaution A11 in appendix
  }
  
  return(ka)                                                # unit: 1/h
}
#ka         <- calc_ka(dose.mg.kg = 10)

# fa was obtained from (U.S.EPA Dec 12, 2017) following oral gavage administration of radiolabeled acibenzolar at doses of 0.5 mg/kg and 100 mg/kg.
# Reference: Acibenzolar manuscript
fa        <- 1                                           # absorption fraction



##########################################################################
###                        excretion L/h/kg                            ###  
##########################################################################
#  renal clearance for TEB is negligible 
CLR              <- 0         # L/h/kg no urinary excretion for teb parent compound.

calc_GFR         <- function(CLR){
  if(CLR == "None"){
    Qgfr       <- subset(df_physio, physiology.data.Parameter == 'GFR')$"physiology.data...grepl.species..names.physiology.data..."/(BW^0.25)/1000*60     # L/h/kg
  }else{
    Qgfr       <- CLR 
  }
  
  return(Qgfr)
}
  
#----------------------------                   End of compound specific parameters         -----------------------------------  



#================================================================================================================
###                                      Define modeling variables                                            ###  
#================================================================================================================
# in vivo experiment
days          <- 3
Time_max      <- days * 24                             # h
StartTime     <- 0                                     # Time_min 
StopTime      <- Time_max
dt            <- 0.01

Times         <- seq(StartTime, StopTime, dt)

#           Dose regimen: Oral         #
# oral dose 
oral_mg.kg        <- 30                  # mg/kg
oral_input        <- oral_mg.kg  * 1000 / MW_parent   # mg/kg -> umol/kg bw

# dose 
iv_mg.kg          <- 50                             # 
iv_input          <- iv_mg.kg  * 1000  / MW_parent   # mg/kg -> umol/kg bw

Infusion_mg.kg    <- 0                                            # 1, 10 or 100; 2017 in vivo rat study
Infusion_input    <- Infusion_mg.kg  * 1000 / MW_parent           # mg/kg -> ug/kg bw
CTinf             <- 1                                            # 1 h infusion time
CRATE             <- Infusion_input/CTinf                         # ug/kg bw /h

CTIMEinf          <- c(0,CTinf,Time_max)
CRATEinf          <- c(CRATE,0,0)

# Define an interpolation function that returns rate when given time -
"const"
Cstep.doseinf     <- approxfun(CTIMEinf, CRATEinf, method = "const")

ka                <- 2              #calc_ka(dose.mg.kg = Oral_mg.kg, Papp = Papp_parent)
kt                <-  kt_rabbit   

#########################    parameters for exposure scenario    #########################
parms <- c( Oral_mg.kg                = 0, 
            Oral_input                = 0,   
            
            pH                        = pH,
            fa                        = fa,             
            
            Rblood2plasma_parent      = Rblood2plasma_parent,
            BW                        = BW,              
            
            fub_parent                = fub_parent,   
            ka                        = ka,
            kt                        = kt,

            # partiton coefficient
            Kgut2pu_parent            = Kgut2pu_parent,      
            Kliver2pu_parent          = Kliver2pu_parent,     
            Klung2pu_parent           = Klung2pu_parent,  
            Kbrain2pu_parent          = Kbrain2pu_parent,
            Kadipose2pu_parent        = Kadipose2pu_parent,
            Krest2pu_parent           = Krest2pu_parent ,     
 
            
            #incubation_hep_parent     = incubation_hep_parent,                                     # for rat metabolic stability test (hepatic)
            Vmax_hep_parent            = Vmax_hep_parent,                                       # unit: nmol/min/10^6 cells; Table 1 -->   umol/h/10^6 cellS; determined at 1000 cells/mL; page 13/90 of frost 2017
            Km_hep_parent              = Km_hep_parent,                                         # unit: uM (umol/L)
            parent_hep_conc_metabolic  = parent_hep_conc_metabolic,                  # hepatocyte metabolic experiment determined at 1000 cells/mL; page 13/90 of frost 2017
            fuinc_hep_parent           = fuinc_hep_parent,                                      # measured at 1000 cells/mL binding experiment
            parent_hep_conc_binding    = parent_hep_conc_binding,   
            #incubation_int_parent     = incubation_int_parent,                                     # for rat intestinal stability test
            Vmax_int_parent           = Vmax_int_parent,                                           # unit: nmol/min/mg protein; Table 1 -->   umol/h/mg protein; Intestinal metabolism determined at 50 ug/mL Frost 2017
            Km_int_parent             = Km_int_parent,                                             # unit: uM (umol/L)
            parent_int_mic_conc_metabolic     = parent_int_mic_conc_metabolic,                     # Intestinal metabolism determined at 50 ug/mL Frost 2017
            fuinc_int_parent          = fuinc_int_parent,                                          # Measured at 1 ug/mL binding experiment; Frost 2017
            parent_int_mic_conc_binding       = parent_int_mic_conc_binding,      
            #incubation_plasma_parent  = incubation_plasma_parent,       
            Clint_plasma_parent       = Clint_plasma_parent,  
            
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
               Aliver_parent       = 0,
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

ndays  <- 0
Dose_events_oral   <- data.frame(var = "Agutlumen_parent",           # Agutlumen_parent; Aven_parent
                                 time = seq(from = 0, to = ndays * 24, by = 24),
                                 value = rep(oral_input,ndays+1),
                                 method = "add")

Dose_events_iv      <- data.frame(var = "Aven_parent",           # Agutlumen_parent; Aven_parent
                                  time = 0,
                                  value = iv_input,
                                  method = "add")
df <- ode(y = initState, 
          times = Times, 
          func = pbpk8cpt, 
          parms = parms,
          atol = 1e-8, rtol = 1e-10,
          events  = list(data = Dose_events_iv),
          method = 'lsoda')
  
df  <- as.data.frame(df)
colnames(df)[1] <- "Time"

write.csv(df, paste('C:/XXX/AI/Propiconazole/PBK_Modeling/Reference compounds/Triadimefon/Results/', species, '_single_iv_', iv_mg.kg, 'mgkg.csv', sep = ''), row.names = FALSE)


####### Plot (oral)
# Theme Reference: https://ggplot2.tidyverse.org/reference/ggtheme.html
# 2017 in vivo data' 1, 10, 100 mg/kg oral; 
path     <- "C:/XXX/AI/Propiconazole/PBK_Modeling/Reference compounds/Triadimefon/Figures/plot.csv"
obs      <- read.csv(path)
model    <- read.csv("C:/XXX/Propiconazole/PBK_Modeling/Reference compounds/Triadimefon/Figures/iv_50mgkg.csv")

# Mass balance
ggplot() +
  geom_line (data = df, aes(Time, Mass_parent_bal+iv_input), col="#00AFBB", lwd=2) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20)) + ylim(-0.5, 0.5) + theme_bw(base_size = 14)

ggplot() +
  geom_line (data = df[-(1:5),], aes(Time, C_blood_parent), col="#00AFBB", lwd=2) + 
  geom_point(data = obs, aes(Time_h, C_blood)) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20)) + theme_bw(base_size = 14) + xlim(0, 24)

ggplot() +
  geom_line (data = df[-(1:2),], aes(Time, log10(C_blood_parent)), col="#00AFBB", lwd=2) + 
  geom_point(data = obs, aes(Time_h, log10(C_blood))) + 
  geom_point(data = model, aes(Time_h, log10(C_blood))) + ylab("Mass balance (umol)") +
  xlab("Time (h)") + theme(text = element_text(size = 20)) + theme_bw(base_size = 14) + xlim(0, 24)

# Concentration
# https://stackoverflow.com/questions/54976769/how-to-use-ggplot2-legend-to-denote-different-geoms
ggplot() +
  geom_line (data = df[-c(1:4),], aes(Time, C_blood_parent,linetype = 'Predicted (This study)'), col="#1C78C0", alpha = 0.8, lwd=1) + 
  ylab(expression("Blood Concentration ("*mu*"mol/L)")) +
  geom_line (data = model, aes(Time_h, C_blood,linetype = 'Predicted (Crowell et al, 2011)'), col="#8c8c8c", lwd=1) + # Jónsdótt
  geom_point(data = obs,aes(Time_h, C_blood, color = "Exp.", shape = 'Exp.'),size=2)+   #  colour = "#FC4E07",
  xlab("Time (h)") + theme(text = element_text(size = 20)) + xlim(0, 24) + theme_bw(base_size = 12)+
  # legend
  scale_color_manual(name = NULL, values = '#d8923d') +
  scale_shape_manual(name=NULL, values=c("Exp."=8))+
  scale_linetype_manual(name = NULL,values=c("Predicted (This study)"="solid","Predicted (Crowell et al, 2011)"="dotdash"))+
  #scale_color_manual(name=NULL, values=c("Predicted (This study)"="#00A1D5FF","Predicted (Jonsdottir et al, 2016)"="#8c8c8c")) +
  theme(legend.position = c(0.7, 0.82),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm')) + 
  # Titles and Labels
  ggtitle(
    label = paste("Rat [i.v. dose of 50 mg/kg bw]"),
    subtitle = paste("Triadimefon")) +
  theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=10))
ggsave(('C:/XXX/Propiconazole/PBK_Modeling/Reference compounds/Triadimefon/Results/Rat_50mgkg_iv_0413.tiff'), width = 5, height = 4, dpi = 600, compression = 'lzw')
