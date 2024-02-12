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
# - Species         : Mouse dietary dose (Maternal only)
# - Author          : Yaoxing Wu
# - Date            : Jun, 2022
# - Structure       : GI tract, Liver, Lungs, Adipose, Brain, Rest, Vein, Arterial

# Clean the environment
rm(list = ls())
#dev.off()

species                 <- 'Mouse'  

source(paste("C:/XXX/General_code/Maternal_8compt_InductionModel_AdLibitum_Parent.R"))
source(paste("C:/XXX/General_code/Species_PhyData.R"))

MPPGL      <- MPPGL_mouse 
MPPGGI     <- MPPGGI_mouse
kdeg       <- kdeg_mouse

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
fub_parent              <-  0.05                                         # Unbound fraction in plasma; measured in sed
Rblood2plasma_parent    <-  0.579                                        # if not available, place "None" here.
Papp_parent             <-  40.8                                         # unit: 10-6 cm/s                 


# in vitro metabolic-related properties   
Vmax_unit                        <- 'umol/h/kg bw'                   # 'umol/h/kg bw' or 'none'

# liver use difenoconazole value
type_hep_clearance               <- 'liver'
incubation_hep_parent            <- 'scaled'                           # for rat metabolic stability test (hepatic)
Km_hep_parent                    <- 4.4                                # unit: uM (umol/L)
BW_scaledfrom                    <- BW_rat 
df_tissue_rat                    <- tissue.data[which(tissue.data$Species == 'Rat'),]
Vliver_rat                       <- subset(df_tissue_rat , variable == "Vol (L/kg)" & 
                                             tolower(Tissue) == 'liver')$value            # unit: uM (umol/L)
Vmax_hep_parent_rat              <- 145.8E-6 * 60 * Km_hep_parent * 17 / 30.6             # L/h/mg * uM ->   umol/h/mg 322
Vmax_hep_parent                  <- Vmax_hep_parent_rat * MPPGL_rat * (Vliver_rat * 1.05 * 1000)         # umol/h/mg protein * mg/g tissue * mL/kg bw -> umol/h/kg bw

parent_hep_conc_metabolic        <- 0.5                                   
fuinc_hep_parent                 <- 9999                                 
parent_hep_conc_binding          <- 9999
Clint_ori_hep_parent             <- 9999
# intestine
type_int_clearance               <- 'intestine'
incubation_int_parent            <- 'scaled'                         
Vmax_int_parent                  <- Vmax_hep_parent / MPPGL * MPPGGI        
Km_int_parent                    <- Km_hep_parent                         # unit: uM (umol/L)
parent_int_mic_conc_metabolic    <- 9999                                  # Intestinal metabolism
fuinc_int_parent                 <- 9999                              
parent_int_mic_conc_binding      <- 9999  
# plasma
type_plasma_clearance            <- 'plasma'
incubation_plasma_parent         <- 'half-life'   
Clint_plasma_parent              <-  1e6                                  # half life (h)
fu_plasma_parent                 <- 'None'

# Induction related parapeters
EC50         <- 7.09                                                      # uM mouse
Emax         <- 12.7                                                      # fold
E0           <- 1 
fmCYP        <- 0.65                                                      # CYP3A4

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
source(paste("C:/XXX/General_code/Maternal_DM_Parent.R"))
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
    #ka  <- exp(1.2121) * dose.mg.kg^(- 0.4476)                   # dose unit: mg/kg; chemical specific and subject to change
    ka      <- exp(0.9930) * dose.mg.kg^(-0.3667) 
  }
  
  return(ka)                                                # unit: 1/h
}
# test: ka         <- calc_ka(dose.mg.kg = 10)


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
days          <- 30
Time_max      <- days * 24                             # h
StartTime     <- 0                                     # Time_min 
StopTime      <- Time_max
dt            <- 0.1
Times         <- seq(StartTime, StopTime, dt)

# oral dose # ad libitum dose
oral_mg.kg        <- 149#578                              # mg/kg
oral_input        <- oral_mg.kg  * 1000 / MW_parent    # mg/kg -> umol/kg bw

ndays             <- 13                                # number of feeding days 0-13
var               <- rep("Agutlumen_parent", 24 * (ndays+1)/0.5)
time              <- seq(0, (24 * (ndays+1)), by = 0.5)
time              <- time[-length(time)]
value_1           <- read_xlsx('C:/XXX/General_code/Data/Mouse_FoodIntake.xlsx', sheet = 'Combined')$'Intake_%'
value_percent     <- c(rep(value_1, (ndays+1)))
value             <- value_percent / 100 * oral_input 

Dose_events_oral_AdLib    <- data.frame(var,
                                        time,
                                        value,
                                        method = "add")

Dose4ka             <- data.frame(time,
                                  value)        # check: sum(Dose4ka$value[1:48]/ 1000 *  MW_parent, na.rm = TRUE)
Times_df            <- data.frame(Times)
colnames(Times_df)  <- c("time")
Dose4ka_df          <- merge(Times_df, Dose4ka, by = "time", all = TRUE)
Dose4ka_df          <- Dose4ka_df %>% fill(value)        # Replacing missing data from top to bottom
Dose4ka_df$value    <- Dose4ka_df$value / 1000 *  MW_parent    # mg/kg


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
kt                <-  0.11                                        #   h-1

#########################    parameters for exposure scenario    #########################
parms <- c( Oral_mg.kg                = oral_mg.kg, 
            Oral_input                = oral_input,   
            
            pH                        = pH,
            fa                        = fa,             
            
            Rblood2plasma_parent      = Rblood2plasma_parent,
            BW                        = BW,              
            
            fub_parent                = fub_parent,   
            kt                        = kt,
            
            EC50                      = EC50,                  # uM
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
          events  = list(data = Dose_events_oral_AdLib),
          method = 'lsoda')
  
df  <- as.data.frame(df)
colnames(df)[1] <- "Time"

write.csv(df, paste('C:/XXX/PPZ/Plots/', species, '_oralAdLib_', oral_mg.kg, 'mgkg.csv', sep = ''), row.names = FALSE)

#df <- read.csv(paste('C:/XXX/PPZ/Plots/', species, '_oralAdLib_', oral_mg.kg, 'mgkg.csv', sep = ''))
####### Plot (oral)
# Theme Reference: https://ggplot2.tidyverse.org/reference/ggtheme.html
path     <- "C:/XXX/Reference compounds/Difenoconazole/Figures/DIFEN.xlsx"
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
  geom_line (data = df, aes(Time/24, Eliver), col="#00AFBB", lwd=2) + 
  #geom_line (data = df, aes(Time/24, Egut ), col="red", lwd=2) + ylab("Net fold increase") +
  xlab("Time (d)") + theme(text = element_text(size = 20))  + theme_bw(base_size = 14) + geom_vline(xintercept=14)#+ ylim(0,1)


# plot
#######################   Mouse Experimental data   #########################
#oral_mg.kg <- 149
#df                  <- read.csv(paste('C:XXX/', species, '_oralAdLib_', oral_mg.kg, 'mgkg.csv', sep = ''))

Cyp3a11_dose        <- c(150,  210,  149, 578)
Cyp3a11_time        <- c(14,     4,   14,   14)
Cyp3a11_response    <- c(3.7, 6.55, 3.66, 5.24)
Cyp3a11_responseSD  <- c(2.0, 3.32,  0.3, 0.52)
class               <- c('150 mg/kg BW', '210 mg/kg BW', '850 ppm Ad-lib (149 mg/kg BW)', '2500 ppm Ad-lib (578 mg/kg BW)')
Cyp3a11             <- data.frame(Cyp3a11_dose, Cyp3a11_time, Cyp3a11_response, Cyp3a11_responseSD, class)
Cyp3a11$class       <- factor(Cyp3a11$class, levels =  c('150 mg/kg BW/d', '210 mg/kg BW/d', '850 ppm Ad-lib (149 mg/kg BW/d)', '2500 ppm Ad-lib (578 mg/kg BW/d)'))

dvalue  <- oral_mg.kg
dose    <- paste(oral_mg.kg, 'mgkg', sep = '')

coeff <- 2 #2 #1/2

ggplot() + 
  geom_point(data = Cyp3a11[Cyp3a11$Cyp3a11_dose == dvalue,],aes(Cyp3a11_time, Cyp3a11_response/coeff, shape = "Cyp3a11 observation"),  color = "#8c8c8c",  size=3.5, alpha = 0.8)+   #  colour = "#FC4E07",
  geom_errorbar(data = Cyp3a11[Cyp3a11$Cyp3a11_dose == dvalue,], aes(Cyp3a11_time, 
                    ymin= Cyp3a11_response/coeff - Cyp3a11_responseSD/coeff, 
                    ymax = Cyp3a11_response/coeff + Cyp3a11_responseSD/coeff),  color = "#8c8c8c", width=0.7, linewidth = 1) +
  geom_line (data =df, aes(Time/24, C_blood_parent, color = 'Predicted blood conc'),  lwd=0.4) + 
  # add data for the second y-axis
  geom_line(data = df, aes(x = Time/24, y = (Eliver /coeff), color = "Predicted Cyp3a11 conc"), lwd=0.6) +
  ylab(expression("Blood Concentration ("*mu*"mol/L)")) + 
  scale_y_continuous(sec.axis = sec_axis(~ . * coeff, name = "Relative Cyp3a11 Concentration [fold]")) + 
  xlab("Time (d)") + theme(text = element_text(size = 20)) + xlim(0, 25) + theme_bw(base_size = 12)+
  # legend
  scale_color_manual(name = NULL,
                     values = c( "#00A1D5FF", '#ffc425'),
                     breaks=c("Predicted blood conc", "Predicted Cyp3a11 conc")) +
  scale_shape_manual(name=NULL, values=c(16)) +
  scale_linetype_manual(name = NULL,values='solid')+ guides(color = guide_legend(override.aes = list(shape = 16)))+
  theme(legend.position = c(0.83, 0.88),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid",
                                         colour ="white"),
        legend.spacing.y = unit(-0.2, 'cm')) +
  # Titles and Labels
  ggtitle(
    label = paste("Mouse [Ad libitum oral dose of ", dvalue, " mg/kg BW/d]", sep = ''),
    subtitle = paste("Propiconazole (PPZ)")) +
  theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=10))

ggsave(paste('C:/XXX/PPZ/Plots/MousewCYP_AdLib_', dose, '_oral_1007.tiff', sep = ''), width = 7.7, height = 5.5, dpi = 600, compression = 'lzw')



# 
