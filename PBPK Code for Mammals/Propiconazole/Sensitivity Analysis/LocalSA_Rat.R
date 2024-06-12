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
library(httk)
library(reshape2)


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

Qrest       ==  Qcardiac - Qliver - Qgut - Qadipose - Qbrain 
#======================================================================================================
#                          Compound physico-chemical & in vitro properties                            #
#======================================================================================================
#######################   Parent  #########################
# physico-chemical properties  
BW_mouse    <- 0.02
BW_rat      <- 0.25
BW_human    <- 70
pH          <- 7.45

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

BW_scaledfrom                    <- BW_mouse
# in vitro metabolic-related properties   
Vmax_unit                        <- 'umol/h/million cells'                   # 'umol/h/kg bw' or 'none'
# liver
type_hep_clearance               <- 'liver'
incubation_hep_parent            <- 'hepatocytes'                           # for rat metabolic stability test (hepatic)
Km_hep_parent                    <- 4.6                                       # unit: uM (umol/L)
Vmax_hep_parent                  <- 164 / 1E6 * 60 * Km_hep_parent           # L/h/millon cells * umol/L  -> umol/h/million cells (rat) 


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
EC50         <- 3.57
Emax         <- 7.5
E0           <- 1 
fmCYP        <- 0.65                   # CYP3A4



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
source(paste("C:/XXX/Mammals/General_code/Maternal_DM_Parent.R"))
# Validation only
tissuelist <- list(liver=c("liver"),lung=c("lung"),gut=c("gut"), brain=c("brain"), adipose=c("adipose"),
                   muscle.bone=c('bone', 'heart', 'skin', 'spleen', 'muscle', 'rest', "kidney"))
round(lump_tissues(PC_parent,tissuelist=tissuelist,species = species)$Kmuscle.bone2pu,1)  == round(Krest2pu_parent, 1) 

# absorption extent
fa                <-  1                                                         # absorption fraction
ka                <-  1   
kt                <-  kt_rat   




####################################################################
###                    clearance (L/d/kg BW)                     ###  
####################################################################
# hepatic and intestinal clearance for parent; hepatic clearance for daughter (rat only)
liver.density           <-  1.05

##################################################################
###                    excretion    L/d/kg BW                  ###  
##################################################################
# Table 3, Scanes el al 2022 Quantitative morphometric 
# Abstract, Gasthuys, Elke 2019 Comparative physiology of glomerular filtration rate by plasma clearance of exogenous creatinine and exo-iohexol in six different avian species
factor            <-  1
Qgfr_parent       <-  2.57 /1000 * 60 * 24                       ## 2.6 / 1000 * 60  * 24 * factor   L/d/kg BW L for laying hen and duck
Qgfr_daughter     <-  2.57 /1000 * 60 * 24  

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

# infusion dose
Infusion_mg.kg    <- 0                                            # 1, 10 or 100; 2017 in vivo rat study
Infusion_input    <- Infusion_mg.kg  * 1000 / MW_parent           # mg/kg -> ug/kg bw
CTinf             <- 1                                            # 1 h infusion time
CRATE             <- Infusion_input/CTinf                         # ug/kg bw /h
CTIMEinf          <- c(0,CTinf,1*24)                              # # needs to be consistant with the overall time frame!!!
CRATEinf          <- c(CRATE,0,0)
# Define an interpolation function that returns rate when given time - "const"
Cstep.doseinf     <- approxfun(CTIMEinf, CRATEinf, method = "const")

#=====================     End of user input    =======================#



#=========================================================================================
#                                 PBTK model equations                                   #
#=========================================================================================

########################################################################################################################
#                                                  PBPK MODEL RUN                                                      #
########################################################################################################################
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




#########################    parameters for exposure scenario ###############
factor               <- 1
fraction_daughter    <- 1

parms <- c( Rblood2plasma_parent      = Rblood2plasma_parent,
            BW                        = BW,              
            
            fub_parent                = fub_parent,   
            ka                        = ka,
            kt                        = kt,
            
            EC50                      = EC50,              # uM
            Emax                      = Emax,                  # fold
            kdeg                      = kdeg,
            
            # partiton coefficient
            Kgut2pu_parent            = Kgut2pu_parent,      
            Kliver2pu_parent          = Kliver2pu_parent,     
            Klung2pu_parent           = Klung2pu_parent,  
            Kbrain2pu_parent          = Kbrain2pu_parent,
            Kadipose2pu_parent        = Kadipose2pu_parent,
            Krest2pu_parent           = Krest2pu_parent ,     
            
            #incubation_hep_parent     = incubation_hep_parent,                                     # for rat metabolic stability test (hepatic)
            Vmax_hep_parent            = Vmax_hep_parent,                                       
            Km_hep_parent              = Km_hep_parent,                                             # unit: uM (umol/L)
            Vmax_int_parent            = Vmax_int_parent ,                                  
            Km_int_parent              = Km_int_parent,                                             # unit: uM (umol/L)
            CLR                        = CLR,
             
            # Parametes for flow 
            Qcardiac                  = Qcardiac, 
            Qgut                      = Qgut,
            Qliver                    = Qliver,
            Qadipose                  = Qadipose,
            Qbrain                    = Qbrain,
            #Qrest                     = Qrest,   # make sure to remove Qrest!!!!!
            
            # volume 
            Vart                      = Vart,
            Vgut                      = Vgut,
            Vliver                    = Vliver,
            Vlung                     = Vlung,
            Vbrain                    = Vbrain,
            Vadipose                  = Vadipose,
            Vrest                     = Vrest,
            Vven                      = Vven)  



Pred  <- function(pars, Qrest.new){
  
  parms_rest <- c(Qrest               = Qrest.new)
  parms      <- c(pars,parms_rest)
  
  # Modeling duration
  days          <- 1
  Time_max      <- days * 24                                  
  StartTime     <- 0                                        # Time_d
  StopTime      <- Time_max                                 
  dt            <- 0.1
  Times         <- seq(StartTime, StopTime, dt)
  #print(Times)
  # dose
  dose          <- 10                                       # mg/kg bw
  Oral_input    <- dose * 1000 / MW_parent                  # umol/kg/d, STUDY 2

  Dose_events    <- data.frame(var    = "Agutlumen_parent",
                               time   =  0,
                               value  = Oral_input,
                               method = "add")   
  
  df  <- ode(y       = initState, 
             times   = Times, 
             func    = pbpk8cpt, 
             parms   = parms,
             atol    = 1e-6, rtol = 1e-8,
             events  = list(data = Dose_events),
             method  = 'lsoda')
  #print(parms)

  df               <- as.data.frame(df)
  colnames(df)[1]  <- "Time"
  outdf  <- cbind.data.frame(Time = df$Time,
                             AUC  = max(df$AUC_Cblood_parent),
                             Cmax = max(df$C_blood_parent ))
  
  
  return (outdf)
}

#df <- Pred(parms)
NSC.AUC              <- matrix(nrow=length(parms),ncol=1)       # number of row is equal to number of parms
NSC.Cmax             <- matrix(nrow=length(parms),ncol=1)
rownames(NSC.AUC)    <- names(parms)
rownames(NSC.Cmax)   <- names(parms)

fold                 <- 1.1

for (i in 1: length(parms)){
  
  Qrest            <- parms[['Qcardiac']] - (parms[['Qliver']] + parms[['Qgut']] + parms[['Qadipose']] + parms[['Qbrain']])
  pars.changed     <- parms[i]  * fold 
  pars.rest        <- parms[-i]
  pars             <- c(pars.changed, pars.rest)    
  Qrest.new        <- pars[['Qcardiac']] - (pars[['Qliver']] + pars[['Qgut']] + pars[['Qadipose']] + pars[['Qbrain']])
  
  cat('i = ', i,  ', Pars: ',  names(parms[i]), ', Original value: ', parms[i], ', Changed: ',  pars.changed, '\n')
  cat('-------- Qrest = ', Qrest, '; Qrest.new = ', Qrest.new, '--------\n' )
  
  if(names(parms[i])!= 'Qrest'){
  delta.AUC        <- (Pred(pars, Qrest.new)$AUC  - Pred(parms, Qrest)$AUC)  / Pred(parms, Qrest)$AUC  / (fold-1) * 100
  delta.Cmax       <- (Pred(pars, Qrest.new)$Cmax - Pred(parms, Qrest)$Cmax) / Pred(parms, Qrest)$Cmax / (fold-1) * 100

  NSC.AUC[i,1]        <- tail(delta.AUC, n = 1)
  NSC.Cmax[i,1]       <- tail(delta.Cmax, n = 1)
  
  colnames(NSC.AUC)[1]   <- 'AUC'
  colnames(NSC.Cmax)[1]  <- 'Cmax'
  }
}

NSC               <- cbind(NSC.AUC, NSC.Cmax)
write.csv(NSC, paste('C:/XXX/Mammals/PPZ/SA/LocalSA_rat.csv'), row.names = TRUE)
#=============================================================================
############                End of PBPK model                    #############
#=============================================================================



