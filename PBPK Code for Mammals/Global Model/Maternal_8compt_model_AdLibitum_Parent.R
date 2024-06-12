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


#===================================================================================================================
#                                                 Metabolism                                                       #
#===================================================================================================================
# General hepatic parameter
HPGL_human              <-  117.5                           # million cells/g liver; Number of hepacytes per gram liver; Acibenzolar manuscript appendix.
MPPGL_human             <-  45                              # mg/g;microsomal protein content per gram liver; Acibenzolar manuscript appendix
MPPGGI_human            <-  3                               # mg/g tissue; microsomal protein content in the GI. # not used: table 2 of TK0648566-01 report

HPGL_rat                <-  108                             # million cells/g liver; Number of hepacytes per gram liver; Acibenzolar manuscript appendix.
MPPGL_rat               <-  45                              # mg/g tissue;microsomal protein content per gram liver; Acibenzolar manuscript appendix
MPPGGI_rat              <-  3                               # mg/g tissue; microsomal protein content in the GI. table 2 of TK0648566-01 report

MPPGL_mouse             <-  45                              # mg/g tissue;Miyoung Yoon 2019;Sakai C 2014
MPPGGI_mouse            <-  3                               # mg/g tissue;

# Induction related parameter
kdeg_mouse             <- 0.0006 * 60                       # kdeg for CYP3A; h-1
kdeg_rat               <- 0.0005 * 60                       # kdeg for CYP3A
kdeg_human             <- 0.00032 * 60                      # kdeg for CYP3A

kt_mouse               <-  0.11                             # h-1
kt_rat                 <-  0.05                             # h-1
kt_rabbit              <-  0.04                             # h-1

BW_mouse    <- 0.02
BW_rat      <- 0.25
BW_human    <- 70
#===================================================================================================================
#                                             PBPK MODEL EQUATIONS                                                 #
#===================================================================================================================
# Reference
# PKNCA : https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk  : https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

# - Target compound : parent & daughter
# - Species         : Mammals (Maternal only)
# - Author          : Yaoxing Wu
# - Date            : Jun, 2022
# - Structure       : GI tract, Liver, Lungs, Adipose, Brain, Rest, Vein, Arterial
#===================================================================================================================

pbpk8cpt <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    #print(class(Agut_parent))
    #print(Agut_parent)
    ## chemical concentrations in tissue compartment
    Cgut_parent         =  Agut_parent / Vgut               # tissue concentration umol/L umolg/kg BW / L/kg BW
    Cliver_parent       =  Aliver_parent / Vliver           # tissue concentration umol/L
    Cbrain_parent       =  Abrain_parent / Vbrain           # tissue concentration umol/L
    Cadipose_parent     =  Aadipose_parent / Vadipose       # tissue concentration umol/L
    Crest_parent        =  Arest_parent / Vrest             # tissue concentration umol/L
    Cven_parent         =  Aven_parent / Vven               # blood concentration in vein umol/L
    Clung_parent        =  Alung_parent / Vlung             # tissue concentration umol/L
    Cart_parent         =  Aart_parent / Vart               # blood concentration in artery umol/L
    
    
    ## absorption (parent)
    #ka                               <- calc_ka(dose.mg.kg = Oral_mg.kg, Papp = Papp_parent)
    
    ## Clearance
    # hepatic, intestinal and plasma clearance rate L/h/BW for parent compound
    fu_hep_parent                    <- calc_hep_mic_fu(pH                   = pH,
                                                        media_conc_metabolic = if(parent_hep_conc_metabolic == 9999){parent_hep_conc_metabolic ='none'}else{parent_hep_conc_metabolic},  #parent_hep_conc_metabolic,
                                                        media_conc_binding   = if(parent_hep_conc_binding == 9999){parent_hep_conc_binding ='none'}else{media_conc_binding},  #parent_hep_conc_binding, 
                                                        fuinc_binding        = if(fuinc_hep_parent == 9999){fuinc_hep_parent ='none'}else{fuinc_binding},  #fuinc_hep_parent, 
                                                        Compound_type        = Compound_type_parent, 
                                                        LogP                 = LogP_parent, 
                                                        pKa                  = pKa_a_parent)
    
    fu_int_parent                    <- calc_hep_mic_fu(pH                   = pH, 
                                                        media_conc_metabolic = if(parent_int_mic_conc_metabolic == 9999){parent_int_mic_conc_metabolic ='none'}else{parent_int_mic_conc_metabolic},
                                                        media_conc_binding   = if(parent_int_mic_conc_binding == 9999){parent_int_mic_conc_binding ='none'}else{parent_int_mic_conc_binding},      #parent_int_mic_conc_binding, 
                                                        fuinc_binding        = if(fuinc_int_parent == 9999){fuinc_int_parent ='none'}else{fuinc_int_parent},                            #fuinc_hep_parent, 
                                                        Compound_type        = Compound_type_parent, 
                                                        LogP                 = LogP_parent, 
                                                        pKa                  = pKa_a_parent)
    
   
    
    # clearance
    Clint_hep_parent_invivo          <- calc_metabolic_clearance(type_clearance       = type_hep_clearance,
                                                                 incubation           = incubation_hep_parent, 
                                                                 fuinc                = fu_hep_parent,
                                                                 Vmax_unit            = Vmax_unit,
                                                                 Vmax = Vmax_hep_parent, Km = Km_hep_parent, 
                                                                 Clint_ori            = if(Clint_ori_hep_parent==9999){Clint_ori_hep_parent = 'none'}else{Clint_ori_hep_parent},
                                                                 C_tissue             = Cliver_parent, Ktissue2pu = Kliver2pu_parent,
                                                                 tissue_specific_volume = Vliver, 
                                                                 MPPG = MPPGL, fub = fub_parent, BW_scaledfrom = BW_scaledfrom) 
    
    Clint_int_parent_invivo          <- calc_metabolic_clearance(type_clearance       = type_int_clearance,
                                                                 incubation           = incubation_int_parent, 
                                                                 fuinc                = fu_hep_parent, 
                                                                 Vmax_unit            = Vmax_unit,
                                                                 Vmax = Vmax_int_parent, Km = Km_int_parent, 
                                                                 Clint_ori            = if(Clint_ori_int_parent==9999){Clint_ori_int_parent = 'none'}else{Clint_ori_int_parent},
                                                                 C_tissue             = Cgut_parent, Ktissue2pu = Kgut2pu_parent,
                                                                 tissue_specific_volume = Vgut, 
                                                                 MPPG = MPPGGI, fub = fub_parent, BW_scaledfrom = BW_scaledfrom) 
    
    Clint_plasma_parent_invivo       <- 0
    
    
   
    
    # Renal clearance of daughter compound
    Qgfr_parent                  <- 0# CLR
    
    
    #############          equations for different compartment (parent)              ##########
    RateC               <-  Cstep.doseinf(t)                                      # Infusion rate
    
    ## absorption (parent)
    Dose_mg.kg                      <- Dose4ka_df[round(Dose4ka_df$time,1) == round(t,1),]$value
    #print(Dose_mg.kg)
    ka                              <- calc_ka(dose.mg.kg = Dose_mg.kg, Papp = 'none')
    
    # Gutlumen
    dAgutlumen_parent    = - ka * Agutlumen_parent - kt * Agutlumen_parent                              # umol/h/kg Bw
    
    # Gut 
    Cgutblood_parent     =  Rblood2plasma_parent / (Kgut2pu_parent * fub_parent) * Cgut_parent 
    dAgut_parent         =  ka * fa * Agutlumen_parent + Qgut * (Cart_parent - Cgutblood_parent) - 
                                Clint_int_parent_invivo  * Cgut_parent / (Kgut2pu_parent * fub_parent)
    
    # Liver
    Cliverblood_parent   =   Rblood2plasma_parent / (Kliver2pu_parent * fub_parent) * Cliver_parent                         # L/h/kg * umol/L -> umol/h/kg bw
    dAliver_parent       =   Qliver * Cart_parent + Qgut * Cgutblood_parent  - (Qliver + Qgut) * Cliverblood_parent - 
                                Clint_hep_parent_invivo * Cliver_parent / Kliver2pu_parent                                  # L/h/kg BW * umol/L  --> umol/h/kg BW
    
    # Brain
    Cbrainblood_parent   =   Rblood2plasma_parent / (Kbrain2pu_parent * fub_parent) * Cbrain_parent
    dAbrain_parent       =   Qbrain * (Cart_parent - Cbrainblood_parent)
    
    # Adipose
    Cadiposeblood_parent  =   Rblood2plasma_parent / (Kadipose2pu_parent * fub_parent) * Cadipose_parent
    dAadipose_parent      =   Qadipose * (Cart_parent - Cadiposeblood_parent)
    
    # Rest of body
    Crestblood_parent    =   Rblood2plasma_parent / (Krest2pu_parent * fub_parent) * Crest_parent
    dArest_parent        =   Qrest * (Cart_parent - Crestblood_parent)
    
    # Venous blood
    dAven_parent         =   RateC + (Qliver + Qgut) * Cliverblood_parent  + Qbrain * Cbrainblood_parent + 
                                Qadipose * Cadiposeblood_parent + Qrest * Crestblood_parent - Qcardiac * Cven_parent - 
                                Clint_plasma_parent_invivo * Vven *  Cven_parent - Qgfr_parent * Cven_parent * fub_parent
    
    # Lung
    Clungblood_parent    =   Rblood2plasma_parent / (Klung2pu_parent * fub_parent) * Clung_parent
    dAlung_parent        =   Qcardiac * (Cven_parent - Clungblood_parent)
    
    # Artery
    dAart_parent         =   Qcardiac * (Clungblood_parent - Cart_parent) - Clint_plasma_parent_invivo  * Vart * Cart_parent 
    
    # AUC
    dAUC_Cplasma_parent  =  Cven_parent / Rblood2plasma_parent
    dAUC_Cblood_parent   =  Cven_parent
    
    dAurine_parent       =  Qgfr_parent * Cven_parent * fub_parent 
    
    ######################         Mass balance check          #######################
    ## Mass balance of aceibenzolar
    dAgut_parent_in       = ka * fa * Agutlumen_parent 
    dAven_parent_in       = RateC
    dAgut_parent_out      = kt * Agutlumen_parent 
    dAliver_parent_out    = Clint_hep_parent_invivo * Cliver_parent / Kliver2pu_parent
    dAven_parent_out      = Clint_plasma_parent_invivo * Vven *  Cven_parent
    dAint_parent_out      =  Clint_int_parent_invivo  * Cgut_parent / (Kgut2pu_parent * fub_parent)
    dAart_parent_out      = Clint_plasma_parent_invivo  * Vart * Cart_parent
    
    Mass_parent_in           = Agut_parent_in + Aven_parent_in                                      # Total Amount Acibenzolar from IV Dose and Absorbed in Gut; umol/kg bw
    Mass_parent_stored       = Agut_parent + Aliver_parent + Abrain_parent + Aadipose_parent +
                                  Arest_parent + Aven_parent + Alung_parent + Aart_parent           # Total Amount of Parent Remaining in the Body (??mol)
    Mass_parent_out          = Aliver_parent_out + Aven_parent_out + Aint_parent_out + Aart_parent_out           # Total Amount of Parent Excreted from the Body (??mol)
    Mass_parent_bal          = Mass_parent_in - Mass_parent_stored - Mass_parent_out
    
   
    ##### output parameters
    
    list(c(dAgutlumen_parent, 
           dAgut_parent, 
           dAliver_parent,
           dAbrain_parent,
           dAadipose_parent,
           dArest_parent,
           dAven_parent,
           dAlung_parent,
           dAart_parent,
           dAUC_Cplasma_parent, 
           dAUC_Cblood_parent,
           dAurine_parent,
           
           dAgut_parent_in, 
           dAven_parent_in,
           dAgut_parent_out,
           dAliver_parent_out,
           dAven_parent_out,
           dAint_parent_out,
           dAart_parent_out
           
    ),
    
    "C_blood_parent"     = Cven_parent,
    "C_plasma_parent"    = Cven_parent / Rblood2plasma_parent, 
    "C_brain_parent"     = Cbrain_parent,
    "C_liver_parent"     = Cliver_parent,
    "C_adipose_parent"   = Cadipose_parent,
    'C_lung_parent'      = Clung_parent,
    'C_gut_parent'       = Cgut_parent,

    'Mass_parent_in'        = Mass_parent_in ,   
    'Mass_parent_stored'    = Mass_parent_stored ,    
    'Mass_parent_out'       = Mass_parent_out ,         
    'Mass_parent_bal'       = Mass_parent_bal  
    ) 
  })
}


#-------------------------          End of PBPK modeling equations        --------------------------------  


