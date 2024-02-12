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

kt_human               <-  0.02                             # h-1
kt_mouse               <-  0.11                             # h-1
kt_rat                 <-  0.054                            # h-1

BW_mouse    <- 0.02
BW_rat      <- 0.25
BW_human    <- 70

#===================================================================================================================
#                                             PBPK MODEL EQUATIONS                                                 #
#===================================================================================================================
# Reference
# PKNCA : https://cran.r-project.org/web/packages/PKNCA/vignettes/AUC-Calculation-with-PKNCA.html
# httk  : https://github.com/USEPA/CompTox-ExpoCast-httk/tree/main/httk

# - Target compound : PPZ
# - Species         : Human (Maternal only)
# - Author          : Yaoxing Wu
# - Date            : Jun, 2022
# - Structure       : GI tract, Liver, Lungs, Adipose, Brain, Rest, Vein, Arterial
#===================================================================================================================

pbpk8cpt <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    #################  function ########################
    ## blood
    if(Rblood2plasma_parent == 'None'){
      Rblood2plasma_parent <- calc_rblood2plasma(chem.cas  = CAS_parent,
                                                 species   = species,
                                                 adjusted.Funbound.plasma = F)
    }else{
      Rblood2plasma_parent <- Rblood2plasma_parent
    }
    
    
    calc_hep_mic_fu        <- function(pH, media_conc_metabolic, media_conc_binding,fuinc_binding,
                                       Compound_type, LogP, pKa){     
      if(tolower(Compound_type) == 'acid'){
        LogP       <- LogP + log10(1/(1+10^(pH-pKa)))           # logD calculation for acid
      }else{ 
        LogP       <- LogP
      }                                     # logP for bases and neutrals 
      
      if(tolower(fuinc_binding) != 'none'){
        fuinc           <- 1 / ( media_conc_metabolic / media_conc_binding * ((1 - fuinc_binding) / fuinc_binding) + 1)     
      }else{
        Kmic            <- 10^(0.072*(LogP)^2+0.067*(LogP)-1.126)
        fuinc           <- 1 / ( 125*0.005*Kmic+1)              # calc_hep_fu(); Ratio of cell volume to incubation volume. Default (0.005) 
      }
      
      return(fuinc)
    }
    
    # hepatic and intestinal clearance for parent; hepatic clearance for daughter (rat only)
    calc_metabolic_clearance  <- function(type_clearance, incubation, fuinc, Vmax_unit,       # if media conc not available, use 1 instead for both
                                          Vmax, Km, Clint_ori, 
                                          C_tissue, Ktissue2pu, tissue_specific_volume, 
                                          MPPG, fub, BW_scaledfrom){
      
      if(tolower(type_clearance) == 'liver'){
        
        if(is.numeric(Vmax)){
          if(Vmax_unit == 'umol/h/kg bw'){
            if(incubation == 'scaled'){
              Vmax_scaledfrom   <- Vmax * BW_scaledfrom                                                            # umol/h; volume per unit time
              Vmax_scaled       <- Vmax_scaledfrom * ((BW / BW_scaledfrom)^0.75) /BW                               # umol/h/kg bw; 0.67 or 0.75
              Clint             <- Vmax_scaled /(Km * fuinc + C_tissue / Ktissue2pu) 
            }else{
              Clint        <- Vmax /(Km * fuinc + C_tissue / Ktissue2pu) 
            }
          }else{
            # Vmax is in the unit of umol/h/mg protein or umol/h/million cells; Km is in the unit of um (umol/L); C_tissue is in the unit of umol/L.
            Clint          <- Vmax /(Km * fuinc + C_tissue / Ktissue2pu) * MPPG * (tissue_specific_volume * 1.05 * 1000)  # unit: L/h/kg BW
          }
        }else if(is.numeric(Clint_ori)){
          if(incubation == 'scaled'){
            # # Clint_scaled from is in the unit of h-1
            Clint_scaledfrom   <- Clint_ori * BW_scaledfrom
            Clint              <- Clint_scaledfrom * ((BW / BW_scaledfrom)^0.75) /BW       
          }else{
            # Clint_ori is in the unit of uL/h/million cell
            Clint          <- Clint_ori / fuinc * 1e-6 * MPPG * tissue_specific_volume * 1.05 * 1000              # L/h/million cells * million cells/g tissue * g/kg  -->  L/h/kg BW
          }
        }
        
      }else if(tolower(type_clearance) == 'intestine'){
        if(is.numeric(Vmax)){
          if(Vmax_unit == 'umol/h/kg bw'){
            if(incubation == 'scaled'){
              Vmax_scaledfrom   <- Vmax * BW_scaledfrom                                                            # umol/h; volume per unit time
              Vmax_scaled       <- Vmax_scaledfrom * ((BW / BW_scaledfrom)^0.75) /BW                               # umol/h/kg bw; 0.67 or 0.75
              Clint             <- Vmax_scaled /(Km * fuinc + C_tissue / (Ktissue2pu * fub))
            }else{
              Clint        <- Vmax /(Km * fuinc + C_tissue / (Ktissue2pu * fub)) 
            }
          }else{
            # Vmax is in the unit of umol/h/mg protein or umol/h/million cells; Km is in the unit of um (umol/L); C_tissue is in the unit of umol/L.
            Clint          <- Vmax /(Km * fuinc + C_tissue / (Ktissue2pu * fub)) * MPPG * (tissue_specific_volume * 1.05 * 1000)                # unit: L/h/kg BW
          }
        }else if(is.numeric(Clint_ori)){
          if(incubation == 'scaled'){
            # # Clint_scaled from is in the unit of h-1
            Clint_scaledfrom   <- Clint_ori * BW_scaledfrom
            Clint              <- Clint_scaledfrom * ((BW / BW_scaledfrom)^0.75) /BW   
          }else if(incubation == 'half-life') {
            # Clint_ori is in the unit of h (half-life)
            # Clint           <- log(2) * 0.2 / 1000 / Clint_ori                                                   # case by case; check the unit; equation A14 in appenxix
            Clint           <-  log(2) / Clint_ori                                                                # h-1         
          }else{
            # Clint_ori is in the unit of uL/h/million cell or uL/h/mg protein
            Clint          <- Clint_ori / fuinc * 1e-6 * MPPG * tissue_specific_volume * 1.05 * 1000              # L/h/million cells * million cells/g tissue * g/kg  -->  L/h/kg BW
          }
        }
      }
      
      return(Clint)                                                                                          # unit: (L/h/kg BW)
    }
    
    
    
    
    
    
    #################     unanalyzed parameters needs to be included in the model    ###################
    BW_mouse    <- 0.02
    BW_rat      <- 0.25
    BW_human    <- 70
    
    Compound_name_parent    <- 'Propiconazole (PPZ)'
    Compound_type_parent    <- 'base'
    CAS_parent              <- "60207-90-1-a"                                # add a to CAS since PPZ has been included in HTTK database
    LogP_parent             <-  3.72                                         # log Kow; from raav
    MW_parent               <-  342.2
    pKa_a_parent            <-  'None'                                       # very weak base without any numerical value; chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://hal-unilim.archives-ouvertes.fr/hal-00929974/document
    pKa_b_parent            <-  1.09                                         # pka accept (base); Teb?H+ -> Teb + H+ 

    #Papp_parent             <-  40.8  
    
    BW_scaledfrom           <- BW_rat
    Vmax_unit               <- 'umol/h/kg bw'                         # 'umol/h/kg bw' or 'none'
    type_hep_clearance      <- 'liver'
    incubation_hep_parent   <- 'microsome'                            # for rat metabolic stability test (hepatic)
        
    E0                      <- 1 
    MPPGL_human             <- 45
    MPPGGI_human            <- 3
    fa                      <- 1
    
    Qcardiac                  = 4.799987
    Qgut                      = 0.9857191
    Qheart                    = 0.2057108
    Qkidney                   = 1.062884
    Qliver                    = 1.242935
    Qadipose                  = 0.2227825
    Qbrain                    = 0.6001021
    Qlung                     = 0.1199997
    Qmuscle                   = 0.6428332
    Qskin                     = 0.2572163
    Qrest                     = 1.748448
    
    Vmax_unit                        <- 'umol/h/kg bw'                   # 'umol/h/kg bw' or 'none'
    type_hep_clearance               <- 'liver'
    incubation_hep_parent            <- 'microsome'                            # for rat metabolic stability test (hepatic)
    fuinc_hep_parent                 <- 9999                                  # binding was not measured in in vitro metabolsim assays
    parent_hep_conc_binding          <- 9999
    Clint_ori_hep_parent             <- 9999
    
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
                                                        media_conc_metabolic = if(parent_hep_conc_metabolic == 9999){parent_hep_conc_metabolic ='none'},  #parent_hep_conc_metabolic,
                                                        media_conc_binding   = if(parent_hep_conc_binding == 9999){parent_hep_conc_binding ='none'},  #parent_hep_conc_binding, 
                                                        fuinc_binding        = if(fuinc_hep_parent == 9999){fuinc_hep_parent ='none'},  #fuinc_hep_parent, 
                                                        Compound_type        = Compound_type_parent, 
                                                        LogP                 = LogP_parent, 
                                                        pKa                  = pKa_a_parent)
    
    fu_int_parent                    <- calc_hep_mic_fu(pH                   = pH, 
                                                        media_conc_metabolic = if(parent_int_mic_conc_metabolic == 9999){parent_int_mic_conc_metabolic ='none'},
                                                        media_conc_binding   = if(parent_int_mic_conc_binding == 9999){parent_int_mic_conc_binding ='none'},      #parent_int_mic_conc_binding, 
                                                        fuinc_binding        = if(fuinc_int_parent == 9999){fuinc_int_parent ='none'},                            #fuinc_hep_parent, 
                                                        Compound_type        = Compound_type_parent, 
                                                        LogP                 = LogP_parent, 
                                                        pKa                  = pKa_a_parent)
    
   
    
    # clearance
    Clint_hep_parent_invivo          <- calc_metabolic_clearance(type_clearance       = type_hep_clearance,
                                                                 incubation           = incubation_hep_parent, 
                                                                 fuinc                = fu_hep_parent,
                                                                 Vmax_unit            = Vmax_unit,
                                                                 Vmax = fmCYP * Vmax_hep_parent * Eliver + (1 - fmCYP) * Vmax_hep_parent, Km = Km_hep_parent, 
                                                                 Clint_ori            = if(Clint_ori_hep_parent==9999){Clint_ori_hep_parent = 'none'},
                                                                 C_tissue             = Cliver_parent, Ktissue2pu = Kliver2pu_parent,
                                                                 tissue_specific_volume = Vliver, 
                                                                 MPPG = MPPGL, fub = fub_parent, BW_scaledfrom = BW_scaledfrom) 
    
    Clint_int_parent_invivo          <- calc_metabolic_clearance(type_clearance       = type_int_clearance,
                                                                 incubation           = incubation_int_parent, 
                                                                 fuinc                = fu_hep_parent, 
                                                                 Vmax_unit            = Vmax_unit,
                                                                 Vmax = Vmax_int_parent * Egut, Km = Km_int_parent, 
                                                                 Clint_ori            = "None",
                                                                 C_tissue             = Cgut_parent, Ktissue2pu = Kgut2pu_parent,
                                                                 tissue_specific_volume = Vgut, 
                                                                 MPPG = MPPGGI, fub = fub_parent, BW_scaledfrom = BW_scaledfrom) 
    
    Clint_plasma_parent_invivo       <- 0
    
    
   
    
    # Renal clearance of daughter compound
    Qgfr_parent                  <- 0# CLR
    
    
    #############          equations for different compartment (parent)              ##########

    # Gutlumen
    dAgutlumen_parent    = - ka * Agutlumen_parent - kt * Agutlumen_parent                              # umol/h/kg Bw
    
    # Gut 
    Cgutblood_parent     =  Rblood2plasma_parent / (Kgut2pu_parent * fub_parent) * Cgut_parent 
    dAgut_parent         =  ka * fa * Agutlumen_parent + Qgut * (Cart_parent - Cgutblood_parent) - 
                                Clint_int_parent_invivo  * Cgut_parent / (Kgut2pu_parent * fub_parent)
    #########  Gut CYP3A Enzyme  #########
    Igut                 =  Cgut_parent / (Kgut2pu_parent * fub_parent)
    dEgut                =  kdeg * E0 + kdeg * E0 * Emax * Igut /(EC50 * fu_hep_parent + Igut) - kdeg * Egut
    
    # Liver
    Cliverblood_parent   =   Rblood2plasma_parent / (Kliver2pu_parent * fub_parent) * Cliver_parent                         # L/h/kg * umol/L -> umol/h/kg bw
    dAliver_parent       =   Qliver * Cart_parent + Qgut * Cgutblood_parent  - (Qliver + Qgut) * Cliverblood_parent - 
                                Clint_hep_parent_invivo * Cliver_parent / Kliver2pu_parent                                  # L/h/kg BW * umol/L  --> umol/h/kg BW
    #########  Liver CYP3A Enzyme  #########
    Iliver               =  Cliver_parent / (Kliver2pu_parent * fub_parent)
    #dEliver              =  kdeg * E0 + kdeg * E0 * Emax * Iliver * fub_parent /(EC50 * fu_hep_parent + Iliver * fub_parent) - kdeg * Eliver
    dEliver              =  kdeg * E0 + kdeg * E0 * Emax * Iliver / (EC50 * fu_hep_parent + Iliver) - kdeg * Eliver
    
    # Brain
    Cbrainblood_parent   =   Rblood2plasma_parent / (Kbrain2pu_parent * fub_parent) * Cbrain_parent
    dAbrain_parent       =   Qbrain * (Cart_parent - Cbrainblood_parent)
    
    # Adipose
    Cadiposeblood_parent  =  Rblood2plasma_parent / (Kadipose2pu_parent * fub_parent) * Cadipose_parent
    dAadipose_parent      =  Qadipose * (Cart_parent - Cadiposeblood_parent)
    
    # Rest of body
    Crestblood_parent    =   Rblood2plasma_parent / (Krest2pu_parent * fub_parent) * Crest_parent
    dArest_parent        =   Qrest * (Cart_parent - Crestblood_parent)
    
    # Venous blood
    dAven_parent         =    (Qliver + Qgut) * Cliverblood_parent  + Qbrain * Cbrainblood_parent + 
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
    ## Mass balance of parent
    dAgut_parent_in       = ka * fa * Agutlumen_parent 
    dAven_parent_in       = 0
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
           dEgut,
           dAliver_parent,
           dEliver,
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
    "C_gut_parent"       = Cgut_parent,
    "C_liver_parent"     = Cliver_parent,
    "C_adipose_parent"   = Cadipose_parent,
    'C_lung_parent'      = Clung_parent,
    'Iliver'             = Iliver,

    'Mass_parent_in'        = Mass_parent_in ,   
    'Mass_parent_stored'    = Mass_parent_stored ,    
    'Mass_parent_out'       = Mass_parent_out ,         
    'Mass_parent_bal'       = Mass_parent_bal  
    ) 
  })
}


#-------------------------          End of PBPK modeling equations        --------------------------------  


