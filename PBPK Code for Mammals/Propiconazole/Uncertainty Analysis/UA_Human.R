library(FME)
library(invgamma) 
library(ggplot2)
library(httk)
library(deSolve)
library(pksensi)
library(ggplot2)
library(PKNCA)
library(sensitivity)
library(ODEsensitivity)
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(tidyr)
library(tidyverse)
library(dplyr)
library(minpack.lm)  # Package for model fitting


#=================================================================================================
##########                            Uncertainty analysis                              ##########
#=================================================================================================
SA_mouse      <- read.csv('C:/XXX/Mammals/PPZ/SA/LocalSA_Mouse.csv')
SA_rat        <- read.csv('C:/XXX/Mammals/PPZ/SA/LocalSA_rat.csv')
SA_human      <- read.csv('C:/XXX/PPZ/SA/LocalSA_Human.csv')

SA_mouse$Species        <- 'Mouse'
SA_rat$Species          <- 'Rat'
SA_human$Species        <- 'Human'

SA                    <- as.data.frame(rbind(SA_mouse, SA_rat, SA_human))
SA                    <- unique(SA)
SA$Species            <- factor(SA$Species, levels=c('Mouse', 'Rat', 'Human'))
colnames(SA)[1]       <- 'parameter'
SA                    <- SA[SA$parameter != 'Qrest',]

SA$parameter <- sub("_parent", "", SA$parameter)
SA$parameter <- sub("fub", "fup", SA$parameter)

SA$AUC                <- (SA$AUC)/100
SA$Cmax               <- (SA$Cmax)/100

SA_AUC                <- subset(SA, select = - abs(Cmax))
SA_Cmax               <- subset(SA, select = - abs(AUC))
SA_AUC                <- SA_AUC[abs(SA_AUC$AUC) >= 0.1, ]
SA_Cmax               <- SA_Cmax[abs(SA_Cmax$Cmax) >= 0.1, ]

SA_AUC_human          <- SA_AUC[SA_AUC$Species == 'Human',]
SA_AUC_human  <- SA_AUC_human  %>%
  mutate(group = case_when(
    abs(AUC) >=0.5 ~ 'high',
    abs(AUC) < 0.5 & abs(AUC) >= 0.2 ~ 'medium',
    abs(AUC) < 0.2 ~ 'low'
  ))

print(SA_AUC_human$parameter)

library(truncnorm) 
library(ggridges)     # Used to create Figure
library(ggplot2)   # ggplot is the basic package for creating plots.
#install.packages("tkrplot")
library(tkrplot)
#install.packages("rriskDistributions")
library(rriskDistributions)
library(EnvStats)
library(httk)
library(mc2d)
library(fitdistrplus)

#rm(list = ls())
set.seed(10) 

# Normal distributions were implemented for physiological parameters; 
# chemical-specific parameters were assumed to be log-normally distributed (Henri et al., 2017; Li et al., 2017; Yang et al., 2015).

################################################################################################
#####################  NOW run local sensitivity file (human) until line 278  ##################
################################################################################################
#dev.off()

species                 <- 'Human'  

source(paste("C:/XXX/General_code/Maternal_8compt_InductionModel_Parent.R"))
source(paste("C:/XXX/General_code/Species_PhyData.R"))
MPPGL      <- MPPGL_human
MPPGGI     <- MPPGGI_human
kdeg       <- kdeg_human

BW_mouse    <- 0.02
BW_rat      <- 0.25
BW_human    <- 70


#### mean value of parameters
Rblood2plasma_parent.mean          <- 0.597
fub_parent.mean                    <- 0.05 
Emax.mean                          <- 5.62
Vmax_hep_parent.mean               <- 71.3
Km_hep_parent.mean                 <- 4.6
Vmax_int_parent.mean               <- 4.753333
Km_int_parent.mean                 <- 4.6
ka.mean                            <- 1#0.3732544
kdeg.mean                          <- 0.0192
Qgut.mean                          <- 0.9857191
Qadipose.mean                      <- 0.2227825

#### std value of parameters

Rblood2plasma_parent.sd      <- 0.3 * Rblood2plasma_parent.mean
fub_parent.sd                <- 0.3 * fub_parent.mean
Vmax_hep_parent.sd           <- 0.7 * Vmax_hep_parent
Km_hep_parent.sd             <- 0.5 * Km_hep_parent.mean 
Vmax_int_parent.sd           <- 0.7 * Vmax_int_parent.mean
Km_int_parent.sd             <- 0.5 * Km_int_parent.mean 
Emax.sd                      <- 0.5 * Emax.mean
ka.sd                        <- 0.5 * ka.mean 
kdeg.sd                      <- 0.5 * kdeg.mean 
Qgut.sd                      <- 0.3 * Qgut.mean
Qadipose.sd                  <- 0.3 * Qadipose.mean




#############     define distribution
# Normal distributions were implemented for physiological parameters; 
m.log.Rblood2plasma_parent             <- log(Rblood2plasma_parent.mean^2/(Rblood2plasma_parent.sd^2+Rblood2plasma_parent.mean^2)^0.5) 
sd.log.Rblood2plasma_parent            <- (log(1+Rblood2plasma_parent.sd^2/Rblood2plasma_parent.mean^2))^0.5 

m.log.fub_parent             <- log(fub_parent.mean^2/(fub_parent.sd^2+fub_parent.mean^2)^0.5) 
sd.log.fub_parent            <- (log(1+fub_parent.sd^2/fub_parent.mean^2))^0.5 

m.log.Vmax_hep_parent        <- log(Vmax_hep_parent.mean^2/(Vmax_hep_parent.sd^2+Vmax_hep_parent.mean^2)^0.5)
sd.log.Vmax_hep_parent       <- (log(1+Vmax_hep_parent.sd^2/Vmax_hep_parent.mean^2))^0.5

m.log.Km_hep_parent          <- log(Km_hep_parent.mean^2/(Km_hep_parent.sd^2+Km_hep_parent.mean^2)^0.5)
sd.log.Km_hep_parent         <- (log(1+Km_hep_parent.sd^2/Km_hep_parent.mean^2))^0.5

m.log.Vmax_int_parent        <- log(Vmax_int_parent.mean^2/(Vmax_int_parent.sd^2+Vmax_int_parent.mean^2)^0.5)
sd.log.Vmax_int_parent       <- (log(1+Vmax_int_parent.sd^2/Vmax_int_parent.mean^2))^0.5

m.log.Km_int_parent          <- log(Km_int_parent.mean^2/(Km_int_parent.sd^2+Km_int_parent.mean^2)^0.5)
sd.log.Km_int_parent         <- (log(1+Km_int_parent.sd^2/Km_int_parent.mean^2))^0.5

m.log.Emax          <- log(Emax.mean^2/(Emax.sd^2+Emax.mean^2)^0.5)
sd.log.Emax         <- (log(1+Emax.sd^2/Emax.mean^2))^0.5

m.log.ka           <- log(ka.mean^2/(ka.sd^2+ka.mean^2)^0.5)
sd.log.ka          <- (log(1+ka.sd^2/ka.mean^2))^0.5

####################
N=1000
set.seed(10) 

idata <- 
  tibble(ID=1:N) %>% 
  mutate(
    Qgut = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = Qgut.mean, sd = Qgut.sd),
      b = qnorm(0.975, mean = Qgut.mean, sd = Qgut.sd),
      mean = Qgut.mean,
      sd = Qgut.sd
    ),
    
    Qadipose = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = Qadipose.mean, sd = Qadipose.sd),
      b = qnorm(0.975, mean = Qadipose.mean, sd = Qadipose.sd),
      mean = Qadipose.mean,
      sd = Qadipose.sd
    ),
    
    ## Chemical properties
    Vmax_hep_parent = rlnormTrunc(
      N,
      meanlog = m.log.Vmax_hep_parent,
      sdlog = sd.log.Vmax_hep_parent,
      min = qlnorm(0.025, meanlog = m.log.Vmax_hep_parent, sdlog = sd.log.Vmax_hep_parent),
      max = qlnorm(0.975, meanlog = m.log.Vmax_hep_parent, sdlog = sd.log.Vmax_hep_parent)
    ),
    Km_hep_parent = rlnormTrunc(
      N,
      meanlog = m.log.Km_hep_parent,
      sdlog = sd.log.Km_hep_parent,
      min = qlnorm(0.025, meanlog = m.log.Km_hep_parent, sdlog = sd.log.Km_hep_parent),
      max = qlnorm(0.975, meanlog = m.log.Km_hep_parent, sdlog = sd.log.Km_hep_parent)
    ),
    
    Vmax_int_parent = rlnormTrunc(
      N,
      meanlog = m.log.Vmax_int_parent,
      sdlog = sd.log.Vmax_int_parent,
      min = qlnorm(0.025, meanlog = m.log.Vmax_int_parent, sdlog = sd.log.Vmax_int_parent),
      max = qlnorm(0.975, meanlog = m.log.Vmax_int_parent, sdlog = sd.log.Vmax_int_parent)
    ),
    
    Km_int_parent = rlnormTrunc(
      N,
      meanlog = m.log.Km_int_parent,
      sdlog = sd.log.Km_int_parent,
      min = qlnorm(0.025, meanlog = m.log.Km_int_parent, sdlog = sd.log.Km_int_parent),
      max = qlnorm(0.975, meanlog = m.log.Km_int_parent, sdlog = sd.log.Km_int_parent)
    ),
    
    Rblood2plasma_parent = rlnormTrunc(
      N,
      meanlog = m.log.Rblood2plasma_parent,
      sdlog = sd.log.Rblood2plasma_parent,
      min = qlnorm(0.025, meanlog = m.log.Rblood2plasma_parent, sdlog = sd.log.Rblood2plasma_parent),
      max = qlnorm(0.975, meanlog = m.log.Rblood2plasma_parent, sdlog = sd.log.Rblood2plasma_parent)
    ),
    
    ka = rlnormTrunc(
      N,
      meanlog = m.log.ka,
      sdlog = sd.log.ka,
      min = qlnorm(0.025, meanlog = m.log.ka, sdlog = sd.log.ka),
      max = qlnorm(0.975, meanlog = m.log.ka, sdlog = sd.log.ka)
    ), 
    
    fub_parent = rlnormTrunc(
      N,
      meanlog = m.log.fub_parent,
      sdlog = sd.log.fub_parent,
      min = qlnorm(0.025, meanlog = m.log.fub_parent, sdlog = sd.log.fub_parent),
      max = 1
    )
    
  
    
  )


N_refined  <- nrow(idata)
N_refined 



###########   Single dose of 66 mg/kg BW  #########
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


#============================================================
######             Single dose population               #####
#============================================================
pars  <- c( Rblood2plasma_parent      = Rblood2plasma_parent,
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
            Vmax_hep_parent            = Vmax_hep_parent,                                       # unit: nmol/min/10^6 cells; Table 1 -->   umol/h/10^6 cellS; determined at 1000 cells/mL; page 13/90 of frost 2017
            Km_hep_parent              = Km_hep_parent,                                         # unit: uM (umol/L)
            Vmax_int_parent            = Vmax_int_parent ,                                           # unit: nmol/min/mg protein; Table 1 -->   umol/h/mg protein; Intestinal metabolism determined at 50 ug/mL Frost 2017
            Km_int_parent              = Km_int_parent,                                             # unit: uM (umol/L)
            CLR                        = CLR,
            
            # Parametes for flow 
            Qcardiac                  = Qcardiac, 
            Qgut                      = Qgut,
            Qliver                    = Qliver,
            Qadipose                  = Qadipose,
            Qbrain                    = Qbrain,
            Qrest                     = Qrest,   
            
            # volume 
            Vart                      = Vart,
            Vgut                      = Vgut,
            Vliver                    = Vliver,
            Vlung                     = Vlung,
            Vbrain                    = Vbrain,
            Vadipose                  = Vadipose,
            Vrest                     = Vrest,
            Vven                      = Vven)  

Pred     <- function(pars){
  
  names(pars) <- names(pars)
  
  # Modeling duration
  days          <- 1
  Time_max      <- days * 24                                    
  StartTime     <- 0                                        # Time_d
  StopTime      <- Time_max                                 
  dt            <- 0.1
  Times         <- seq(StartTime, StopTime, dt)
  # dose
  dose          <- 10                                       # mg/kg bw
  Oral_input    <- dose * 1000 / MW_parent                  # umol/kg/d, STUDY 2
  
  Dose_events    <- data.frame(var    = "Agutlumen_parent",
                               time   =  0,
                               value  = Oral_input,
                               method = "add")   
  
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
  
  
  out                         <- ode(y       = initState, 
                                     times   = Times, 
                                     func    = pbpk8cpt, 
                                     parms   = parms,
                                     atol    = 1e-6, rtol = 1e-8,
                                     events  = list(data = Dose_events),
                                     method  = 'lsoda')  
  
  out                         <- as.data.frame(out)
  colnames(out)[1]            <- "Time"
  out                         <- cbind.data.frame(Time = out$Time,
                                                  AUC  = (out$AUC_Cblood_parent),
                                                  Cmax = max(out$C_blood_parent ))       # ug/L 
  
  return(list("out"  = out))
}



Newtime.r         = Pred(pars)$out$Time
nrwo.r            = length (Newtime.r)

# Create the matrix 


print(SA_AUC_human$parameter)

UA_list  <- c('Rblood2plasma_parent',
              'fub_parent', 
              'Vmax_hep_parent',
              'Km_hep_parent',
              'Vmax_int_parent',
              'Km_int_parent',
              'ka',
              'Qgut',
              'Qadipose')



MC.AUC     = matrix(nrow = N_refined , ncol = length(UA_list))
MC.Cmax    = matrix(nrow = N_refined,  ncol = length(UA_list))
colnames(MC.AUC)   <- UA_list
colnames(MC.Cmax)  <- UA_list

for (j in seq(UA_list)){
  
  p = UA_list[j]
  
  for (i in 1:N_refined){
    
    cat("Parameter = ", p, "; iteration = ", i , "\n")
    
    new_pars       <- pars
    pars_UA        <- idata[i,]%>% dplyr :: select(-ID)
    
    if( p == 'Qadipose'){
      Qrest                 <-  Qcardiac - Qliver - Qgut - pars_UA[[p]] - Qbrain
      new_pars[['Qrest']]   <- Qrest
    }
    if( p == 'Qgut'){
      Qrest                 <-  Qcardiac - Qliver - Qadipose - pars_UA[[p]] - Qbrain
      new_pars[['Qrest']]   <- Qrest
    }
    
    new_pars[[p]]  <- pars_UA[[p]]
    
    
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
    
    days          <- 1
    Time_max      <- days * 24                                    
    StartTime     <- 0                                        # Time_d
    StopTime      <- Time_max                                 
    dt            <- 0.1
    Times         <- seq(StartTime, StopTime, dt)
    
    dose          <- 10                                       # mg/kg bw
    Oral_input    <- dose * 1000 / MW_parent                  # umol/kg/d, STUDY 2
    
    Dose_events    <- data.frame(var    = "Agutlumen_parent",
                                 time   =  0,
                                 value  = Oral_input,
                                 method = "add")   
    
    
    
    df     <- ode(y       = initState, 
                  times   = Times, 
                  func    = pbpk8cpt, 
                  parms   = new_pars,
                  atol    = 1e-6, rtol = 1e-8,
                  events  = list(data = Dose_events),
                  method  = 'lsoda')
    
    df                   <- as.data.frame(df)
    colnames(df)[1]      <- "Time"
    
    MC.AUC[i,j]     <- max(df$AUC_Cblood_parent)
    MC.Cmax[i,j]    <- max(df$C_blood_parent)
    
    
  }
}

# Calculate the 95th percentile for each column
percentile_95 <- apply(MC.AUC, 2, function(x) quantile(x, probs = 0.95))

# Calculate the median for each column
medians <- apply(MC.AUC, 2, median)

# Create a new data frame with the 95th percentile and median
new_df <- rbind(percentile_95, medians)

stat  <- t(new_df)
stat   <- as.data.frame(stat)
stat$ratio  <- round((stat$percentile_95 /stat$medians), 2)

stat <- stat %>%
  mutate(group = case_when(
    ratio >=2 ~ 'high',
    ratio < 2 & ratio >= 1.3 ~ 'medium',
    ratio < 1.3 ~ 'low'
  ))

# end



write.csv(SA_AUC_human, "C:/XXX/PPZ/UA/Results/SA_AUC_human.csv")
write.csv(stat, "C:/XXX/UA/Results/UA_AUC_human.csv")








