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
library(reshape2)    # Package for melt function to reshape the table
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

#===================================================================================================================
#                                             PBPK MODEL EQUATIONS                                                 #
#===================================================================================================================
########             load the data                ##########
#############################################################
# Theme Reference: https://ggplot2.tidyverse.org/reference/ggtheme.html
# in vivo data; 
path     <- "C:/xxx/OneDrive - Syngenta/AI/Propiconazole/PBK_Modeling/Reference compounds/Difenoconazole/Figures/DIFEN.xlsx"
obs      <- read_xlsx(path,  sheet = 'Sheet1')
obs      <- obs[obs$Species == 'CD1Mice'| obs$Species == 'C57BL/6Mice' | obs$Species == 'WTMice',]
obs$Time            <- as.numeric(obs$Time_h)
obs$C_blood_parent  <- as.numeric(obs$Conc_ugl)/MW_parent
obs                 <- na.omit(obs)

oral_15mg.kg        <- as.data.frame(obs[obs$Dose_mg_kg == 15,] %>% dplyr::select(Time, C_blood_parent))
oral_45mg.kg        <- as.data.frame(obs[obs$Dose_mg_kg == 45,] %>% dplyr::select(Time, C_blood_parent))
oral_150mg.kg       <- as.data.frame(obs[obs$Dose_mg_kg == 150,] %>% dplyr::select(Time, C_blood_parent))
iv_1mg.kg           <- as.data.frame(obs[obs$Method == 'iv',] %>% dplyr::select(Time, C_blood_parent))




Pred   <- function(pars){
  ###############################################################################
  # Shared parameters/states
  names(pars) <- names(pars)
  
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
  
  # shared modeling scenarios
  days          <- 7
  Time_max      <- days * 24                             # h
  StartTime     <- 0                                     # Time_min 
  StopTime      <- Time_max
  dt            <- 0.1
  Times         <- seq(StartTime, StopTime, dt)
  ndays         <- 6
  
  ###############################################################################
  ## Exposure scenario 1: 15 mg/kg/d oral
  Oral_45mg.kg           <- 45                            
  Oral_45mg.kg_input     <- Oral_45mg.kg  * 1000 / MW_parent                                # mg/kg -> umol/kg bw
  Dose_events_15mg.kg_oral <- data.frame(var    = "Agutlumen_parent",                       # Agutlumen_parent; Aven_parent
                                         time   = seq(from = 0, to = ndays * 24, by = 24),
                                         value  = rep(Oral_45mg.kg_input, ndays+1),
                                         method = "add")
  
  out_45mg.kg.oral             <- ode(y = initState, 
                                       times = Times, 
                                       func  = pbpk8cpt, 
                                       parms = pars,                      # Has to be 'pars' here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                       atol = 1e-8, rtol = 1e-10,
                                       events  = list(data = Dose_events_15mg.kg_oral),
                                       method = 'lsoda')
  
  out_45mg.kg.oral               <- as.data.frame(out_45mg.kg.oral)
  colnames(out_45mg.kg.oral)[1]  <- "Time"
  out_45mg.kg.oral               <- cbind.data.frame (Time             = out_45mg.kg.oral$Time,
                                                      C_blood_parent   = out_45mg.kg.oral$C_blood_parent)    # umol/L
  

  return(list("out_45mg.kg.oral "     = out_45mg.kg.oral))
}


## Cost fuction (from FME pckage) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  
  out               <-  Pred(pars)
  #print(out_5mgkg.oral.2005 )
  
  cost  <- modCost(model=out$out_45mg.kg.oral,   obs= oral_45mg.kg,  weight='std',   x="Time")               # weight std can only be used if there is more than 1 data point

  return(cost)
}

#==============================================
#######   tune for selected parameters  #######
#==============================================
pars  <- c(ka                    = 1,
           Km_hep_parent         = 4.9,      
           Vmax_hep_parent       = 63,     
           EC50                  = 2.97,     # uM
           Emax                  = 10 )  

system.time(Fit<- modFit(f = MCcost, 
                         p = pars, 
                         method  = "Marq" ,                             # or "Nelder-Mead"
                         lower   = c(0.5, 3,  30,  1.5, 6),            # range definement matters
                         upper   = c(3,   6,  70,  4.5,  18),  
                         control = nls.lm.control(nprint=1)))

summary(Fit)
Fit$par

'            ka   Km_hep_parent Vmax_hep_parent            EC50            Emax 
      0.6122179       3.1543364      69.5628047       1.9375712      17.7385217'

Sim.fit.out_45mgkg        <- Pred(Fit$par)$out_45mg.kg.oral  
df.Sim.fit.out_45mgkg     <- cbind.data.frame (Time=Sim.fit.out_45mgkg$Time, C_blood_parent   = Sim.fit.out_45mgkg$C_blood_parent)

# plot
ggplot() +
  geom_line (data = df.Sim.fit.out_45mgkg, aes(Time, (C_blood_parent)), col="#00AFBB", lwd=2) + 
  geom_point(data = oral_45mg.kg, aes(Time, C_blood_parent), col = 'red', size=2) + ylab("Conc (umol/L)") +
  xlab("Time (h)") + theme(text = element_text(size = 20)) + theme_bw(base_size = 14) #+ xlim(120, 160)
