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
library(janitor)
library(gridExtra)
library(deSolve)
library(pksensi)
library(ggplot2)
library(sensitivity)
library(ODEsensitivity)
library('stringr')
library(pracma)
library(janitor)
library(gridExtra)

# Clean the environment
rm(list = ls())
#dev.off()

species                 <- 'Rat'  
Compound_name_parent    <- 'Propiconazole (PPZ)'
Compound_type_parent    <- 'base'
CAS_parent              <- "60207-90-1-a"                                # add a to CAS since PPZ has been included in HTTK database
LogP_parent             <-  3.72                                         # log Kow; from raav
MW_parent               <-  342.2
pKa_a_parent            <-  'None'                                       # very weak base without any numerical value; chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://hal-unilim.archives-ouvertes.fr/hal-00929974/document
pKa_b_parent            <-  1.09                                         # pka accept (base); Teb?H+ -> Teb + H+ 
fub_parent              <-  0.05
Rblood2plasma_parent    <-  0.597


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


source(paste("C:/XXX/SA/SA_Maternal_8compt_InductionModel_Rat.R"))
source(paste("C:/XXX/Species_PhyData.R"))
source(paste("C:/XXX/Maternal_DM_Parent.R"))

BW_scaledfrom    <- BW_rat
MPPGL      <- MPPGL_rat
MPPGGI     <- MPPGGI_rat
kdeg       <- kdeg_rat

BW_mouse    <- 0.02
BW_rat      <- 0.25
BW_human    <- 70


oral_mg.kg    <- 10                                        # 
oral_input    <- oral_mg.kg  * 1000  / MW_parent          # mg/kg -> ug/kg bw

initState <-  c(Agutlumen_parent    = oral_input, 
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



#=========================================================================================================================
##########################                           sensitivity                            ##############################
#=========================================================================================================================
# https://nanhung.rbind.io/project/fda-pbpk/
# https://rstudio-pubs-static.s3.amazonaws.com/500619_2a8dc150b1b744088bd651fd3e2d93fe.html
# https://rpubs.com/Nanhung/Suppl_190516
## http://rstudio-pubs-static.s3.amazonaws.com/498762_32e363e09a714c339d047fccbc6d039b.html  ##

path_sens <- "C:/XXX/Mammals/PPZ/SA/"
time_int  <- 1           # % Time points to examine
duration  <- seq(1, time_int, by = 1)
p         <- vector('list', time_int)

LVpars    <- c( 'Km_hep_parent',                  
                'Vmax_hep_parent',
                
                'ka',
                'kt',
                'kdeg',
                
                'fub_parent',
                'Rblood2plasma_parent',
                'Papp_parent',
                
                'EC50',       
                'Emax',   
                'fmCYP', 
                                                                # vary flow will make model unstable, so do not check morris of flow rate!!!
                'Vart' ,
                'Vgut',
                'Vliver',
                'Vlung' ,
                'Vbrain' ,
                'Vadipose',
                'Vrest',
                'Vven' , 
 
                'Kgut2pu_parent',           
                'Kliver2pu_parent',        
                'Klung2pu_parent' ,       
                'Kadipose2pu_parent',   
                'Kbrain2pu_parent',   
                'Krest2pu_parent' )

length(LVpars)
limit      <- c( rep(1.5, 11),  rep(1.3, 14))

baseline   <- c(Km_hep_parent     = 4.6,                  
                Vmax_hep_parent   = 0.045264,
                ka                = 1.123657,
                kt                = 0.054,
                kdeg              = 0.03,
                
                fub_parent            = 0.05,
                Rblood2plasma_parent  = 0.758,
                Papp_parent           = 40.8,
                
                EC50              = 3.57,       
                Emax              = 7.5,   
                fmCYP             = 0.65, 
                
                Vart       = Vart,
                Vgut       = Vgut,
                Vliver     = Vliver,  # 0.03486
                Vlung      = Vlung,
                Vbrain     = Vbrain,
                Vadipose   = Vadipose,
                Vrest      = Vrest,
                Vven       = Vven,
                
                Kgut2pu_parent       = 236.6 ,     
                Kliver2pu_parent     = 106.8,     
                Klung2pu_parent      = 180,   
                Kadipose2pu_parent   = 4141, 
                Kbrain2pu_parent     = 487.5,
                Krest2pu_parent      = 122.1511)

LVbinf     <- baseline/limit
LVbsup     <- baseline*limit


dat        <- list()
morris     <- list()
morris_p   <- list()
morris_i   <- list()
morris_i_p <- list()
morrisM     <- list()
morrisM_p   <- list()
morrisM_i   <- list()
morrisM_i_p <- list()

ii = 24

for (jj in duration) {
  
  ii = jj * 24 
  
  cat('ii = ', ii, '\n')
  morris_test <- ODEmorris(mod   = pbpk8cpt, 
                           pars  = LVpars,
                           state_init = initState,    ## Setting of the initial values of the state variables
                           times = ii,                ## time v https://rdrr.io/cran/ODEsensitivity/man/ODEmorris.ODEnetwork.html
                           binf  = LVbinf, 
                           bsup  = LVbsup, 
                           r = 500,                   #  number of trajectories (r) https://gsa-module.readthedocs.io/en/stable/implementation/morris_screening_method.html
                           design = list(type = "oat", 
                                         levels = 5, grid.jump = 1),
                           scale = TRUE,
                           ode_method = "lsoda",
                           parallel_eval = TRUE,
                           parallel_eval_ncores = 2) 
  

  morris_results2      <- as.data.frame(morris_test$AUC_Cblood_parent)
  morris_results2_df   <- rbind(data.frame(time1 = 'AUC_Cblood_parent'), morris_results2)

  names            <- rownames(morris_results2)
  dat              <- cbind(dat,  unlist(morris_results2_df))
  
  # https://rdrr.io/cran/ODEsensitivity/man/plot.ODEmorris.html
  
  title           <- paste0("Time = ", ii, "hr")
  

  ###########         Cmax_blood_parent          ####################################
  morris_resultsM_p <- as.data.frame(morris_test$AUC_Cblood_parent)         # select target state variable
  dfmM_p            <- cbind(category = rownames(morris_resultsM_p), morris_resultsM_p)
  names(dfmM_p)[2]  <- 'value'
  rownames(dfmM_p)  <- NULL
  muM_p             <- dfmM_p[grep("mu.star", dfmM_p$category), ]
  sigmaM_p          <- dfmM_p[grep("sigma",   dfmM_p$category), ]
  
  muM_df_p     <- separate(data = muM_p,    col = category, into = c("category", "parameter"), sep = "star")
  sigmaM_df_p  <- separate(data = sigmaM_p, col = category, into = c("category", "parameter"), sep = "sigma")
  
  muM_df_p$parameter    <- gsub("^.{0,1}", "", muM_df_p$parameter)
  sigmaM_df_p$parameter <- gsub("^.{0,1}", "", sigmaM_df_p$parameter)
  
  morrisM_df_p           <- merge(muM_df_p, sigmaM_df_p, by = 'parameter')
  colnames(morrisM_df_p) <- c('parameter', 'mu', 'mu_value', 'sigma', 'sigma_value')
  morrisM_df_p$sigma     <- 'sigma'
  morrisM_df_p$logmu_value    <- log(morrisM_df_p$mu_value)
  morrisM_df_p$logsigma_value <- log(morrisM_df_p$sigma_value)
  
  #############################################################
  
  #Save to file
  name             <- 'Propiconazole'
  timing           <- paste0(ii, "hr")
  
  morrisM_df_p$Group  <- name
  morrisM_df_p$Time   <- timing
  morrisM_p           <- rbind(morrisM_p, morrisM_df_p)
  morrisM_i_p         <- rbind(morrisM_i_p, morrisM_df_p)
  
  
  mu_value_max        <- max(morrisM_df_p$logmu_value)
  sigma_value_max     <- max(morrisM_df_p$logsigma_value)
  mu_value_3          <- mu_value_max/2
  sigma_value_3       <- sigma_value_max/2
  
  print(paste0("Hour = ", ii, " hr"))
  #print(paste0("sigma_value_max = ", sigma_value_max))

  
}

dat <- cbind(names,dat)
write.csv(morrisM_i_p,  paste("C:/XXX/Mammals/PPZ/SA/morris_AUC_blood_parent_rat_", oral_mg.kg,"_mg.kg.csv", sep = ''), row.names=FALSE)

#########################################################################################
#################     calculate global index for Cmax_blood_parent_m     ################  
#########################################################################################
morrisM_i_p$index      <- (morrisM_i_p$mu_value^2 + morrisM_i_p$sigma_value^2)^0.5
morrisM_i_p$logindex   <- log10(morrisM_i_p$index)
morrisM_i_p            <- morrisM_i_p [morrisM_i_p$mu_value != 0, ]
morrisM_i_p$Time       <- factor(morrisM_i_p$Time, levels=c('1hr', '6hr', '12hr', '24hr', '48hr', '72hr'))
morrisM_i_p            <- morrisM_i_p[!is.na(morrisM_i_p$parameter), ]

##########                 bar plot                    ########## 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")


p  <- ggplot(morrisM_i_p, aes(x = reorder(parameter, -index), y = ((index)) , fill= Time)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7, color="black",size=0.4) +
  scale_x_discrete(guide = guide_axis(angle = 60))+
  #scale_fill_manual(values=cbPalette) +                                # http://rstudio-pubs-static.s3.amazonaws.com/5312_98fc1aba2d5740dd849a5ab797cc2c8d.html
  scale_fill_brewer(palette="GnBu")+
  #scale_fill_hue(l=40) + 
  #geom_text(aes(label = round(Percent,0)), vjust=1.6, color="white",
  #position = position_dodge(0.65), size=5) + 
  theme(aspect.ratio = 1.5, 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 12))+
  theme_bw(base_size = 14)+
  labs(y = "Morris Index") + labs(x = "") +
  ggtitle(expression(paste('Sensitivity Analysis - AUC'[C[blood]], '')))

p

ggsave(paste("Morris_AUCCblood_parent_rat_", oral_mg.kg, "_mg.kg.tiff", sep= ''),scale = 1,
       plot = p,
       path = "C:/XXX/HTTK/Mammals/PPZ/SA/",
       width = 20, height = 15, units = "cm", dpi=300, compression = "lzw")

