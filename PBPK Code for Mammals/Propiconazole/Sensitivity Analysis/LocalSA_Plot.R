library(readxl)

rm(list = ls())


SA_mouse      <- read.csv('C:/XXX/PPZ/SA/LocalSA_Mouse.csv')
SA_rat        <- read.csv('C:/XXX/PPZ/SA/LocalSA_rat.csv')
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

##########                 bar plot                    ########## 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

p  <- ggplot(SA_AUC, aes(x = reorder(parameter, -AUC), y = ((AUC)) , fill= Species)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7, color="black",size=0.4) +
  scale_x_discrete(guide = guide_axis(angle = 60))+
  scale_fill_manual(values=cbPalette) +                                # http://rstudio-pubs-static.s3.amazonaws.com/5312_98fc1aba2d5740dd849a5ab797cc2c8d.html
  #scale_fill_brewer(palette="GnBu")+
  #scale_fill_hue(l=40) + 
  #geom_text(aes(label = round(Percent,0)), vjust=1.6, color="white",
  #position = position_dodge(0.65), size=5) + 
  theme(aspect.ratio = 1.5, 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 12))+
  theme_bw(base_size = 16)+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14),
        plot.margin = margin(l = 10, unit = "mm"))+
  labs(y = "Absolute NSC Value [%]") + labs(x = "") +
  #ggtitle(expression(paste('Sensitivity Analysis - AUC'[C[blood]], '')))
  ggtitle(expression(paste('Sensitivity Analysis - AUC')))



library(tidyverse)

# Ensure all combinations of 'parameter' and 'Species' are present
SA_AUC <- SA_AUC %>%
  complete(parameter, Species, fill = list(AUC = 0))

# Then plot
p  <- ggplot(SA_AUC, aes(x = reorder(parameter, -abs(AUC)), y = ((AUC)) , fill= Species)) +
  geom_hline(yintercept = 0.5,  linetype="dashed", color = "#FC4E07", alpha = 0.7, size = 0.4)  +           # Add horizontal line at y=0.5
  geom_hline(yintercept = -0.5, linetype="dashed", color = "#FC4E07", alpha = 0.7,size = 0.4)  +           # Add horizontal line at y=-0.5
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.2, ymax = 0.2, alpha = 0.6, fill = "grey") + #
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.5, ymax = 0.5, alpha = 0.2, fill = "grey") + #
  geom_bar(stat="identity", position=position_dodge(width = 0.7), width = 0.7, color="black",size=0.4) +
  scale_x_discrete(guide = guide_axis(angle = 60))+
  scale_y_continuous(breaks = seq(-4, 2, by = 1))+
  scale_fill_manual(values=cbPalette) +
  theme(aspect.ratio = 1.5, 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 12))+
  theme_bw(base_size = 16)+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14),
        plot.margin = margin(l = 10, unit = "mm"))+
  labs(y = "Sensitivity Ratio") + labs(x = "") +
  ggtitle(expression(paste('Sensitivity Analysis - AUC')))

p

ggsave(paste("LocalSA_AUC.tiff", sep= ''),scale = 1,
       plot = p,
       path = "C:/XXX/Mammals/PPZ/SA/",
       width = 24, height = 15, units = "cm", dpi=600, compression = "lzw")


cbPalette <- c("#009E73", "#F0E442", "#0072B2")

p  <- ggplot(SA_Cmax, aes(x = reorder(parameter, -abs(Cmax)), y = ((Cmax)) , fill= Species)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7, color="black",size=0.4) +
  scale_x_discrete(guide = guide_axis(angle = 60))+
  scale_fill_manual(values=cbPalette) +                                # http://rstudio-pubs-static.s3.amazonaws.com/5312_98fc1aba2d5740dd849a5ab797cc2c8d.html
  #scale_fill_brewer(palette="GnBu")+
  #scale_fill_hue(l=40) + 
  #geom_text(aes(label = round(Percent,0)), vjust=1.6, color="white",
  #position = position_dodge(0.65), size=5) + 
  theme(aspect.ratio = 1.5, 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 12))+
  theme_bw(base_size = 16)+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14),
        plot.margin = margin(l = 10, unit = "mm"))+
  labs(y = "Absolute NSC Value [%]") + labs(x = "") +
  ggtitle(expression(paste('Sensitivity Analysis - C'[max], '')))

p

ggsave(paste("LocalSA_Cmax.tiff", sep= ''),scale = 1,
       plot = p,
       path = "C:/XXX/Mammals/PPZ/SA/",
       width = 24, height = 15, units = "cm", dpi=600, compression = "lzw")
