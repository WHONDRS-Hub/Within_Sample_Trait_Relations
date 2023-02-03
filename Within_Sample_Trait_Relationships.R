# ==============================================================================
#
# Quantifying relationships among thermodynamic traits at a within-sample level
#
# Status: In progress
# 
# ==============================================================================
# To do: figure out directory (remove full.names = T) , update samples variable
# why are my slopes different than james, probably should change gibbs colnames early on to prevent confusion
#
# ==============================================================================
#
# Author: James Stegen, Brieanne Forbes
# 8 Dec 2022
#
# ==============================================================================

rm(list=ls())

library(tidyverse)
library(glue)
library(viridisLite)
library(ggpubr)
library(cowplot)

# ================================= User inputs ================================

study <- 'ECA'

# ================================ theme =======================================

theme_set(
  theme(
    text = element_text(family = 'serif'),
    axis.title = element_text(size = 16, face = 'bold'),
    axis.text = element_text(size = 14, face = 'bold'),
    line = element_line(size = 0.05),
    axis.line = element_line(size = 0.5),
    panel.background = element_rect(color = 'white'),
    panel.border = element_rect(
      colour = 'black',
      fill = NA,
      size = 0.5
    ),
    plot.title = element_text(size = 25, face = 'bold'),
    legend.position = 'none'
  )
)
# =============================== find and read data ===========================
mol_file <- list.files(pattern = paste0(study, '.+Mol'), full.names = T)

dat_file <- list.files(pattern = paste0(study, '.+Data'), full.names = T)

mol <- read_csv(mol_file)%>%
  rename(Mass=1) %>%
  select(Mass,'delGd','lamO2','delGcoxPerCmol')

dat <- read_csv(dat_file)%>%
  rename(Mass=1)


# =========================== remove peaks with extreme values =================

mol <- mol %>%
  filter(mol$lamO2 < 1, 
         mol$lamO2 > 0)

dat <- dat %>%
  filter(Mass %in% mol$Mass)

# hist(mol$delGd)
# hist(mol$lamO2); range(mol$lamO2,na.rm = T); length(which(mol$lamO2 < 0)); length(which(mol$lamO2 > 1))
# hist(mol$delGcoxPerCmol); range(mol$delGcoxPerCmol,na.rm = T)

# ======================== turning mol variables into z-scores =================

mol <- mol %>%
  mutate(delGd = (delGd-mean(delGd,na.rm = T))/sd(delGd,na.rm = T),
         lamO2 = (lamO2-mean(lamO2,na.rm = T))/sd(lamO2,na.rm = T),
         delGcoxPerCmol = (delGcoxPerCmol-mean(delGcoxPerCmol,na.rm = T))/sd(delGcoxPerCmol,na.rm = T))

# the printing should be 0 1 for each variable
print(c(round(mean(mol$delGd,na.rm = T),digits = 10),sd(mol$delGd,na.rm = T)))
print(c(round(mean(mol$lamO2,na.rm = T),digits = 10),sd(mol$lamO2,na.rm = T)))
print(c(round(mean(mol$delGcoxPerCmol,na.rm = T),digits = 10),sd(mol$delGcoxPerCmol,na.rm = T)))

# ==== get R2 and slope of each peak within sample (Gibbs vs lambda) ===========

sample_stats_combine = tibble(Sample_Name = as.character(),
                      Number_of_Formulas = as.numeric(),
                      Comp_mol_lambda_R2 = as.numeric(), 
                      C_mol_lambda_R2 = as.numeric(), 
                      Comp_mol_lambda_slope = as.numeric(), 
                      C_mol_lambda_slope = as.numeric(),
                      Comp_mol_lambda_eq = as.character(), 
                      C_mol_lambda_eq = as.character(),
                      Comp_mol_C_mol_Rsq_Ratio = as.numeric(),
                      Comp_mol_C_mol_Slope_Ratio = as.numeric() 
                      )

all_sample_peaks <- tibble('Sample_Name' = as.character(),
                       "delGd" = as.numeric(),
                       "lamO2" = as.numeric(),
                       "delGcoxPerCmol" = as.numeric()
                       )

all_samples <- colnames(dat) %>%
  str_remove('Mass')

random_samples <- sample(all_samples, 50)



for (i in random_samples) {
  
  dat_temp <- dat %>%
    select('Mass', i)%>%
    filter(get(i) == 1)
    
  mol_temp <- mol %>%
    filter(!is.na(lamO2))
  
  combine <- dat_temp %>%
    left_join(mol_temp)%>%
    filter(!is.na(delGd) & !is.na(lamO2) & !is.na(delGcoxPerCmol))
  
  combine_filter <- combine %>%
    select("delGd","lamO2","delGcoxPerCmol") %>%
    add_column(Sample_Name = i, .before = 'delGd')
  
  all_sample_peaks <- all_sample_peaks %>%
    add_row(combine_filter)
  
  Comp_mol_lm <- summary(lm(combine$lamO2 ~ combine$delGd))
  Comp_mol_R2 <- Comp_mol_lm$r.squared # R2 for Gibbs per Comp mol (pH 7)
  Comp_mol_slope <- Comp_mol_lm$coefficients[2,1] # slope for Gibbs per Comp mol (pH 7)
  Comp_mol_eq <- glue('y = {round(Comp_mol_lm$coefficients[2],3)}x + {round(Comp_mol_lm$coefficients[1], 3)}')
  
  C_mol_lm <- summary(lm(combine$lamO2 ~ combine$delGcoxPerCmol)) 
  C_mol_R2 <- C_mol_lm$r.squared # R2 for Gibbs per C mol (pH 7)
  C_mol_slope <- C_mol_lm$coefficients[2,1] # slope for Gibbs per C mol (pH 7)
  C_mol_eq <- glue('y = {round(C_mol_lm$coefficients[2],3)}x + {round(C_mol_lm$coefficients[1], 3)}')
  
    sample_stats <- tibble(
      Sample_Name = i,
      Number_of_Formulas = nrow(combine),
      Comp_mol_lambda_R2 = Comp_mol_R2,
      C_mol_lambda_R2 = C_mol_R2,
      Comp_mol_lambda_slope = Comp_mol_slope,
      C_mol_lambda_slope = C_mol_slope,
      Comp_mol_lambda_eq = Comp_mol_eq,
      C_mol_lambda_eq = C_mol_eq,
      Comp_mol_C_mol_Rsq_Ratio = (Comp_mol_R2/C_mol_R2),
      Comp_mol_C_mol_Slope_Ratio = (Comp_mol_slope/C_mol_slope)
    )
  
  sample_stats_combine <- sample_stats_combine %>%
    add_row(sample_stats)%>%
    arrange(Sample_Name)
  # ================================  plot =======================================
  
  # colors <- c("#2A788EFF", "#7AD151FF")
  
  c_mol <- ggplot(data = combine_filter, aes(x = delGcoxPerCmol, y = lamO2, color = "#2A788EFF")) +
    geom_point(size = .75)+
    geom_smooth(method=lm, size = .75, color = 'black')+
    labs(x = 'Gibbs (C-mol) ', y = 'Lambda')+ 
    scale_color_manual(values = "#2A788EFF")
  
  c_mol <- ggdraw(c_mol) + 
    draw_label(bquote(R^2 == . (round(sample_stats$C_mol_lambda_R2, 2))), x = 0.3, y = 0.9, fontfamily = 'serif')+
    draw_label(sample_stats$C_mol_lambda_eq, x = 0.35, y = 0.85, fontfamily = 'serif')
  c_mol
  
  comp_mol<- ggplot(data = combine_filter, aes(x = delGd, y = lamO2, color = "#2A788EFF")) +
    geom_point(size = .75)+
    geom_smooth(method=lm, size = .75, color = 'black')+
    labs(x = 'Gibbs (Comp-mol) ', y = 'Lambda')+
    # annotate('text', x = quantile(combine_filter$delGd, probs = 0.95, na.rm = T), y = quantile(combine_filter$lamO2, probs = 0.95, na.rm = T), label = paste('R^2 ==',round(sample_stats$Comp_mol_lambda_R2, 2)), parse = T, size = 4, family = 'serif')+
    # annotate('text', x = quantile(combine_filter$delGd, probs = 0.95, na.rm = T), y = quantile(combine_filter$lamO2, probs = 0.95, na.rm = T), label = sample_stats$Comp_mol_lambda_eq, size = 4, family = 'serif')+
    scale_color_manual(values = "#2A788EFF")
  
  comp_mol <- ggdraw(comp_mol) + 
    draw_label(bquote(R^2 == . (round(sample_stats$Comp_mol_lambda_R2, 2))), x = 0.3, y = 0.9, fontfamily = 'serif')+
    draw_label(sample_stats$Comp_mol_lambda_eq, x = 0.35, y = 0.85, fontfamily = 'serif')
  comp_mol
  
  combine_plots <- ggarrange(
    c_mol,
    comp_mol,
    ncol = 2,
    nrow = 1,
    widths = c(4),
    heights = c(4, 4)
  )
  
  combine_plots <-
    annotate_figure(combine_plots, top = text_grob(
      i,
      size = 14,
      family = 'serif',
      face = 'bold'
    ))
  
  combine_plots_out <- paste0(here::here(),'/ECA_Within_Sample_Peaks_Plots/', i, '.pdf' )
  
  ggsave(
    combine_plots_out,
    combine_plots,
    device = 'pdf',
    width = 10,
    height = 5,
    units = 'in',
    dpi = 300
  )
}






# ECA2_0033_03.a_p08 <-  ECA2_0033_03.a_p08 %>%
#   select("delGd","lamO2","delGcoxPerCmol") %>%
#   add_column(Sample_Name = 'ECA2_0033_03.a_p08', .before = 'delGd')
# 
# ECA2_0054_09.a_p15 <- ECA2_0054_09.a_p15 %>%
#   add_column(Sample_Name = 'ECA2_0054_09.a_p15')


# sample.reg.stats = as.data.frame(sample.reg.stats)
# colnames(sample.reg.stats) = c('Sample_ID','Num_of_formulas','Comp.mol.R2','Comp.mol.slope','C.mol.R2','C.mol.slope')
# sample.reg.stats[1:ncol(sample.reg.stats)] = lapply(sample.reg.stats[1:ncol(sample.reg.stats)],as.character)
# sample.reg.stats[which(colnames(sample.reg.stats) != 'Sample_ID')] = lapply(sample.reg.stats[which(colnames(sample.reg.stats) != 'Sample_ID')],as.numeric)
# 
# #### read in file for bad calibrations and remove those samples
# bad.cal = read.csv("S19S_Water_Sed_Hawkes_Poorly_Calibrated_Samples.csv",stringsAsFactors = F)
# # changing hyphens to periods because that happened during the processing of the data, so this is needed to match sample names
# bad.cal$samples = gsub(pattern = "-",replacement = ".",x = bad.cal$samples)
# # this is a check and should be zero
# length(which(bad.cal$samples %in% sample.reg.stats$Sample_ID)) - nrow(bad.cal)
# 
# # this is for a check
# orig.nrow = nrow(sample.reg.stats)
# 
# # now remove bal cal samples
# sample.reg.stats = sample.reg.stats[-which(sample.reg.stats$Sample_ID %in% bad.cal$samples),]
# 
# # this is a check and should be zero
# orig.nrow - nrow(sample.reg.stats) - nrow(bad.cal)
# 
# #### make histograms to look at variation in R2 and slopes
# 
# hist(sample.reg.stats$Num_of_formulas)
# range(sample.reg.stats$Num_of_formulas)
# 
# hist(sample.reg.stats$Comp.mol.R2)
# hist(sample.reg.stats$C.mol.R2)
# 
# pdf("S19S_R2_v_R2.pdf")
# 
#   par(pty="s")
#   plot(sample.reg.stats$Comp.mol.R2 ~ sample.reg.stats$C.mol.R2,typ="n",ylim=c(0,1),xlim=c(0,1),cex=0.5,xlab="R2: Lambda v. Gibbs (C-mol)",ylab="R2: Lambda v. Gibbs (Comp-mol)",cex.lab=2,cex.axis=1.5)
#   abline(v=0.5,lwd=2,col=8,lty=2)
#   abline(h=0.5,lwd=2,col=8,lty=2)
#   abline(0,1,lwd=3,col=3,lty=3)
#   
#   points(sample.reg.stats$Comp.mol.R2[grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)] ~ sample.reg.stats$C.mol.R2[grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)],cex=0.5,col=2)
#   points(sample.reg.stats$Comp.mol.R2[-grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)] ~ sample.reg.stats$C.mol.R2[-grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)],cex=0.5,col=4)
#   
#   
# dev.off()
# 
# hist(sample.reg.stats$Comp.mol.slope)
# hist(sample.reg.stats$C.mol.slope)
# pdf("S19S_Slope_v_Slope.pdf")
# 
#   par(pty="s")
#   plot(sample.reg.stats$Comp.mol.slope ~ sample.reg.stats$C.mol.slope,xlim=c(0,2),ylim=c(0,2),typ="n",cex=0.5,xlab="Slope: Lambda v. Gibbs (C-mol)",ylab="Slope: Lambda v. Gibbs (Comp-mol)",cex.lab=2,cex.axis=1.5)
#   #abline(v=0.75,lwd=2,col=8,lty=2)
#   #abline(h=0.75,lwd=2,col=8,lty=2)
#   abline(0,1,lwd=3,col=3,lty=3)
#   
#   points(sample.reg.stats$Comp.mol.slope[grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)] ~ sample.reg.stats$C.mol.slope[grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)],cex=0.5,col=2)
#   points(sample.reg.stats$Comp.mol.slope[-grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)] ~ sample.reg.stats$C.mol.slope[-grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)],cex=0.5,col=4)
#   
# dev.off()
# 
# sample.reg.stats$Comp.R2_C.R2_Ratio = sample.reg.stats$Comp.mol.R2/sample.reg.stats$C.mol.R2
# sample.reg.stats$Comp.Slope_C.Slope_Ratio = sample.reg.stats$Comp.mol.slope/sample.reg.stats$C.mol.slope
# 
# hist(sample.reg.stats$Comp.R2_C.R2_Ratio)
# hist(sample.reg.stats$Comp.Slope_C.Slope_Ratio)
# plot(sample.reg.stats$Comp.R2_C.R2_Ratio ~ sample.reg.stats$Comp.Slope_C.Slope_Ratio)
# 
# write.csv(sample.reg.stats,"S19S_Trait_Reg_Stats.csv",row.names = F,quote = F)

