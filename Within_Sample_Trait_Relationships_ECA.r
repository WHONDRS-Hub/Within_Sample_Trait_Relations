# this is for quantifying relationships among thermodynamic traits at a within-sample level
# the questions:
# 1) how variable are the trait relationships in terms of slope and R2?
#     a) if variation is low, this indicates universal constraints that are unique to sediment systems (note surface water show a weak relation)
#     b) if variation is high, this indicates the trait relationships are influenced by local conditions as well as some global constraints (note the sediment relationships were consistently strong in chemogeography)  
# 2) Is variation in the trait relationships larger between sites that within sites?
#     a) if yes, this indicates something at the reach scale constrains these relationships
#     b) if no, this indicates that there are very small scale, sample-specific features/processes that impact the relationships
# 3) Does sediment moisture content relate to the slope or R2 of the trait relationships?
#     a) if yes, this indicates that gas-phase oxygen supply or redox more generally could be a key factor
#     b) if no, this indicates that oxygen/redox may not be a key influence over the relatioships
# 4) Is one relationship consistently stronger/steeper than the other?
#     a) if yes, this points to a real/true/valid/important bias in how energy content of molecules influences how microbes use those molecules to build biomass
#     b) if no, then energy content per C and per compound may be equally important for how microbes use orgnaic molecules to build biomass

rm(list=ls())

#### read in the mol and data files (these are early career data from 2020)
mol = read.csv(file = "Processed_ECA2_2_out_of_14_Clean_Mol.csv",stringsAsFactors = F,row.names = 1)
mol = mol[,c('delGd','lamO2','delGcoxPerCmol')]
dat = read.csv(file = "Processed_ECA2_2_out_of_14_Clean_Data.csv",stringsAsFactors = F,row.names = 1)

# there are a few peaks that have extreme values for lambda. not sure why, but removing them
peaks.to.remove = c(
  
  rownames(mol)[which(mol$lamO2 < 0)],
  rownames(mol)[which(mol$lamO2 > 1)]
  
)

mol = mol[-which(rownames(mol) %in% peaks.to.remove),]
dat = dat[-which(rownames(dat) %in% peaks.to.remove),]

hist(mol$delGd)
hist(mol$lamO2); range(mol$lamO2,na.rm = T); length(which(mol$lamO2 < 0)); length(which(mol$lamO2 > 1))
hist(mol$delGcoxPerCmol); range(mol$delGcoxPerCmol,na.rm = T)

# turning mol variables into z-scores. the printing should be 0 1 for each variable
mol$delGd = (mol$delGd-mean(mol$delGd,na.rm = T))/sd(mol$delGd,na.rm = T); print(c(round(mean(mol$delGd,na.rm = T),digits = 10),sd(mol$delGd,na.rm = T)))
mol$lamO2 = (mol$lamO2-mean(mol$lamO2,na.rm = T))/sd(mol$lamO2,na.rm = T); print(c(round(mean(mol$lamO2,na.rm = T),digits = 10),sd(mol$lamO2,na.rm = T)))
mol$delGcoxPerCmol = (mol$delGcoxPerCmol-mean(mol$delGcoxPerCmol,na.rm = T))/sd(mol$delGcoxPerCmol,na.rm = T); print(c(round(mean(mol$delGcoxPerCmol,na.rm = T),digits = 10),sd(mol$delGcoxPerCmol,na.rm = T)))

#### loop through each sample and quantify the R2 and slope of the linear regression each Gibbs vs. lambda
#### compile outcomes into a dataframe

sample.reg.stats = numeric()

for (i in 1:ncol(dat)) {
  
  mol.temp = mol[which(rownames(mol) %in% rownames(dat)[which(dat[,i] > 0)]),]
  mol.temp = mol.temp[-which(is.na(mol.temp$lamO2)==T),]
  
  if (length(which(is.na(mol.temp) == T)) > 0) {
  
    print("Error: Not all NAs removed")  
    break()
    
  }
  print(c(colnames(dat)[i],range(mol.temp),"last 2 must be numbers"))
  
  Comp.mol.R2 = summary(lm(mol.temp$lamO2 ~ mol.temp$delGd))$r.squared # R2 for Gibbs per Comp mol (pH 7)
  Comp.mol.slope = summary(lm(mol.temp$lamO2 ~ mol.temp$delGd))$coefficients[2,1] # slope for Gibbs per Comp mol (pH 7)
  
  C.mol.R2 = summary(lm(mol.temp$lamO2 ~ mol.temp$delGcoxPerCmol))$r.squared # R2 for Gibbs per C mol (pH 7)
  C.mol.slope = summary(lm(mol.temp$lamO2 ~ mol.temp$delGcoxPerCmol))$coefficients[2,1] # slope for Gibbs per C mol (pH 7)
  
  sample.reg.stats = rbind(sample.reg.stats,c(
    
    colnames(dat)[i],
    nrow(mol.temp),
    Comp.mol.R2,
    Comp.mol.slope,
    C.mol.R2,
    C.mol.slope
    
  ))
  
}

sample.reg.stats = as.data.frame(sample.reg.stats)
colnames(sample.reg.stats) = c('Sample_ID','Num_of_formulas','Comp.mol.R2','Comp.mol.slope','C.mol.R2','C.mol.slope')
sample.reg.stats[1:ncol(sample.reg.stats)] = lapply(sample.reg.stats[1:ncol(sample.reg.stats)],as.character)
sample.reg.stats[which(colnames(sample.reg.stats) != 'Sample_ID')] = lapply(sample.reg.stats[which(colnames(sample.reg.stats) != 'Sample_ID')],as.numeric)
sample.reg.stats$Sample_ID = substring(sample.reg.stats$Sample_ID,first = 1,last = 12)

#### make histograms to look at variation in R2 and slopes

hist(sample.reg.stats$Num_of_formulas)
range(sample.reg.stats$Num_of_formulas)

hist(sample.reg.stats$Comp.mol.R2)
hist(sample.reg.stats$C.mol.R2)
pdf("R2_v_R2.pdf")
  par(pty="s")
  plot(sample.reg.stats$Comp.mol.R2 ~ sample.reg.stats$C.mol.R2,ylim=c(0,1),xlim=c(0,1),cex=0.5,xlab="R2: Lambda v. Gibbs (C-mol)",ylab="R2: Lambda v. Gibbs (Comp-mol)",cex.lab=2,cex.axis=1.5)
  abline(v=0.5,lwd=2,col=8,lty=2)
  abline(h=0.5,lwd=2,col=8,lty=2)
  abline(0,1,lwd=3,col=4,lty=3)
dev.off()

hist(sample.reg.stats$Comp.mol.slope)
hist(sample.reg.stats$C.mol.slope)
pdf("Slope_v_Slope.pdf")
  par(pty="s")
  plot(sample.reg.stats$Comp.mol.slope ~ sample.reg.stats$C.mol.slope,cex=0.5,xlim=c(0,1.5),ylim=c(0,1.5),xlab="Slope: Lambda v. Gibbs (C-mol)",ylab="Slope: Lambda v. Gibbs (Comp-mol)",cex.lab=2,cex.axis=1.5)
  abline(v=0.75,lwd=2,col=8,lty=2)
  abline(h=0.75,lwd=2,col=8,lty=2)
  abline(0,1,lwd=3,col=4,lty=3)
dev.off()

sample.reg.stats$Comp.R2_C.R2_Ratio = sample.reg.stats$Comp.mol.R2/sample.reg.stats$C.mol.R2
sample.reg.stats$Comp.Slope_C.Slope_Ratio = sample.reg.stats$Comp.mol.slope/sample.reg.stats$C.mol.slope

hist(sample.reg.stats$Comp.R2_C.R2_Ratio)
hist(sample.reg.stats$Comp.Slope_C.Slope_Ratio)
plot(sample.reg.stats$Comp.R2_C.R2_Ratio ~ sample.reg.stats$Comp.Slope_C.Slope_Ratio)

write.csv(sample.reg.stats,"ECA_Trait_Reg_Stats.csv",row.names = F,quote = F)

#### bring in sediment moisture

moisture = read.csv("merged_weights.csv")
moisture$Site = substring(moisture$ID_scanned,first = 1,last = 9)
moisture$ID_scanned = substring(moisture$ID_scanned,first = 1,last = 12)
head(moisture)

#### merge sediment with trait relationship data

sample.reg.stats = merge(x = sample.reg.stats,y = moisture,by.x = 'Sample_ID',by.y = 'ID_scanned',all.x = T)
dim(sample.reg.stats)

#### regressions against moisture

plot(sample.reg.stats$Comp.R2_C.R2_Ratio ~ sample.reg.stats$percent_water_content_dry)
plot(sample.reg.stats$Comp.R2_C.R2_Ratio ~ sample.reg.stats$percent_water_content_wet)

plot(sample.reg.stats$Comp.mol.R2 ~ sample.reg.stats$percent_water_content_dry)
plot(sample.reg.stats$C.mol.R2 ~ sample.reg.stats$percent_water_content_dry)

pdf("R2_v_Wet_Moisture.pdf",height=14)

  par(pty="s",mfrow=c(2,1))
  mod = sample.reg.stats$Comp.mol.R2 ~ sample.reg.stats$percent_water_content_wet
  plot(mod,cex=0.5,xlab="% Moisture per Wet Mass",ylab="R2: Lambda v. Gibbs (Comp-mol)",cex.lab=2,cex.axis=1.5)
  # mod.sum = summary(lm(mod))
  #abline(mod.sum,lwd=2,col=2)
  #mod.sum
  # Multiple R-squared:  0.01851,	Adjusted R-squared:  0.01632 
  # F-statistic: 8.465 on 1 and 449 DF,  p-value: 0.003799

  mod = sample.reg.stats$C.mol.R2 ~ sample.reg.stats$percent_water_content_wet
  plot(mod,cex=0.5,xlab="% Moisture per Wet Mass",ylab="R2: Lambda v. Gibbs (C-mol)",cex.lab=2,cex.axis=1.5)
  # mod.sum = summary(lm(mod))
  #abline(mod.sum,lwd=2,col=2)
  #mod.sum
  # Multiple R-squared:  0.003267,	Adjusted R-squared:  0.001048 
  # F-statistic: 1.472 on 1 and 449 DF,  p-value: 0.2257

dev.off()

#### quantify variance within each site

within.site.variance.R2 = numeric()

unique.sites = unique(sample.reg.stats$Site)

for (i in unique.sites) {
  
  if (length(sample.reg.stats$Comp.mol.R2[which(sample.reg.stats$Site == i)]) >= 8) {
  
    Comp.R2.var = 1000*var(sample.reg.stats$Comp.mol.R2[which(sample.reg.stats$Site == i)])
    C.R2.var = 1000*var(sample.reg.stats$C.mol.R2[which(sample.reg.stats$Site == i)])
    within.site.variance.R2 = rbind(within.site.variance.R2,c(
      
      i,
      Comp.R2.var,
      C.R2.var,
      length(sample.reg.stats$Comp.mol.R2[which(sample.reg.stats$Site == i)])
      
    ))
  
  } else {print(i)}
  
}

within.site.variance.R2 = as.data.frame(within.site.variance.R2)
colnames(within.site.variance.R2) = c('Site','Comp.R2.var','C.R2.var','Num.of.Samples')
within.site.variance.R2[1:ncol(within.site.variance.R2)] = lapply(within.site.variance.R2[1:ncol(within.site.variance.R2)],as.character)
within.site.variance.R2[which(colnames(within.site.variance.R2) != 'Site')] = lapply(within.site.variance.R2[which(colnames(within.site.variance.R2) != 'Site')],as.numeric)

#### calc between site variance by randomly sampling 10 samples from any site

between.site.variance.R2 = numeric()

for (i in 1:nrow(within.site.variance.R2)) {

rand.samples = numeric()

sampled.sites = sample(x = as.vector(na.omit(unique.sites)),size = 10,replace = F)

for (curr.site in sampled.sites) {
  
  rand.samples = c(rand.samples,sample(x = sample.reg.stats$Sample_ID[which(sample.reg.stats$Site == curr.site)],size = 1))

}

Comp.R2.var = 1000*var(sample.reg.stats$Comp.mol.R2[which( sample.reg.stats$Sample_ID %in% rand.samples )])
C.R2.var = 1000*var(sample.reg.stats$C.mol.R2[which( sample.reg.stats$Sample_ID %in% rand.samples )])

between.site.variance.R2 = rbind(between.site.variance.R2,c(Comp.R2.var,C.R2.var))

}

between.site.variance.R2 = as.data.frame(between.site.variance.R2)
colnames(between.site.variance.R2) = c('Comp.R2.var','C.R2.var')

#### overlay density plots for variances

btw.Comp.R2.var.density = density(between.site.variance.R2$Comp.R2.var,from = 0)
btw.C.R2.var.density = density(between.site.variance.R2$C.R2.var,from = 0)

wth.Comp.R2.var.density = density(within.site.variance.R2$Comp.R2.var,from = 0)
wth.C.R2.var.density = density(within.site.variance.R2$C.R2.var,from = 0)

pdf("Wth_Btw_R2_Var.pdf",height = 14)

  par(pty="s",mfrow=c(2,1))
  plot(btw.Comp.R2.var.density,typ="l",lwd=2,ylim=c(0,0.5),xlab="R2: Lambda v Gibbs Variance (Comp-mol)",cex.lab=2,cex.axis=1.5,main="")
  points(wth.Comp.R2.var.density,typ="l",col=2,lwd=2)
  legend(x = 2.5,y = 0.5,cex=1.5,legend = c("Within Site","Between Site"),col = c(2,1),lwd=2)
  wilcox.test(x = between.site.variance.R2$Comp.R2.var,y = within.site.variance.R2$Comp.R2.var)
  # W = 1416, p-value = 1.255e-05

  plot(btw.C.R2.var.density,typ="l",lwd=2,ylim=c(0,0.2),xlab="R2: Lambda v Gibbs Variance (C-mol)",cex.lab=2,cex.axis=1.5,main="")
  points(wth.C.R2.var.density,typ="l",col=2,lwd=2)
  #legend(x = 3.3,y = 0.5,legend = c("Within Site","Between Site"),col = c(2,1),lwd=2)
  wilcox.test(x = between.site.variance.R2$C.R2.var,y = within.site.variance.R2$C.R2.var)
  # W = 1514, p-value = 9.481e-08

dev.off()
