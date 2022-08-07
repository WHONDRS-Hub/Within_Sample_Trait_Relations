# this is for quantifying relationships among thermodynamic traits at a within-sample level, for S19S water and sediments
# the questions:
# 1) how variable are the trait relationships in terms of slope and R2?
#     a) if variation is low, this indicates universal constraints that are unique to sediment systems (note surface water show a weak relation)
#     b) if variation is high, this indicates the trait relationships are influenced by local conditions as well as some global constraints (note the sediment relationships were consistently strong in chemogeography)  
# 2) Is one relationship consistently stronger/steeper than the other?
#     a) if yes, this points to a real/true/valid/important bias in how energy content of molecules influences how microbes use those molecules to build biomass
#     b) if no, then energy content per C and per compound may be equally important for how microbes use orgnaic molecules to build biomass

rm(list=ls())

#### read in the mol and data files (these are S19S data)
mol = read.csv(file = "Processed_S19S_Water_Sed_Hawkes_Mol.csv",stringsAsFactors = F,row.names = 1)
mol = mol[,c('delGd','lamO2','delGcoxPerCmol')]
# README: this file is too large to store on github, so need to just have it stored in local version of the repo
dat = read.csv(file = "Processed_S19S_Water_Sed_Hawkes_Data.csv",stringsAsFactors = F,row.names = 1)

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

#### read in file for bad calibrations and remove those samples
bad.cal = read.csv("S19S_Water_Sed_Hawkes_Poorly_Calibrated_Samples.csv",stringsAsFactors = F)
# changing hyphens to periods because that happened during the processing of the data, so this is needed to match sample names
bad.cal$samples = gsub(pattern = "-",replacement = ".",x = bad.cal$samples)
# this is a check and should be zero
length(which(bad.cal$samples %in% sample.reg.stats$Sample_ID)) - nrow(bad.cal)

# this is for a check
orig.nrow = nrow(sample.reg.stats)

# now remove bal cal samples
sample.reg.stats = sample.reg.stats[-which(sample.reg.stats$Sample_ID %in% bad.cal$samples),]

# this is a check and should be zero
orig.nrow - nrow(sample.reg.stats) - nrow(bad.cal)

#### make histograms to look at variation in R2 and slopes

hist(sample.reg.stats$Num_of_formulas)
range(sample.reg.stats$Num_of_formulas)

hist(sample.reg.stats$Comp.mol.R2)
hist(sample.reg.stats$C.mol.R2)

pdf("S19S_R2_v_R2.pdf")

  par(pty="s")
  plot(sample.reg.stats$Comp.mol.R2 ~ sample.reg.stats$C.mol.R2,typ="n",ylim=c(0,1),xlim=c(0,1),cex=0.5,xlab="R2: Lambda v. Gibbs (C-mol)",ylab="R2: Lambda v. Gibbs (Comp-mol)",cex.lab=2,cex.axis=1.5)
  abline(v=0.5,lwd=2,col=8,lty=2)
  abline(h=0.5,lwd=2,col=8,lty=2)
  abline(0,1,lwd=3,col=3,lty=3)
  
  points(sample.reg.stats$Comp.mol.R2[grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)] ~ sample.reg.stats$C.mol.R2[grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)],cex=0.5,col=2)
  points(sample.reg.stats$Comp.mol.R2[-grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)] ~ sample.reg.stats$C.mol.R2[-grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)],cex=0.5,col=4)
  
  
dev.off()

hist(sample.reg.stats$Comp.mol.slope)
hist(sample.reg.stats$C.mol.slope)
pdf("S19S_Slope_v_Slope.pdf")

  par(pty="s")
  plot(sample.reg.stats$Comp.mol.slope ~ sample.reg.stats$C.mol.slope,xlim=c(0,2),ylim=c(0,2),typ="n",cex=0.5,xlab="Slope: Lambda v. Gibbs (C-mol)",ylab="Slope: Lambda v. Gibbs (Comp-mol)",cex.lab=2,cex.axis=1.5)
  #abline(v=0.75,lwd=2,col=8,lty=2)
  #abline(h=0.75,lwd=2,col=8,lty=2)
  abline(0,1,lwd=3,col=3,lty=3)
  
  points(sample.reg.stats$Comp.mol.slope[grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)] ~ sample.reg.stats$C.mol.slope[grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)],cex=0.5,col=2)
  points(sample.reg.stats$Comp.mol.slope[-grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)] ~ sample.reg.stats$C.mol.slope[-grep(pattern = "Sed",x = sample.reg.stats$Sample_ID)],cex=0.5,col=4)
  
dev.off()

sample.reg.stats$Comp.R2_C.R2_Ratio = sample.reg.stats$Comp.mol.R2/sample.reg.stats$C.mol.R2
sample.reg.stats$Comp.Slope_C.Slope_Ratio = sample.reg.stats$Comp.mol.slope/sample.reg.stats$C.mol.slope

hist(sample.reg.stats$Comp.R2_C.R2_Ratio)
hist(sample.reg.stats$Comp.Slope_C.Slope_Ratio)
plot(sample.reg.stats$Comp.R2_C.R2_Ratio ~ sample.reg.stats$Comp.Slope_C.Slope_Ratio)

write.csv(sample.reg.stats,"S19S_Trait_Reg_Stats.csv",row.names = F,quote = F)

