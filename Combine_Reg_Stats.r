# this is for combining the OM trait regression stats from across multiple datasets and then making plots
# the plots are based on the R2 and slope defined 2D spaces
# currently combining ECA and S19S

eca = read.csv("ECA_Trait_Reg_Stats.csv",stringsAsFactors = F)
s19s = read.csv("S19S_Trait_Reg_Stats.csv",stringsAsFactors = F)

# a check, must be TRUE
identical(colnames(eca),colnames(s19s))

stats.comp = rbind(eca,s19s)
str(stats.comp)

#### plotting

pdf("Combined_R2_v_R2.pdf")

  par(pty="s")
  plot(stats.comp$Comp.mol.R2 ~ stats.comp$C.mol.R2,typ="n",ylim=c(0,1),xlim=c(0,1),cex=0.5,
       xlab=expression(R^{2} ~ "for" ~ lambda ~ vs. ~ Delta*G[cox] ~ (C-mol^{-1})),
       ylab=expression(R^{2} ~ "for" ~ lambda ~ vs. ~ Delta*G[cox] ~ (Comp-mol^{-1})),
       cex.lab=2,cex.axis=1.5)
  abline(v=0.5,lwd=2,col=8,lty=2)
  abline(h=0.5,lwd=2,col=8,lty=2)
  abline(0,1,lwd=3,col=8,lty=3)

  points(stats.comp$Comp.mol.R2[grep(pattern = "Sed",x = stats.comp$Sample_ID)] ~ stats.comp$C.mol.R2[grep(pattern = "Sed",x = stats.comp$Sample_ID)],cex=0.5,col=2)
  points(stats.comp$Comp.mol.R2[grep(pattern = "ECA",x = stats.comp$Sample_ID)] ~ stats.comp$C.mol.R2[grep(pattern = "ECA",x = stats.comp$Sample_ID)],cex=0.5,col=1)
  points(stats.comp$Comp.mol.R2[-grep(pattern = "Sed|ECA",x = stats.comp$Sample_ID)] ~ stats.comp$C.mol.R2[-grep(pattern = "Sed|ECA",x = stats.comp$Sample_ID)],cex=0.5,col=4)

  legend(x = 0.6,y = 0.4,legend = c("Inland Sediment","Inland Water","Parafluvial Sediment"),col=c(2,4,1),pch=1)
  
dev.off()

axis.max = ceiling(max(range(stats.comp$Comp.mol.slope),range(stats.comp$C.mol.slope)))

pdf("Combined_Slope_v_Slope.pdf")

par(pty="s")
plot(stats.comp$Comp.mol.slope ~ stats.comp$C.mol.slope,typ="n",ylim=c(0,axis.max),xlim=c(0,axis.max),cex=0.5,
     xlab=expression("Slope" ~ "for" ~ lambda ~ vs. ~ Delta*G[cox] ~ (C-mol^{-1})),
     ylab=expression("Slope" ~ "for" ~ lambda ~ vs. ~ Delta*G[cox] ~ (Comp-mol^{-1})),
     cex.lab=2,cex.axis=1.5)
  abline(0,1,lwd=3,col=8,lty=3)

  points(stats.comp$Comp.mol.slope[grep(pattern = "Sed",x = stats.comp$Sample_ID)] ~ stats.comp$C.mol.slope[grep(pattern = "Sed",x = stats.comp$Sample_ID)],cex=0.5,col=2)
  points(stats.comp$Comp.mol.slope[grep(pattern = "ECA",x = stats.comp$Sample_ID)] ~ stats.comp$C.mol.slope[grep(pattern = "ECA",x = stats.comp$Sample_ID)],cex=0.5,col=1)
  points(stats.comp$Comp.mol.slope[-grep(pattern = "Sed|ECA",x = stats.comp$Sample_ID)] ~ stats.comp$C.mol.slope[-grep(pattern = "Sed|ECA",x = stats.comp$Sample_ID)],cex=0.5,col=4)

  legend(x = 1,y = 0.5,legend = c("Inland Sediment","Inland Water","Parafluvial Sediment"),col=c(2,4,1),pch=1)

dev.off()

