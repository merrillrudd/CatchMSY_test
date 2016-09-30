rm(list=ls())

library(catchMSY)
library(LIME)

# main_dir <- "F:\\Merrill\\Git_Projects\\CatchMSY_test"
main_dir <- "C:\\Git_Projects\\CatchMSY_test"
source(file.path(main_dir, "R", "test_functions.R"))

## setup results directory
sim_dir <- file.path(main_dir, "cross_test")
## remove results (for testing)
# unlink(sim_dir, TRUE)
## create directory
dir.create(sim_dir, showWarnings=FALSE)

## setup directory for figures
fig_dir <- file.path(main_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

##### --------------- testing options ------------------ #####

## test across different life histories 
lh_vec <- c("HAKE")#, "CRSNAP", "SIGSUT")
lh <- lapply(lh_vec, function(x) choose_lh_list(species=x, selex="asymptotic", param_adjust=c("R0"), val=c(1000), start_ages=1))
names(lh) <- lh_vec

## test across different true fishing mortality patterns 
Fdyn_set <- c("Endogenous", "Constant", "Ramp")

## test across different true recruitment patterns
Rdyn_set <- c("Constant")#, "Pulsed", "BH")

## test across different true levels of recruitment variation
SigmaR_set <- c(0, 0.3, 0.6, 0.9)

################################################
## Run cmsy models
################################################

## data availability scenarios -- LC currently not working, remove for now. 
avail_set <- c("catch", "catch_index", "catch_bsurvey", "catch_ML", "catch_index_ML", "catch_bsurvey_ML") 
avail_set_LC <- c("catch_LC", "catch_index_LC", "catch_bsurvey_LC")

## create combos
cmsy_modcombos <- as.matrix(expand.grid("Model"="CMSY", "Data_avail"=avail_set, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "SigmaR"=paste0("SigmaR_",SigmaR_set), "LH"=paste0("LH_", lh_vec)))

## transform model combinations into directories
cmsy_dir_vec <- model_paths(modcombos=cmsy_modcombos, res_dir=sim_dir)



#### no other data types included
### increasing fishing mortality, zero recruitment variation
model=1
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_catch_Fendog_SigmaR0.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()

png(file.path(fig_dir, "example_datagen_SigmaR0.6.png"), height=10, width=8, res=200, units="in")
par(mfrow=c(3,1), mar=c(0,0,0,0), omi=c(1,1.5,1,1))
  true1 <- readRDS(file.path(cmsy_dir_vec[37], 1, "True.rds"))
  true2 <- readRDS(file.path(cmsy_dir_vec[43], 1, "True.rds"))
  # true3 <- readRDS(file.path(cmsy_dir_vec[13], 1, "True.rds"))

plot(true1$F_t, ylim=c(0,0.5), type="l", lwd=4, xaxt="n", yaxt="n")
lines(true2$F_t, lwd=4, col="red", lty=2)
# lines(true3$F_t, lwd=4, col="red", lty=3)
axis(2, cex.axis=2, las=2)
mtext(side=2, "Fishing\nmortality", cex=2, line=4)
plot(true1$R_t, ylim=c(0,2000), type="l", lwd=4, xaxt="n", yaxt="n")
lines(true2$R_t, lwd=4, col="red", lty=2)
# lines(true3$R_t, lwd=4, col="red", lty=3)
axis(2, cex.axis=2, las=2)
mtext(side=2, "Recruitment", cex=2, line=6)
plot(true1$D_t, ylim=c(0,1), type="l", lwd=4, cex.axis=2, yaxt="n")
lines(true2$D_t, lwd=4, col="red", lty=2)
# lines(true3$D_t, lwd=4, col="red", lty=2)
axis(2, cex.axis=2, las=2)
mtext(side=2, "Depletion\n(B/K)", cex=2, line=4)
dev.off()




#### no other data types
### increasing fishing mortality, 0.6 recruitment variation
model=37
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_catch_Fendog_SigmaR0.6.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()

#### catch + mean length
### increasing fishing mortality, no recruitment variation
model=4
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_ML_Fendog_SigmaR0.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()

#### catch + index
model <- 2
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_index_Fendog_SigmaR0.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()

#### catch + biomass
model <- 3
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_biomass_Fendog_SigmaR0.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()

#### catch + index + mean length
model <- 5
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_indexML_Fendog_SigmaR0.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()

#### catch + index
## 0.6 recruitment variation
model <- 38
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_index_Fendog_SigmaR0.6.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()

#### catch + mean length
### 0.6 recruitment variation
model=40
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_ML_Fendog_SigmaR0.6.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()

#### no other data types included
### constant fishing mortality, zero recruitment variation
model=7
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_catch_Fconstant_SigmaR0.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()




model <- 1
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
plot(x=1, y=1, type="n", ylim=c(0,1), xlim=c(0,22))
for(i in out$idx){
  lines(as.numeric(out$ps.dt[i,]), col="red")
}
model <- 2
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
for(i in out$idx){
  lines(as.numeric(out$ps.dt[i,]), col="blue")
}








png(file.path(dir, "sampling_space_scatterhist.png"), height=10, width=12, res=200, units="in")
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  x <- out$S[,"msy"]
  y <- out$S[,"fmsy"]
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1), omi=c(0.5,0.5,1,1))
  plot(x,y, xlab="", ylab="", pch=19, col="gray", xaxs="i", yaxs="i", cex.lab=2, cex.axis=2)
  par(mar=c(0,3,1,1))
  plot(x=1, y=1, xlim=c(out$dfPriorInfo$par1[3], out$dfPriorInfo$par2[3]), ylim=c(0,top), type="n", axes=F, ann=F, xaxs="i", yaxs="i")
  # barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  polygon(x=c(out$dfPriorInfo$par1[3], out$dfPriorInfo$par2[3], out$dfPriorInfo$par2[3], out$dfPriorInfo$par1[3]), y=c(0,0,top/3,top/3), col="gray", border=NA)
  par(mar=c(3,0,1,1))
  plot(x=1, y=1, ylim=c(out$dfPriorInfo$par1[2], out$dfPriorInfo$par2[2]), xlim=c(0,top), type="n", axes=F, ann=F, xaxs="i", yaxs="i")
  # barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  polygon(x=c(0,0,top/3,top/3), y=c(out$dfPriorInfo$par1[2], out$dfPriorInfo$par2[2], out$dfPriorInfo$par2[2], out$dfPriorInfo$par1[2]), col="gray", border=NA)
  par(oma=c(3,3,0,0))
  mtext("MSY", side=1, line=1, outer=TRUE, adj=0, 
    at=.8 * (mean(x) - min(x))/(max(x)-min(x)), cex=2)
  mtext("Fmsy", side=2, line=1, outer=TRUE, adj=0, 
    at=(.8 * (mean(y) - min(y))/(max(y) - min(y))), cex=2)
dev.off()

png(file.path(dir, "narrow_scatterhist.png"), height=10, width=12, res=200, units="in")
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  x <- out$S[,"msy"]
  y <- out$S[,"fmsy"]
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  xhist2 <- hist(out$S[out$idx,"msy"], plot=FALSE)
  yhist2 <- hist(out$S[out$idx,"fmsy"], plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  
  par(mar=c(3,3,1,1), omi=c(0.5,0.5,1,1))
  out_cols <- rep("gray", nrow(out$S))
  out_cols[which(as.numeric(unlist(out$code))==0)] <- "black"
  plot(x,y, xlab="", ylab="", pch=19, col=out_cols, xaxs="i", yaxs="i", cex.lab=2, cex.axis=2)
  
  par(mar=c(0,3,1,1))
  plot(x=1, y=1, xlim=c(out$dfPriorInfo$par1[3], out$dfPriorInfo$par2[3]), ylim=c(0,top), type="n", axes=F, ann=F, xaxs="i", yaxs="i")
  polygon(y=c(0,0,top/3,top/3), x=c(out$dfPriorInfo$par1[3], out$dfPriorInfo$par2[3], out$dfPriorInfo$par2[3], out$dfPriorInfo$par1[3]), col="gray", border=NA)
  par(new=TRUE)
  barplot(c(rep(0,length(seq((xhist2$breaks[2]-xhist2$breaks[1]),min(xhist2$breaks)-(xhist2$breaks[2]-xhist2$breaks[1]), by=(xhist2$breaks[2]-xhist2$breaks[1])))),xhist2$counts), axes=FALSE, ylim=c(0, top), space=0, xaxs="i", yaxs="i", col="black")
  par(mar=c(3,0,1,1))
  plot(x=1, y=1, ylim=c(out$dfPriorInfo$par1[2], out$dfPriorInfo$par2[2]), xlim=c(0,top), type="n", axes=F, ann=F, xaxs="i", yaxs="i")
  polygon(x=c(0,0,top/3,top/3), y=c(out$dfPriorInfo$par1[2], out$dfPriorInfo$par2[2], out$dfPriorInfo$par2[2], out$dfPriorInfo$par1[2]), col="gray", border=NA)
  par(new=TRUE)
  barplot(c(rep(0,length(seq(0,0.02,by=0.01))), yhist2$counts, rep(0,length(seq(0.16,0.24,by=0.01)))), axes=FALSE, xlim=c(0,top), space=0, xaxs="i", yaxs="i", col="black", horiz=TRUE)

  par(oma=c(3,3,0,0))
  mtext("MSY", side=1, line=1, outer=TRUE, adj=0, 
    at=.8 * (mean(x) - min(x))/(max(x)-min(x)), cex=2)
  mtext("Fmsy", side=2, line=1, outer=TRUE, adj=0, 
    at=(.8 * (mean(y) - min(y))/(max(y) - min(y))), cex=2)
dev.off()


png(file.path(dir, "depl_catchonly.png"), height=10, width=12, res=200, units="in")
par(mfrow=c(1,1), mar=c(6,6,2,2))
plot(x=1,y=1,type="n",xlim=c(0,20),ylim=c(0,1.1),xaxs="i",yaxs="i",cex.lab=2,cex.axis=2,xlab="Year",ylab="Depletion (B/K)")
for(i in 1:nrow(out$ps.dt)){
	lines(as.numeric(out$ps.dt[i,]), col="red", lwd=2)
}
for(i in out$idx){
	lines(as.numeric(out$ps.dt[i,]), lwd=2)
}
dev.off()

png(file.path(dir, "depl_catchindex.png"), height=10, width=12, res=200, units="in")
par(mfrow=c(1,1), mar=c(6,6,2,2))
plot(x=1,y=1,type="n",xlim=c(0,20),ylim=c(0,1.1),xaxs="i",yaxs="i",cex.lab=2,cex.axis=2,xlab="Year",ylab="Depletion (B/K)")
for(i in 1:nrow(out_ind$ps.dt)){
	lines(as.numeric(out_ind$ps.dt[i,]), col="red", lwd=2)
}
for(i in out_ind$idx){
	lines(as.numeric(out_ind$ps.dt[i,]), lwd=2)
}
dev.off()

par(mfrow=c(1,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
plot(true_finc$F_t/max(true_finc$F_t), type="l", lwd=4, yaxt="n", xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", ylim=c(0,1.2))
plot(true_fcon$F_t/max(true_fcon$F_t), type="l", lwd=4, yaxt="n", xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", ylim=c(0,1.2))
plot(true_framp$F_t/max(true_framp$F_t), type="l", lwd=4, yaxt="n", xaxs="i", yaxs="i", cex.axis=2, xlab="", ylab="", ylim=c(0,1.2))
mtext(side=1, "Year", cex=2, line=4, outer=TRUE)
mtext(side=2, "Relative fishing mortality", cex=2, line=2, outer=TRUE)



plot((out$mlobs - out$mlexp)^2)