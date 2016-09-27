rm(list=ls())

##############################################################
##### ------------- install packages ------------------- #####
##############################################################

devtools::install_github("merrillrudd/catchMSY", build.vignettes=TRUE, dependencies=TRUE)
library(catchMSY)
library(dirichlet)

##############################################################
##### --------------- directories -----------------------#####
##############################################################

main_dir <- "F:\\Merrill\\Git_Projects\\CatchMSY_test"
# main_dir <- "C:\\Git_Projects\\CatchMSY_test"
source(file.path(main_dir, "R", "test_functions.R"))

## location of external data
data_dir <- file.path(main_dir, "inst", "extdata")

## setup directory for figures
fig_dir <- file.path(main_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)



##############################################################
##### -----------  demonstrations 	   ----------------- #####
##### -----------  self-testing 		   ------------- #####
##############################################################
## from hake demo
nsamp <- 1000
ncores <- 1

## test that model runs
catchMSYModel(hake)

# Set parameter sampling frame
hake$dfPriorInfo$dist[1] = "lnorm"
hake$dfPriorInfo$par1[1] = log(hake$m)
hake$dfPriorInfo$par2[1] = 0.05 * hake$m
hake$dfPriorInfo$dist[2] = "unif"
hake$dfPriorInfo$par1[2] = 0.20 * hake$m
hake$dfPriorInfo$par2[2] = 1.50 * hake$m
hake$dfPriorInfo$dist[3] = "unif"
hake$dfPriorInfo$par1[3] = quantile(hake$data$catch,0.05)
hake$dfPriorInfo$par2[3] = quantile(hake$data$catch,0.95)

# Change selectivity
hake$sel50 <- 4.0
hake$sel95 <- 5.0

## change recruitment variation
hake$sigma_r <- 0.6

## change ages
hake$age <- 1:30


## age at 50% and 95% selectivity
selexPriorInfo <- data.frame("id"=4, "dist"="lnorm", "par1"=log(hake$sel50), "par2"=0.1*hake$sel50, "log"=TRUE, "stringAsFactors"=FALSE)
hake$dfPriorInfo <- rbind.data.frame(hake$dfPriorInfo, selexPriorInfo)

# Generate ransome samples from dfPriorInfo
hake <- sample.sid(sID=hake, selex=TRUE, n=nsamp)
colnames(hake$S) <- c("M", "Fmsy", "MSY", "sel50")

## simulate mean length and length composition
OM <- hake
OM$la.cv <- 0.10
OM_new <- catchMSYModel(OM)
LC <- t(OM_new$LF)
ML <- OM_new$ML


png(file.path(fig_dir, "sampling_space_scatter.png"), width=10, height=8, units="in", res=200)
pairs(hake$S, gap=0, pch=20, cex.axis=1.3)
dev.off()

## histogram
png(file.path(fig_dir, "MSY_sampling_space_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hake$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hake$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
legend("topright", legend=c("Sampling space"), pch=15, col="gray")
dev.off()

png(file.path(fig_dir, "sel50_sampling_space_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hake$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hake$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
legend("topright", legend=c("Sampling space"), pch=15, col="gray")
dev.off()

## boxplot
png(file.path(fig_dir, "MSY_boxplot_1.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"MSY"], col="gray", lwd=2, xlim=c(0,8), ylim=c(hake$dfPriorInfo$par1[3], hake$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_1.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"sel50"], col="gray", lwd=2, xlim=c(0,8), ylim=c(0,exp(hake$dfPriorInfo$par1[4]*2)))
dev.off()


#### ------------- catch-only method -------------------#####
M0 <- hake

# year and catch data only
M0$data <- M0$data[,c("year","catch")]

# run model with each sample
M0      <- sir.sid(M0, selex=FALSE, ncores)

# Get MSY statistics
M0$msy.stats <- summary(M0$S[M0$idx,3])

# Narrow down samples
M0_cols <- rep("black", nsamp)
M0_cols[which(as.numeric(unlist(M0$code))>0)] <- "goldenrod"
png(file.path(fig_dir, "catch_only_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M0$S, gap=0, col=M0_cols,pch=20)
dev.off()

## histogram
png(file.path(fig_dir, "MSY_catchonly_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hake$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hake$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only"), pch=15, col=c("gray", "#AA000050"))
dev.off()

png(file.path(fig_dir, "sel50_catchonly_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hake$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hake$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only"), pch=15, col=c("gray", "#AA000050"))
dev.off()

## boxplot
png(file.path(fig_dir, "MSY_boxplot_2.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"MSY"], M0$S[M0$idx,"MSY"], col=c("gray","goldenrod"), lwd=2, xlim=c(0,8), ylim=c(hake$dfPriorInfo$par1[3], hake$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_2.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"sel50"], M0$S[M0$idx,"sel50"], col=c("gray","goldenrod"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hake$dfPriorInfo$par1[4]*2)))
dev.off()



#### ------------- catch + index ------------------------#####
M1 <- hake

# year, catch, and index
M1$data <- M1$data[,c("year","catch","index","index.lse")]

# run model with each sample
M1 <- sir.sid(M1, selex=TRUE, ncores)

# Get MSY statistics
M1$msy.stats <- summary(M1$S[M1$idx,3])

# Narrow down samples
M1_cols <- rep("black", nsamp)
M1_cols[which(1:nsamp %in% unique(M1$idx)==FALSE)] <- "steelblue"
M1_cols[which(as.numeric(unlist(M0$code))>0)] <- "goldenrod"
png(file.path(fig_dir, "catch_index_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M1$S, gap=0, col=M1_cols,pch=20)
dev.off()

## histogram
png(file.path(fig_dir, "MSY_catchindex_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hake$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hake$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"MSY"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index"), pch=15, col=c("gray", "#AA000050", "#0000AA50"))
dev.off()

png(file.path(fig_dir, "sel50_catchindex_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hake$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hake$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"sel50"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index"), pch=15, col=c("gray", "#AA000050", "#0000AA50"))
dev.off()

## boxplot
png(file.path(fig_dir, "MSY_boxplot_3.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"MSY"], M0$S[M0$idx,"MSY"], M1$S[M1$idx,"MSY"], col=c("gray","goldenrod","steelblue"), lwd=2, xlim=c(0,8), ylim=c(hake$dfPriorInfo$par1[3], hake$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_3.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"sel50"], M0$S[M0$idx,"sel50"], M1$S[M1$idx,"sel50"], col=c("gray","goldenrod","steelblue"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hake$dfPriorInfo$par1[4]*2)))
dev.off()

##############################################################
##### -----------  further development ----------------- #####
##### -----------  demonstrate new routines------------- #####
##############################################################


#### ------------- catch + mean length ------------------#####

M2 <- hake
M2$data <- cbind(M2$data[,c("year","catch")], "meanlength"=ML, "meanlength.lse"=rep(0.6, length(ML)))

M2 <- sir.sid(M2,selex=TRUE,ncores)

# Get MSY statistics
M2$msy.stats <- summary(M2$S[M2$idx,3])

# Narrow down samples
M2_cols <- rep("black", nsamp)
M2_cols[which(1:nsamp %in% unique(M2$idx)==FALSE)] <- "violet"
M2_cols[which(as.numeric(unlist(M0$code))>0)] <- "goldenrod"
png(file.path(fig_dir, "catch_meanlen_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M2$S, gap=0, col=M2_cols,pch=20)
dev.off()

## histogram
png(file.path(fig_dir, "MSY_catchmeanlen_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hake$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hake$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M2$S[M2$idx,"MSY"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+MeanLength"), pch=15, col=c("gray", "#AA000050", "#00AA0050"))
dev.off()

png(file.path(fig_dir, "sel50_catchmeanlen_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hake$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hake$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M2$S[M2$idx,"sel50"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+MeanLength"), pch=15, col=c("gray", "#AA000050", "#00AA0050"))
dev.off()


## boxplot
png(file.path(fig_dir, "MSY_boxplot_4.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"MSY"], M0$S[M0$idx,"MSY"], M1$S[M1$idx,"MSY"], M2$S[M2$idx,"MSY"], col=c("gray","goldenrod","steelblue", "violet"), lwd=2, xlim=c(0,8), ylim=c(hake$dfPriorInfo$par1[3], hake$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_4.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"sel50"], M0$S[M0$idx,"sel50"], M1$S[M1$idx,"sel50"], M2$S[M2$idx,"sel50"], col=c("gray","goldenrod","steelblue", "violet"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hake$dfPriorInfo$par1[4]*2)))
dev.off()

#### ------------- catch + mean length + index ------------------#####

M3 <- hake

# year, catch, and index
M3$data <- cbind(M3$data[,c("year","catch","index","index.lse")], "meanlength"=ML, "meanlength.lse"=rep(0.6, length(ML))) 

# run model with each sample
M3 <- sir.sid(M3, ncores)

# Get MSY statistics
M3$msy.stats <- summary(M3$S[M3$idx,3])

# Narrow down samples
M3_cols <- rep("black", nsamp)
M3_cols[which(1:nsamp %in% unique(M3$idx)==FALSE)] <- "purple"
M3_cols[which(1:nsamp %in% unique(M1$idx)==FALSE)] <- "steelblue"
M3_cols[which(as.numeric(unlist(M0$code))>0)] <- "goldenrod"
png(file.path(fig_dir, "catch_indexmeanlen_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M3$S, gap=0, col=M3_cols,pch=20)
dev.off()

## histogram
png(file.path(fig_dir, "MSY_catchindexmeanlen_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hake$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hake$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"MSY"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M3$S[M3$idx,"MSY"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index", "Catch+Index+MeanLength"), pch=15, col=c("gray", "#AA000050", "#0000AA50", "#00AA0050"))
dev.off()

png(file.path(fig_dir, "sel50_catchindexmeanlen_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hake$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hake$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"sel50"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M3$S[M3$idx,"sel50"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index", "Catch+Index+MeanLength"), pch=15, col=c("gray", "#AA000050", "#0000AA50", "#00AA0050"))
dev.off()

## boxplot
png(file.path(fig_dir, "MSY_boxplot_5.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"MSY"], M0$S[M0$idx,"MSY"], M1$S[M1$idx,"MSY"], M2$S[M2$idx,"MSY"], M3$S[M3$idx,"MSY"], col=c("gray","goldenrod","steelblue", "violet", "purple"), lwd=2, xlim=c(0,8), ylim=c(hake$dfPriorInfo$par1[3], hake$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_5.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"sel50"], M0$S[M0$idx,"sel50"], M1$S[M1$idx,"sel50"], M2$S[M2$idx,"sel50"], M3$S[M3$idx,"sel50"], col=c("gray","goldenrod","steelblue", "violet", "purple"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hake$dfPriorInfo$par1[4]*2)))
dev.off()

#### ---------- catch + length composition --------------#####
### ess cannot go beyond about 10
M4 <- hake
M4$data <- cbind(M4$data[,c("year","catch")], LC)

# run model with each sample
M4 <- sir.sid(M4,selex=FALSE,ncores)

# Get MSY statistics
M4$msy.stats <- summary(M4$S[M4$idx,3])

# Narrow down samples
M4_cols <- rep("black", nsamp)
M4_cols[which(1:nsamp %in% unique(M4$idx)==FALSE)] <- "green"
M4_cols[which(as.numeric(unlist(M0$code))>0)] <- "goldenrod"
png(file.path(fig_dir, "catch_lencomp_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M4$S, gap=0, col=M4_cols,pch=20)
dev.off()

## histogram
png(file.path(fig_dir, "MSY_catchlencomp_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hake$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hake$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M4$S[M4$idx,"MSY"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+LengthComp"), pch=15, col=c("gray", "#AA000050","#00AA0050"))
dev.off()

png(file.path(fig_dir, "sel50_catchlencomp_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hake$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hake$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M4$S[M4$idx,"sel50"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+LengthComp"), pch=15, col=c("gray", "#AA000050","#00AA0050"))
dev.off()


## boxplot
png(file.path(fig_dir, "MSY_boxplot_6.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"MSY"], M0$S[M0$idx,"MSY"], M1$S[M1$idx,"MSY"], M2$S[M2$idx,"MSY"], M3$S[M3$idx,"MSY"], M4$S[M4$idx,"MSY"], col=c("gray","goldenrod","steelblue", "violet", "purple", "lawngreen"), lwd=2, xlim=c(0,8), ylim=c(hake$dfPriorInfo$par1[3], hake$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_6.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"sel50"], M0$S[M0$idx,"sel50"], M1$S[M1$idx,"sel50"], M2$S[M2$idx,"sel50"], M3$S[M3$idx,"sel50"], M4$S[M4$idx,"sel50"], col=c("gray","goldenrod","steelblue", "violet", "purple", "lawngreen"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hake$dfPriorInfo$par1[4]*2)))
dev.off()
#### ---------- catch + length composition + index -----------#####

M5 <- hake
M5$data <- cbind(M5$data[,c("year","catch", "index", "index.lse")], LC)

# run model with each sample
M5 <- sir.sid(M5,selex=TRUE,ncores)

# Get MSY statistics
M5$msy.stats <- summary(M5$S[M5$idx,3])

# Narrow down samples
M5_cols <- rep("black", nsamp)
M5_cols[which(1:nsamp %in% unique(M5$idx)==FALSE)] <- "forestgreen"
M5_cols[which(1:nsamp %in% unique(M1$idx)==FALSE)] <- "steelblue"
M5_cols[which(as.numeric(unlist(M0$code))>0)] <- "goldenrod"
png(file.path(fig_dir, "catch_indexlencomp_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M5$S, gap=0, col=M5_cols,pch=20)
dev.off()

## histogram
png(file.path(fig_dir, "MSY_catchindexlencomp_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hake$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hake$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"MSY"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M5$S[M5$idx,"MSY"], col="#00AAAA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index", "Catch+Index+LengthComp"), pch=15, col=c("gray", "#AA000050","#0000AA50", "#00AAAA50"))
dev.off()

png(file.path(fig_dir, "sel50_catchindexlencomp_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hake$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hake$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"sel50"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M5$S[M5$idx,"sel50"], col="#00AAAA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index", "Catch+Index+LengthComp"), pch=15, col=c("gray", "#AA000050","#0000AA50", "#00AAAA50"))
dev.off()

## boxplot
png(file.path(fig_dir, "MSY_boxplot_7.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"MSY"], M0$S[M0$idx,"MSY"], M1$S[M1$idx,"MSY"], M2$S[M2$idx,"MSY"], M3$S[M3$idx,"MSY"], M4$S[M4$idx,"MSY"], M5$S[M5$idx,"MSY"], col=c("gray","goldenrod","steelblue", "violet", "purple", "lawngreen", "forestgreen"), lwd=2, xlim=c(0,8), ylim=c(hake$dfPriorInfo$par1[3], hake$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_7.png"), width=10, height=6, units="in", res=200)
boxplot(hake$S[,"sel50"], M0$S[M0$idx,"sel50"], M1$S[M1$idx,"sel50"], M2$S[M2$idx,"sel50"], M3$S[M3$idx,"sel50"], M4$S[M4$idx,"sel50"], M5$S[M5$idx,"sel50"], col=c("gray","goldenrod","steelblue", "violet", "purple", "lawngreen", "forestgreen"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hake$dfPriorInfo$par1[4]*2)))
dev.off()

#### qqplots
png(file.path(fig_dir, "qqMSY.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(2,2))
qqplot(hake$S[,"MSY"], M0$S[M0$idx,"MSY"], xlab="Sampling space", ylab="Posterior", xlim=c(100,400), ylim=c(100,400), main="Catch")
lines(x=100:400, y=100:400, col="blue", lwd=3)

qqplot(hake$S[,"MSY"], M1$S[M1$idx,"MSY"], xlab="Sampling space", ylab="Posterior", xlim=c(100,400), ylim=c(100,400), main="+Index")
lines(x=100:400, y=100:400, col="blue", lwd=3)

qqplot(hake$S[,"MSY"], M2$S[M2$idx,"MSY"], xlab="Sampling space", ylab="Posterior", xlim=c(100,400), ylim=c(100,400), main="+Mean Length")
lines(x=100:400, y=100:400, col="blue", lwd=3)

qqplot(hake$S[,"MSY"], M3$S[M3$idx,"MSY"], xlab="Sampling space", ylab="Posterior", xlim=c(100,400), ylim=c(100,400), main="+Index + Mean Length")
lines(x=100:400, y=100:400, col="blue", lwd=3)
dev.off()
