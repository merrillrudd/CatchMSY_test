rm(list=ls())

##############################################################
##### ------------- install packages ------------------- #####
##############################################################

devtools::install_github("merrillrudd/catchMSY", build.vignettes=TRUE, dependencies=TRUE)
library(catchMSY)

##############################################################
##### --------------- directories -----------------------#####
##############################################################

# main_dir <- "F:\\Merrill\\Git_Projects\\CatchMSY_test"
main_dir <- "C:\\Git_Projects\\CatchMSY_test"
source(file.path(main_dir, "R", "test_functions.R"))

## location of external data
data_dir <- file.path(main_dir, "inst", "extdata")

## setup directory for figures
fig_dir <- file.path(main_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)

##############################################################
##### -----------  setup hake demo 	   ----------------- #####
##############################################################
## from hake demo
nsamp <- 5000
ncores <- 1

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
hake$sel1 <- 4.0
hake$sel2 <- 5.0

## change recruitment variation
# hake$sigma_r <- 0.6

## age at 50% and 95% selectivity
selexPriorInfo <- data.frame("id"=4, "dist"="lnorm", "par1"=log(hake$sel1), "par2"=0.1*hake$sel1, "log"=TRUE, "stringAsFactors"=FALSE)
hake$dfPriorInfo <- rbind.data.frame(hake$dfPriorInfo, selexPriorInfo)

# hake$smodel <- "dome"

## test that model runs
hake_test <- catchMSYModel(hake)

##############################################################
##### -----------  self-testing 	   ----------------- #####
##############################################################
## simulate mean length and length composition
LC0 <- t(hake_test$LF)
ML0 <- hake_test$ML

set.seed(123)
hakeOM1 <- hake
hakeOM1$la.cv <- 0.07
hakeOM1$sigma_r <- 0.6
hake1 <- catchMSYModel(hakeOM1)
LC1 <- t(hake1$LF)
ML1 <- hake1$ML

set.seed(123)
hakeOM2 <- hake
hakeOM2$la.cv <- 0.14
hakeOM2$sigma_r <- 0.6
hake2 <- catchMSYModel(hakeOM2)
LC2 <- t(hake2$LF)
ML2 <- hake2$ML

##############################################################
##### -----------  model to use 	   ----------------- #####
##### growth variation = 0.14		   ----------------- #####
##### recruitment variation  = 0.6	   ----------------- #####
##### observation error mean lenth = 0.6 --------------- #####
##############################################################
hakeOM <- hakeOM2
hakeOM$sigma_r <- 0
hake_init <- catchMSYModel(hakeOM)
LC <- t(hake_init$LF)
ML <- hake_init$ML

# Generate random samples from dfPriorInfo
hakeOM <- sample.sid(sID=hakeOM, selex=TRUE, n=nsamp)
colnames(hakeOM$S) <- c("M", "Fmsy", "MSY", "sel50")

#### ------------- catch-only method -------------------#####
M0 <- hakeOM

# year and catch data only
M0$data <- M0$data[,c("year","catch")]

# run model with each sample
M0      <- sir.sid(M0, selex=FALSE, ncores)
M0_noSX <- sir.sid(M0, selex=FALSE, ncores)

# Get MSY statistics
M0$msy.stats <- summary(M0$S[M0$idx,3])
M0_noSX$msy.stats <- summary(M0_noSX$S[M0_noSX$idx,3])

#### ------------- catch + index ------------------------#####
M1 <- hakeOM

# year, catch, and index
M1$data <- M1$data[,c("year","catch","index","index.lse")]

# run model with each sample
M1 <- sir.sid(M1, selex=TRUE, ncores)
M1_noSX <- sir.sid(M1, selex=FALSE, ncores)

# Get MSY statistics
M1$msy.stats <- summary(M1$S[M1$idx,3])
M1_noSX$msy.stats <- summary(M1_noSX$S[M1_noSX$idx,3])

##############################################################
##### -----------  further development ----------------- #####
##### -----------  demonstrate new routines------------- #####
##############################################################


#### ------------- catch + mean length ------------------#####

M2 <- hakeOM
M2$data <- cbind(M2$data[,c("year","catch")], "avgSize"=ML, "avgSize.lse"=rep(0.2, length(ML)))

M2 <- sir.sid(M2,selex=TRUE,ncores)
M2_noSX <- sir.sid(M2, selex=FALSE, ncores)

# Get MSY statistics
M2$msy.stats <- summary(M2$S[M2$idx,3])

#### ------------- catch + mean length ------------------#####
### higher observation error!
M2_v2 <- hakeOM
M2_v2$data <- cbind(M2_v2$data[,c("year","catch")], "avgSize"=ML, "avgSize.lse"=rep(0.6, length(ML)))

M2_v2 <- sir.sid(M2_v2,selex=TRUE,ncores)
M2_v2_noSX <- sir.sid(M2_v2, selex=FALSE, ncores)

# Get MSY statistics
M2_v2$msy.stats <- summary(M2_v2$S[M2_v2$idx,3])


#### ------------- catch + mean length + index ------------------#####

M3 <- hakeOM

# year, catch, and index
M3$data <- cbind(M3$data[,c("year","catch","index","index.lse")], "avgSize"=ML, "avgSize.lse"=rep(0.6, length(ML))) 

# run model with each sample
M3 <- sir.sid(M3, selex=TRUE, ncores)
M3_noSX <- sir.sid(M3, selex=FALSE, ncores)

# Get MSY statistics
M3$msy.stats <- summary(M3$S[M3$idx,3])
M3_noSX <- summary(M3_noSX$S[M3_noSX$idx,3])

#### ---------- catch + length composition --------------#####
### ess cannot go beyond about 10
M4 <- hakeOM
M4$data <- cbind(M4$data[,c("year","catch")], "lencomp_sd"=rep(0.2, length(M0$data$year)), LC)

# run model with each sample
M4 <- sir.sid(M4,selex=TRUE,ncores)
M4_noSX <- sir.sid(M4, selex=FALSE, ncores)

# Get MSY statistics
M4$msy.stats <- summary(M4$S[M4$idx,3])

#### ---------- catch + length composition + index -----------#####

M5 <- hakeOM
M5$data <- cbind(M5$data[,c("year","catch", "index", "index.lse")], LC)

# run model with each sample
M5 <- sir.sid(M5,selex=TRUE,ncores)

# Get MSY statistics
M5$msy.stats <- summary(M5$S[M5$idx,3])


##############################################################
##### -----------  figures 	   ----------------- #####
##############################################################

### for report -- all models in boxplot
png(file.path(fig_dir, "boxplots_distr.png"), height=8, width=15, res=200, units="in")
boxplot(hakeOM$S[,"MSY"], M0$S[M0$idx,"MSY"], M1$S[M1$idx,"MSY"], M2$S[M2$idx,"MSY"], M3$S[M3$idx,"MSY"], M4$S[M4$idx,"MSY"], col=c("gray","red", "blue", "orange", "turquoise", "forestgreen"), lwd=2, xlim=c(0.5,6.5), ylim=c(hakeOM$dfPriorInfo$par1[3], hakeOM$dfPriorInfo$par2[3]*1.4), cex.label=1.5)
axis(1, at=1:6, labels=c("Prior", "Catch only", "Index", "Mean Length", "Mean Length+Index", "Length comp"), cex.axis=1.5)
text(x=1, y=500, paste0(round(median(M0$S[,"MSY"]),0), "\n (", round(quantile(M0$S[,"MSY"], 0.05),0), ",", round(quantile(M0$S[,"MSY"], 0.975),0), ")"), col=gray(0.2), font=2, cex=1.6)
text(x=2, y=500, paste0(round(median(M0$S[M0$idx,"MSY"]),0), "\n (", round(quantile(M0$S[M0$idx,"MSY"], 0.05),0), ",", round(quantile(M0$S[M0$idx,"MSY"], 0.975),0), ")"), col="red", font=2, cex=1.6)
text(x=3, y=500, paste0(round(median(M1$S[M1$idx,"MSY"]),0), "\n (", round(quantile(M1$S[M1$idx,"MSY"], 0.05),0), ",", round(quantile(M1$S[M1$idx,"MSY"], 0.975),0), ")"), col="blue", font=2, cex=1.6)
text(x=4, y=500, paste0(round(median(M2$S[M2$idx,"MSY"]),0), "\n (", round(quantile(M2$S[M2$idx,"MSY"], 0.05),0), ",", round(quantile(M2$S[M2$idx,"MSY"], 0.975),0), ")"), col="orange", font=2, cex=1.6)
text(x=5, y=500, paste0(round(median(M3$S[M3$idx,"MSY"]),0), "\n (", round(quantile(M3$S[M3$idx,"MSY"], 0.05),0), ",", round(quantile(M3$S[M3$idx,"MSY"], 0.975),0), ")"), col="turquoise", font=2, cex=1.6)
text(x=6, y=500, paste0(round(median(M4$S[M4$idx,"MSY"]),0), "\n (", round(quantile(M4$S[M4$idx,"MSY"], 0.05),0), ",", round(quantile(M4$S[M4$idx,"MSY"], 0.975),0), ")"), col="forestgreen", font=2, cex=1.6)
mtext("MSY", side=2, line=2.75, cex=1.5)
mtext("Data included in model", side=1, line=3, cex=1.5)
dev.off()

### for report -- all models in boxplot
png(file.path(fig_dir, "boxplots_distr_sel50.png"), height=8, width=15, res=200, units="in")
boxplot(hakeOM$S[,"sel50"], M0$S[M0$idx,"sel50"], M1$S[M1$idx,"sel50"], M2$S[M2$idx,"sel50"], M3$S[M3$idx,"sel50"], M4$S[M4$idx,"sel50"], col=c("gray","red", "blue", "orange", "turquoise", "forestgreen"), lwd=2, xlim=c(0.5,6.5), ylim=c(0, 20), cex.label=1.5)
axis(1, at=1:6, labels=c("Prior", "Catch only", "Index", "Mean Length", "Mean Length+Index", "Length comp"), cex.axis=1.5)
text(x=1, y=18, paste0(round(median(M0$S[,"sel50"]),2), "\n (", round(quantile(M0$S[,"sel50"], 0.05),2), ",", round(quantile(M0$S[,"sel50"], 0.975),2), ")"), col=gray(0.2), font=2, cex=1.6)
text(x=2, y=18, paste0(round(median(M0$S[M0$idx,"sel50"]),2), "\n (", round(quantile(M0$S[M0$idx,"sel50"], 0.05),2), ",", round(quantile(M0$S[M0$idx,"sel50"], 0.975),2), ")"), col="red", font=2, cex=1.6)
text(x=3, y=18, paste0(round(median(M1$S[M1$idx,"sel50"]),2), "\n (", round(quantile(M1$S[M1$idx,"sel50"], 0.05),2), ",", round(quantile(M1$S[M1$idx,"sel50"], 0.975),2), ")"), col="blue", font=2, cex=1.6)
text(x=4, y=18, paste0(round(median(M2$S[M2$idx,"sel50"]),2), "\n (", round(quantile(M2$S[M2$idx,"sel50"], 0.05),2), ",", round(quantile(M2$S[M2$idx,"sel50"], 0.975),2), ")"), col="orange", font=2, cex=1.6)
text(x=5, y=18, paste0(round(median(M3$S[M3$idx,"sel50"]),2), "\n (", round(quantile(M3$S[M3$idx,"sel50"], 0.05),2), ",", round(quantile(M3$S[M3$idx,"sel50"], 0.975),2), ")"), col="turquoise", font=2, cex=1.6)
text(x=6, y=18, paste0(round(median(M4$S[M4$idx,"sel50"]),2), "\n (", round(quantile(M4$S[M4$idx,"sel50"], 0.05),2), ",", round(quantile(M4$S[M4$idx,"sel50"], 0.975),2), ")"), col="forestgreen", font=2, cex=1.6)
mtext("Age at 50% selectivity", side=2, line=2.75, cex=1.5)
mtext("Data included in model", side=1, line=3, cex=1.5)
dev.off()


### plot data
png(file.path(fig_dir, "hake_data_plot.png"), width=12, height=8, res=200, units="in")
par(mfrow=c(2,2), mar=c(4,5,0,0), omi=c(1,1,1,1))
plot(hake$data$year, hake$data$catch, type="l", lwd=3, ylim=c(0, max(hake$data$catch)*1.2), cex.axis=1.5, xaxs="i", yaxs="i", xlab="", ylab="")
mtext(side=2, "Catch", cex=1.5, line=3)
plot(hake$data$year, hake$data$index, type="l", lwd=3, ylim=c(0, max(hake$data$index)*1.2), cex.axis=1.5, xaxs="i", yaxs="i", xlab="", ylab="")
mtext(side=2, "Index", cex=1.5, line=3)

plot(ML, ylim=c(0, max(ML)*1.2), xaxs="i", yaxs="i", cex=1.6, lwd=2, type="o", pch=17, xlab="Year", ylab="", cex.lab=1.5, cex.axis=1.5, xaxt="n")
axis(1, at=seq(1,25,by=5), labels=c(1965, 1970, 1975, 1980, 1985), cex.axis=1.5)
mtext(side=2, "Mean Length", cex=1.5, line=3)

plot(LC[nrow(LC),], pch=17, cex=1.3, cex.axis=1.5, xaxs="i", yaxs="i", xlab="", ylab="")
mtext(side=1, "Length bin (cm)", cex=1.5, line=3)
mtext(side=2, "Frequency", cex=1.5, line=3)
dev.off()



## residuals
## index per year
plot(x=1, y=1, type="n", xlim=c(0,(ncol(M1$index_resid))), ylim=c(min(M1$index_resid),max(M1$index_resid)), xaxs="i", yaxs="i")
for(i in 1:ncol(M1$index_resid)){
	points(x=rep(i,nrow(M1$index_resid)), y=M1$index_resid[,i])
}
abline(h=0, col="red")

## residuals - index from all years
matplot(M1$index_resid, pch=19, ylim=c(-1.5,3), col="#0000AA30")
abline(h=0, col="red", lwd=5)

## residuals - length comp
find_lc <- M4$lc_resid[M4$idx]
choose_lc <- find_lc[[1]]

### one sample
plot(x=1, y=1, type="n", xlim=c(0,ncol(choose_lc)), ylim=c(0,nrow(choose_lc)+2), xaxs="i", yaxs="i")
for(i in 1:nrow(choose_lc)){
	for(j in 1:ncol(choose_lc)){
		col <- ifelse(choose_lc[i,j]<0, "#0000AA50", "#AA000050")
		points(x=j,y=i,cex=abs(as.numeric(choose_lc[i,j]))/15, col=col, pch=19)
	}
}

## median of all samples
lc_citers <- matrix(NA, nrow=nrow(choose_lc), ncol=ncol(choose_lc))
for(i in 1:nrow(choose_lc)){
	for(j in 1:ncol(choose_lc))
		lc_citers[i,j] <- median(sapply(1:length(find_lc), function(x) as.numeric(find_lc[[x]][i,j])))
}
plot(x=1, y=1, type="n", xlim=c(0,ncol(lc_citers)), ylim=c(0,nrow(lc_citers)+2), xaxs="i", yaxs="i", xlab="", ylab="")
for(i in 1:nrow(lc_citers)){
	for(j in 1:ncol(lc_citers)){
		col <- ifelse(lc_citers[i,j]<0, "#0000AA50", "#AA000050")
		points(x=j,y=i,cex=abs(as.numeric(lc_citers[i,j]))/15, col=col, pch=19)
	}
}
mtext(side=1, "Year", cex=1.3, line=3)
mtext(side=2, "Length bin", cex=1.3, line=3)

## across years
lc_cyears <- matrix(NA, nrow=nrow(choose_lc), ncol=length(find_lc))
for(i in 1:nrow(choose_lc)){
	for(j in 1:length(find_lc)){
		lc_cyears[i,j] <- median(unlist(find_lc[[j]][i,]))
	}
}
plot(x=1, y=1, type="n", xlim=c(0,1000), ylim=c(0,nrow(lc_cyears)+2), xaxs="i", yaxs="i", xlab="", ylab="")
for(i in 1:1000){
	for(j in 1:ncol(lc_cyears)){
		col <- ifelse(lc_cyears[i,j]<0, "#0000AA50", "#AA000050")
		points(x=j,y=i,cex=abs(as.numeric(lc_cyears[i,j]))/50, col=col, pch=19)
	}
}
mtext(side=1, "Iteration", cex=1.3, line=3)
mtext(side=2, "Length bin", cex=1.3, line=3)

### combine index and lc residual plots
par(mfrow=c(1,2), mar=c(1,4,1,1), omi=c(1,1,1,1))
matplot(M1$index_resid, pch=19, ylim=c(-1.5,3), col="#33333350", xlab="", ylab="")
abline(h=0, col="red", lwd=5)
mtext(side=1, "Iteration", cex=1.5, line=3)
mtext(side=2, "Residual", cex=1.5, line=3)
plot(x=1, y=1, type="n", xlim=c(0,ncol(choose_lc)), ylim=c(0,nrow(choose_lc)+2), xaxs="i", yaxs="i", xlab="", ylab="")
for(i in 1:nrow(choose_lc)){
	for(j in 1:ncol(choose_lc)){
		col <- ifelse(choose_lc[i,j]<0, "#0000AA50", "#AA000050")
		points(x=j,y=i,cex=abs(as.numeric(choose_lc[i,j]))/15, col=col, pch=19)
	}
}
mtext(side=1, "Year", cex=1.5, line=3)
mtext(side=2, "Length bin", cex=1.5, line=3)




### prior distributions
par(mfrow=c(1,1))
hist(rlnorm(nsamp, log(hake$m), 0.05*hake$m), col="gray", main="", xlim=c(0.14,0.16), xaxs="i", yaxs="i")
text(x=0.142, y=800, "M", font=2, cex=3)

par(mfrow=c(1,1))
hist(runif(nsamp, 0.20*hake$m, 1.5*hake$m), col="gray", main="", xaxs="i", yaxs="i", xlim=c(0,0.25))
text(x=0.025, y=300, "Fmsy", font=2, cex=3, xpd=NA)

par(mfrow=c(1,1))
hist(runif(nsamp, quantile(hake$data$catch,0.05), quantile(hake$data$catch,0.95)), col="gray", main="", xaxs="i", yaxs="i", xlim=c(0,500))
text(x=50, y=350, "MSY", font=2, cex=3, xpd=NA)

par(mfrow=c(1,1))
hist(rlnorm(nsamp, log(hake$sel1), 0.1*hake$sel1), col="gray", main="", xaxs="i", yaxs="i")
text(x=15, y=1200, "sel50", font=2, cex=3, xpd=NA)



## self-generated mean length - no recruitment variation
png(file.path(fig_dir, "mean_length_lowCVgrowth.png"), height=6, width=8, res=200, units="in")
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(ML0, ylim=c(0, max(ML1)*1.2), xaxs="i", yaxs="i", cex=1.6, lwd=2, type="o", pch=17, xlab="Year", ylab="Mean length", cex.lab=1.5, cex.axis=1.3)
lines(ML1, type="o", pch=19, col="blue", lwd=2)
dev.off()

## self-generated length comp - no recruitment variation
png(file.path(fig_dir, "length_comp_lowCVgrowth.png"), height=10, width=12, res=200, units="in")
par(mfrow=c(5,5), mar=c(0,0,0,0), omi=c(1,1,0.2,0.2))
for(i in 1:nrow(LC1)){
	plot(LC0[i,], xaxt="n", yaxt="n", ylim=c(0,50), pch=19)
	lines(LC1[i,],lwd=2, col="blue")
	if(i %in% seq(1,23,by=5)) axis(2, las=2)
	if(i %in% 21:23) axis(1)
}
mtext(side=1, "Length bin (cm)", cex=1.5, outer=TRUE, line=3)
mtext(side=2, "Frequency", cex=1.5, outer=TRUE, line=3)
dev.off()

## self-generated mean length - recruitment variation
png(file.path(fig_dir, "mean_length_SigmaR.png"), height=6, width=8, res=200, units="in")
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(ML0, ylim=c(0, max(ML2)*1.2), xaxs="i", yaxs="i", cex=1.6, lwd=2, type="o", pch=17, xlab="Year", ylab="Mean length", cex.lab=1.5, cex.axis=1.3)
lines(ML2, type="o", pch=19, col="blue", lwd=2)
dev.off()

## self-generated length comp - recruitment variation
png(file.path(fig_dir, "length_comp_SigmaR.png"), height=10, width=12, res=200, units="in")
par(mfrow=c(5,5), mar=c(0,0,0,0), omi=c(1,1,0.2,0.2))
for(i in 1:nrow(LC2)){
	plot(LC0[i,], xaxt="n", yaxt="n", ylim=c(0,50), pch=19)
	lines(LC2[i,],lwd=2, col="blue")
	if(i %in% seq(1,23,by=5)) axis(2, las=2)
	if(i %in% 21:23) axis(1)
}
mtext(side=1, "Length bin (cm)", cex=1.5, outer=TRUE, line=3)
mtext(side=2, "Frequency", cex=1.5, outer=TRUE, line=3)
dev.off()


## sampling space scatterplot
png(file.path(fig_dir, "sampling_space_scatter.png"), width=10, height=8, units="in", res=200)
pairs(hakeOM$S, gap=0, pch=20, cex.axis=1.3)
dev.off()

## sampling space histogram
png(file.path(fig_dir, "MSY_sampling_space_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
legend("topright", legend=c("Sampling space"), pch=15, col="gray")
dev.off()

png(file.path(fig_dir, "sel50_sampling_space_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
legend("topright", legend=c("Sampling space"), pch=15, col="gray")
dev.off()

## sampling space boxplot
png(file.path(fig_dir, "MSY_boxplot_1.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"MSY"], col="gray", lwd=2, xlim=c(0,8), ylim=c(hakeOM$dfPriorInfo$par1[3], hakeOM$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_1.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"sel50"], col="gray", lwd=2, xlim=c(0,8), ylim=c(0,exp(hakeOM$dfPriorInfo$par1[4]*2)))
dev.off()


# catch-only scatterplot 
M0_cols <- rep("black", nsamp)
M0_cols[which(as.numeric(unlist(M0$code))>0)] <- "red"
png(file.path(fig_dir, "catch_only_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M0$S, gap=0, col=M0_cols,pch=20)
dev.off()

## catch-only histogram
png(file.path(fig_dir, "MSY_catchonly_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only"), pch=15, col=c("gray", "#AA000050"))
dev.off()

png(file.path(fig_dir, "MSY_catchonlyNoSX_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0_noSX$S[M0_noSX$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only"), pch=15, col=c("gray", "#AA000050"))
dev.off()

png(file.path(fig_dir, "sel50_catchonly_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only"), pch=15, col=c("gray", "#AA000050"))
dev.off()

png(file.path(fig_dir, "sel50_catchonlyNoSX_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0_noSX$S[M0_noSX$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only"), pch=15, col=c("gray", "#AA000050"))
dev.off()

## catch-only boxplot
png(file.path(fig_dir, "MSY_boxplot_2.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"MSY"], M0$S[M0$idx,"MSY"], col=c("gray","red"), lwd=2, xlim=c(0,8), ylim=c(hakeOM$dfPriorInfo$par1[3], hakeOM$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_2.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"sel50"], M0$S[M0$idx,"sel50"], col=c("gray","red"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hakeOM$dfPriorInfo$par1[4]*2)))
dev.off()


# catch + index scatterplot
M1_cols <- rep("black", nsamp)
M1_cols[which(1:nsamp %in% unique(M1$idx)==FALSE)] <- "steelblue"
M1_cols[which(as.numeric(unlist(M0$code))>0)] <- "red"
png(file.path(fig_dir, "catch_index_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M1$S, gap=0, col=M1_cols,pch=20)
dev.off()

## catch + index histogram
png(file.path(fig_dir, "MSY_catchindex_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"MSY"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index"), pch=15, col=c("gray", "#AA000050", "#0000AA50"))
dev.off()

png(file.path(fig_dir, "MSY_catchindexNoSX_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0_noSX$S[M0_noSX$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1_noSX$S[M1_noSX$idx,"MSY"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index"), pch=15, col=c("gray", "#AA000050", "#0000AA50"))
dev.off()

png(file.path(fig_dir, "sel50_catchindex_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"sel50"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index"), pch=15, col=c("gray", "#AA000050", "#0000AA50"))
dev.off()

png(file.path(fig_dir, "sel50_catchindexNoSX_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0_noSX$S[M0_noSX$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1_noSX$S[M1_noSX$idx,"sel50"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index"), pch=15, col=c("gray", "#AA000050", "#0000AA50"))
dev.off()

## catch + index boxplot
png(file.path(fig_dir, "MSY_boxplot_3.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"MSY"], M0$S[M0$idx,"MSY"], M1$S[M1$idx,"MSY"], col=c("gray","red","steelblue"), lwd=2, xlim=c(0,8), ylim=c(hakeOM$dfPriorInfo$par1[3], hakeOM$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_3.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"sel50"], M0$S[M0$idx,"sel50"], M1$S[M1$idx,"sel50"], col=c("gray","red","steelblue"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hakeOM$dfPriorInfo$par1[4]*2)))
dev.off()


# catch + mean length scatterplot
M2_cols <- rep("black", nsamp)
M2_cols[which(1:nsamp %in% unique(M2$idx)==FALSE)] <- "violet"
M2_cols[which(as.numeric(unlist(M0$code))>0)] <- "red"
png(file.path(fig_dir, "catch_meanlen_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M2$S, gap=0, col=M2_cols,pch=20)
dev.off()

## catch + mean length histogram
png(file.path(fig_dir, "MSY_catchmeanlen_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M2$S[M2$idx,"MSY"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+MeanLength"), pch=15, col=c("gray", "#AA000050", "#00AA0050"))
dev.off()

png(file.path(fig_dir, "sel50_catchmeanlen_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M2$S[M2$idx,"sel50"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+MeanLength"), pch=15, col=c("gray", "#AA000050", "#00AA0050"))
dev.off()


## catch + mean length boxplot
png(file.path(fig_dir, "MSY_boxplot_4.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"MSY"], M0$S[M0$idx,"MSY"], M2$S[M2$idx,"MSY"], col=c("gray","red","forestgreen"), lwd=2, xlim=c(0,8), ylim=c(hakeOM$dfPriorInfo$par1[3], hakeOM$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_4.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"sel50"], M0$S[M0$idx,"sel50"], M2$S[M2$idx,"sel50"],col=c("gray","red","forestgreen"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hakeOM$dfPriorInfo$par1[4]*2)))
dev.off()


# catch + mean length low obs error scatterplot
M2_v2_cols <- rep("black", nsamp)
M2_v2_cols[which(1:nsamp %in% unique(M2_v2$idx)==FALSE)] <- "violet"
M2_v2_cols[which(as.numeric(unlist(M0$code))>0)] <- "red"
png(file.path(fig_dir, "catch_meanlen_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M2_v2$S, gap=0, col=M2_v2_cols,pch=20)
dev.off()

## catch + mean length low obs error histogram
png(file.path(fig_dir, "MSY_catchmeanlen_LowObsErr_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M2_v2$S[M2_v2$idx,"MSY"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+MeanLength"), pch=15, col=c("gray", "#AA000050", "#00AA0050"))
dev.off()

png(file.path(fig_dir, "sel50_catchmeanlen_LowObsErr_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M2_v2$S[M2_v2$idx,"sel50"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+MeanLength"), pch=15, col=c("gray", "#AA000050", "#00AA0050"))
dev.off()


## catch + mean length low obs error boxplot
png(file.path(fig_dir, "MSY_boxplot_4_v2.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"MSY"], M0$S[M0$idx,"MSY"], M2_v2$S[M2_v2$idx,"MSY"], col=c("gray","red","forestgreen"), lwd=2, xlim=c(0,8), ylim=c(hakeOM$dfPriorInfo$par1[3], hakeOM$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_4_v2.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"sel50"], M0$S[M0$idx,"sel50"], M2_v2$S[M2_v2$idx,"sel50"],col=c("gray","red","forestgreen"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hakeOM$dfPriorInfo$par1[4]*2)))
dev.off()


# catch + avgSize + index scatterplot
M3_cols <- rep("black", nsamp)
M3_cols[which(1:nsamp %in% unique(M3$idx)==FALSE)] <- "purple"
M3_cols[which(1:nsamp %in% unique(M1$idx)==FALSE)] <- "steelblue"
M3_cols[which(as.numeric(unlist(M0$code))>0)] <- "red"
png(file.path(fig_dir, "catch_indexmeanlen_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M3$S, gap=0, col=M3_cols,pch=20)
dev.off()

## catch + avgSize + index histogram
png(file.path(fig_dir, "MSY_catchindexmeanlen_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"MSY"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M3$S[M3$idx,"MSY"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space","Catch+Index", "Catch+Index+MeanLength"), pch=15, col=c("gray",  "#0000AA50", "#00AA0050"))
dev.off()

png(file.path(fig_dir, "sel50_catchindexmeanlen_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"sel50"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M3$S[M3$idx,"sel50"], col="#00AA0050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space","Catch+Index", "Catch+Index+MeanLength"), pch=15, col=c("gray",  "#0000AA50", "#00AA0050"))
dev.off()

## catch + avgSize + index boxplot
png(file.path(fig_dir, "MSY_boxplot_5.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"MSY"], M1$S[M1$idx,"MSY"],M3$S[M3$idx,"MSY"], col=c("gray","steelblue", "forestgreen"), lwd=2, xlim=c(0,8), ylim=c(hakeOM$dfPriorInfo$par1[3], hakeOM$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_5.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"sel50"], M1$S[M1$idx,"sel50"], M3$S[M3$idx,"sel50"], col=c("gray","steelblue","forestgreen"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hakeOM$dfPriorInfo$par1[4]*2)))
dev.off()


# catch + length comp scatterplot
M4_cols <- rep("black", nsamp)
M4_cols[which(1:nsamp %in% unique(M4$idx)==FALSE)] <- "green"
M4_cols[which(as.numeric(unlist(M0$code))>0)] <- "red"
png(file.path(fig_dir, "catch_lencomp_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M4$S, gap=0, col=M4_cols,pch=20)
dev.off()

## catch + length comp histogram
png(file.path(fig_dir, "MSY_catchlencomp_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M4$S[M4$idx,"MSY"], col="#AA00FF50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+LengthComp"), pch=15, col=c("gray", "#AA000050","#AA00FF50"))
dev.off()

png(file.path(fig_dir, "MSY_catchlencompNoSX_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M4_noSX$S[M4_noSX$idx,"MSY"], col="#AA00FF50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+LengthComp"), pch=15, col=c("gray", "#AA000050","#AA00FF50"))
dev.off()

png(file.path(fig_dir, "sel50_catchlencomp_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M4$S[M4$idx,"sel50"], col="#AA00FF50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+LengthComp"), pch=15, col=c("gray", "#AA000050","#AA00FF50"))
dev.off()


## catch + length comp boxplot
png(file.path(fig_dir, "MSY_boxplot_6.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"MSY"], M0$S[M0$idx,"MSY"], M4$S[M4$idx,"MSY"], col=c("gray","red","purple"), lwd=2, xlim=c(0,8), ylim=c(hakeOM$dfPriorInfo$par1[3], hakeOM$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_6.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"sel50"], M0$S[M0$idx,"sel50"], M4$S[M4$idx,"sel50"], col=c("gray","red","purple"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hakeOM$dfPriorInfo$par1[4]*2)))
dev.off()

# catch + lengthcomp + index scatterplot
M5_cols <- rep("black", nsamp)
M5_cols[which(1:nsamp %in% unique(M5$idx)==FALSE)] <- "forestgreen"
M5_cols[which(1:nsamp %in% unique(M1$idx)==FALSE)] <- "steelblue"
M5_cols[which(as.numeric(unlist(M0$code))>0)] <- "red"
png(file.path(fig_dir, "catch_indexlencomp_scatter.png"), width=10, height=8, units="in", res=200)
pairs(M5$S, gap=0, col=M5_cols,pch=20)
dev.off()

## catch + lengthcomp + index histogram
png(file.path(fig_dir, "MSY_catchindexlencomp_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, hakeOM$dfPriorInfo$par2[3]*1.2)
ylim <- c(0, nsamp/10)
hist(hakeOM$S[,"MSY"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="MSY", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"MSY"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"MSY"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M5$S[M5$idx,"MSY"], col="#00AAAA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index", "Catch+Index+LengthComp"), pch=15, col=c("gray", "#AA000050","#0000AA50", "#00AAAA50"))
dev.off()

png(file.path(fig_dir, "sel50_catchindexlencomp_hist.png"), width=10, height=6, units="in", res=200)
xlim <- c(0, exp(hakeOM$dfPriorInfo$par1[4])*3)
ylim <- c(0, nsamp)
hist(hakeOM$S[,"sel50"], col="gray", lty="blank", xlim=xlim, ylim=ylim, xlab="sel50", ylab="Frequency", main="")
par(new=TRUE)
hist(M0$S[M0$idx,"sel50"], col="#AA000050", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M1$S[M1$idx,"sel50"], col="#0000AA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
par(new=TRUE)
hist(M5$S[M5$idx,"sel50"], col="#00AAAA50", lty="blank", xlim=xlim, ylim=ylim, xlab="", ylab="", main="")
legend("topright", legend=c("Sampling space", "Catch only", "Catch+Index", "Catch+Index+LengthComp"), pch=15, col=c("gray", "#AA000050","#0000AA50", "#00AAAA50"))
dev.off()

## catch + lengthcomp + index boxplot
png(file.path(fig_dir, "MSY_boxplot_7.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"MSY"], M0$S[M0$idx,"MSY"], M1$S[M1$idx,"MSY"], M2$S[M2$idx,"MSY"], M3$S[M3$idx,"MSY"], M4$S[M4$idx,"MSY"], M5$S[M5$idx,"MSY"], col=c("gray","red","steelblue", "violet", "purple", "lawngreen", "forestgreen"), lwd=2, xlim=c(0,8), ylim=c(hakeOM$dfPriorInfo$par1[3], hakeOM$dfPriorInfo$par2[3]))
dev.off()

png(file.path(fig_dir, "sel50_boxplot_7.png"), width=10, height=6, units="in", res=200)
boxplot(hakeOM$S[,"sel50"], M0$S[M0$idx,"sel50"], M1$S[M1$idx,"sel50"], M2$S[M2$idx,"sel50"], M3$S[M3$idx,"sel50"], M4$S[M4$idx,"sel50"], M5$S[M5$idx,"sel50"], col=c("gray","red","steelblue", "violet", "purple", "lawngreen", "forestgreen"), lwd=2, xlim=c(0,8), ylim=c(0,exp(hakeOM$dfPriorInfo$par1[4]*2)))
dev.off()

#### qqplots
png(file.path(fig_dir, "qqMSY.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(2,2))
qqplot(hakeOM$S[,"MSY"], M0$S[M0$idx,"MSY"], xlab="Sampling space", ylab="Posterior", xlim=c(100,400), ylim=c(100,400), main="Catch")
lines(x=100:400, y=100:400, col="blue", lwd=3)

qqplot(hakeOM$S[,"MSY"], M1$S[M1$idx,"MSY"], xlab="Sampling space", ylab="Posterior", xlim=c(100,400), ylim=c(100,400), main="+Index")
lines(x=100:400, y=100:400, col="blue", lwd=3)

qqplot(hakeOM$S[,"MSY"], M2$S[M2$idx,"MSY"], xlab="Sampling space", ylab="Posterior", xlim=c(100,400), ylim=c(100,400), main="+Mean Length")
lines(x=100:400, y=100:400, col="blue", lwd=3)

qqplot(hakeOM$S[,"MSY"], M3$S[M3$idx,"MSY"], xlab="Sampling space", ylab="Posterior", xlim=c(100,400), ylim=c(100,400), main="+Index + Mean Length")
lines(x=100:400, y=100:400, col="blue", lwd=3)
dev.off()
