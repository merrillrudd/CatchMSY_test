rm(list=ls())

##############################################################
##### ------------- install packages ------------------- #####
##############################################################
devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)

devtools::install_github("merrillrudd/catchMSY", build.vignettes=TRUE, dependencies=TRUE)
library(catchMSY)
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", dep=TRUE)
library(TMBhelper)
library(foreach)
library(doParallel)

##############################################################
##### --------------- directories -----------------------#####
##############################################################

main_dir <- "F:\\Merrill\\Git_Projects\\CatchMSY_test"
# main_dir <- "C:\\Git_Projects\\CatchMSY_test"
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


##############################################################
##### --------- simulation-testing    ------------------ #####
##############################################################

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
SigmaR_set <- c(0, 0.6)

################################################
## Run cmsy models
################################################

## data availability scenarios -- LC currently not working, remove for now. 
avail_set <- c("catch", "catch_index", "catch_bsurvey", "catch_ML", "catch_index_ML", "catch_bsurvey_ML") 
avail_set_LC <- c("catch_LC", "catch_index_LC", "catch_bsurvey_LC")

da <- list("Nyears"=20, "Nyears_comp"=20, "comp_sample"=1000) 

## create combos
cmsy_modcombos <- as.matrix(expand.grid("Model"="CMSY", "Data_avail"=avail_set, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "SigmaR"=paste0("SigmaR_",SigmaR_set), "LH"=paste0("LH_", lh_vec)))
lc_modcombos <- as.matrix(expand.grid("Model"="CMSY", "Data_avail"=avail_set_LC, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "SigmaR"=paste0("SigmaR_",SigmaR_set), "LH"=paste0("LH_", lh_vec)))


## transform model combinations into directories
cmsy_dir_vec <- model_paths(modcombos=cmsy_modcombos, res_dir=sim_dir)
lc_dir_vec <- model_paths(modcombos=lc_modcombos, res_dir=sim_dir)

itervec <- 1:50


##--------------------- setup parallel ----------------------
registerDoParallel(cores=5)

## ------------------ simulate data -------------------------

start_datagen <- Sys.time()

## create true population and generated data into directories
foreach(loop=1:length(cmsy_dir_vec), .packages=c('LIME','catchMSY')) %dopar% generateData(modpath=cmsy_dir_vec[loop], itervec=itervec, spatial=TRUE, Fdynamics=strsplit(cmsy_modcombos[loop,"Fdyn"],"_")[[1]][2], Rdynamics=strsplit(cmsy_modcombos[loop,"Rdyn"],"_")[[1]][2], LType=1, write=TRUE, lh_list=lh, data_avail_list=da, modname=paste0(cmsy_modcombos[loop,"Model"],"_",cmsy_modcombos[loop,"Data_avail"]), rewrite=TRUE, param_adjust="SigmaR", val=as.numeric(strsplit(cmsy_modcombos[loop,"SigmaR"],"_")[[1]][2]))

end_datagen <- Sys.time() - start_datagen

## -------------- run estimation models -----------------------

## catchMSY
start_run <- Sys.time()

foreach(loop=1:length(cmsy_dir_vec), .packages=c('LIME', 'catchMSY')) %dopar% tryCatch(run_cmsy(modpath=cmsy_dir_vec[loop], itervec=itervec, lh_list=lh, data_avail=cmsy_modcombos[loop,"Data_avail"], nyears=20, rewrite=FALSE), error=function(e) print(paste0("issue with ", cmsy_dir_vec[loop])))

modpath=cmsy_dir_vec[loop]
itervec=itervec
lh_list=lh
data_avail=cmsy_modcombos[loop,"Data_avail"]
nyears=20
rewrite=FALSE

end_run <- Sys.time() - start_run

######## length comp models #################################
## ------------------ simulate data -------------------------

start_datagen <- Sys.time()

## create true population and generated data into directories
foreach(loop=1:length(lc_dir_vec), .packages=c('LIME','catchMSY')) %dopar% generateData(modpath=lc_dir_vec[loop], itervec=itervec, spatial=TRUE, Fdynamics=strsplit(lc_modcombos[loop,"Fdyn"],"_")[[1]][2], Rdynamics=strsplit(lc_modcombos[loop,"Rdyn"],"_")[[1]][2], LType=1, write=TRUE, lh_list=lh, data_avail_list=da, modname=paste0(lc_modcombos[loop,"Model"],"_",lc_modcombos[loop,"Data_avail"]), rewrite=TRUE, param_adjust="SigmaR", val=as.numeric(strsplit(lc_modcombos[loop,"SigmaR"],"_")[[1]][2]))

end_datagen <- Sys.time() - start_datagen

## -------------- run estimation models -----------------------

## catchMSY
start_run <- Sys.time()

foreach(loop=1:length(lc_dir_vec), .packages=c('LIME', 'catchMSY')) %dopar% tryCatch(run_cmsy(modpath=lc_dir_vec[loop], itervec=itervec, lh_list=lh, data_avail=lc_modcombos[loop,"Data_avail"], nyears=20, rewrite=TRUE), error=function(e) print(paste0("issue with ", lc_dir_vec[loop])))

end_run <- Sys.time() - start_run




## -------------- figures  -----------------------
compare_re(dir_vec=cmsy_dir_vec, mod_names=avail_set, Fdyn_vec=Fdyn_set, Rdyn_vec=Rdyn_set, SigmaR_vec=c(0,0.6), lh_num=lh_vec, save=TRUE, fig_name="RE_cmsy")
compare_re(dir_vec=cmsy_dir_vec[which(grepl("SigmaR_0/", cmsy_dir_vec))], mod_names=avail_set, Fdyn_vec="Endogenous", Rdyn_vec=Rdyn_set, SigmaR_vec=0, lh_num=lh_vec, save=TRUE, fig_name="RE_cmsy_Fendog_Sigma0")
compare_re(dir_vec=cmsy_dir_vec[which(grepl("SigmaR_0/", cmsy_dir_vec))], mod_names=avail_set, Fdyn_vec=Fdyn_set, Rdyn_vec=Rdyn_set, SigmaR_vec=0, lh_num=lh_vec, save=TRUE, fig_name="RE_cmsy_Fcompare_Sigma0")



model <- 1
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fendog_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 2
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_index_Fendog_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 4
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_ml_Fendog_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()


model <- 7
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fconstant_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 10
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_ml_Fconstant_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()


model <- 19
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fendog_sigR0.6.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 25
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fconstant_sigR0.6.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 28
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_ml_Fconstant_sigR0.6.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()




## biomass dynamics
model <- 1
reoutbd <- sapply(itervec, function(x) re_calc_bd(bd=bdcmsy_dir_vec[model], iter=x))
covoutbd <- sapply(itervec, function(x) cover_calc_bd(bd=bdcmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catchBD_Fendogenous_sigR0.png"), height=10, width=15, res=200, units="in")
plot(x=1,y=1,type="n",ylim=c(-2,15),xlim=c(0,21), cex.lab=2, xlab="Iteration", ylab="Relative Error")
for(i in 1:ncol(reoutbd)){
	segments(x0=i, y0=as.numeric(reoutbd["relcl",i]), y1=as.numeric(reoutbd["reucl",i]), col="gray", lwd=4)
}
points(x=1:ncol(reoutbd), y=reoutbd["re",], col=covoutbd, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4), pch=19, cex=1.5)
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reoutbd)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covoutbd)=="gray"))/length(covoutbd),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 1
reout <- sapply(1:50, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(1:50, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
reoutbd <- sapply(1:20, function(x) re_calc_bd(bd=bdcmsy_dir_vec[model], iter=x))
covoutbd <- sapply(1:20, function(x) cover_calc_bd(bd=bdcmsy_dir_vec[model], iter=x))
png(file.path(fig_dir,"MSY_iters_compareAgeBD_catch_Fendog_sigR0.png"), height=12, width=15, res=200, units="in")
par(mfrow=c(2,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,15), xlim=c(0,21), xaxt="n")
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.8,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.8,0.8), cex=1.5)
plot(x=1,y=1,type="n",ylim=c(-2,15),xlim=c(0,21), cex.lab=2, xlab="Iteration", ylab="Relative Error")
for(i in 1:ncol(reoutbd)){
	segments(x0=i, y0=as.numeric(reoutbd["relcl",i]), y1=as.numeric(reoutbd["reucl",i]), col="gray", lwd=4)
}
points(x=1:ncol(reoutbd), y=reoutbd["re",], col=covoutbd, xlab="Iteration", ylab="Relative Error", cex.lab=2, pch=19, cex=1.5)
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reoutbd)),2)), xy=c(0.8,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covoutbd)=="gray"))/length(covoutbd),2)), xy=c(0.8,0.8), cex=1.5)
mtext("Iteration", side=1, cex=2, line=3)
mtext("Relative error", side=2, cex=2, line=3, outer=TRUE)
dev.off()

model <- 2
reout <- sapply(itervec, function(x) re_calc_bd(bd=bdcmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc_bd(bd=bdcmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_indexBD_Fendogenous_sigR0.png"), height=10, width=15, res=200, units="in")
plot(x=1,y=1,type="n",ylim=c(-2,15),xlim=c(0,21), cex.lab=2, xlab="Iteration", ylab="Relative Error")
for(i in 1:ncol(reout)){
	segments(x0=i, y0=as.numeric(reout["relcl",i]), y1=as.numeric(reout["reucl",i]), col="gray", lwd=4)
}
points(x=1:ncol(reout), y=reout["re",], col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4), pch=19, cex=1.5)
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 2
reout <- sapply(1:50, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(1:50, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
reoutbd <- sapply(1:20, function(x) re_calc_bd(bd=bdcmsy_dir_vec[model], iter=x))
covoutbd <- sapply(1:20, function(x) cover_calc_bd(bd=bdcmsy_dir_vec[model], iter=x))
png(file.path(fig_dir,"MSY_iters_compareAgeBD_index_Fendog_sigR0.png"), height=12, width=15, res=200, units="in")
par(mfrow=c(2,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,10), xlim=c(0,21), xaxt="n")
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.8,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.8,0.8), cex=1.5)
plot(x=1,y=1,type="n",ylim=c(-2,10),xlim=c(0,21), cex.lab=2, xlab="Iteration", ylab="Relative Error")
for(i in 1:ncol(reoutbd)){
	segments(x0=i, y0=as.numeric(reoutbd["relcl",i]), y1=as.numeric(reoutbd["reucl",i]), col="gray", lwd=4)
}
points(x=1:ncol(reoutbd), y=reoutbd["re",], col=covoutbd, xlab="Iteration", ylab="Relative Error", cex.lab=2, pch=19, cex=1.5)
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reoutbd)),2)), xy=c(0.8,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covoutbd)=="gray"))/length(covoutbd),2)), xy=c(0.8,0.8), cex=1.5)
mtext("Iteration", side=1, cex=2, line=3)
mtext("Relative error", side=2, cex=2, line=3, outer=TRUE)
dev.off()


################################################
## Compare to biomass dynamic method
################################################
library(R2jags)  # Interface with JAGS
library(coda) 
library("parallel")
library("foreach")
library("doParallel")
library("gplots")

ncores <- 5
FullSchaefer <- F    # initialize variable; automatically set to TRUE if enough abundance data are available
n.chains     <- ifelse(ncores > 2,3,2) # set 3 chains in JAGS if more than 2 cores are available
ncores_for_computation <- ncores # cores to be used for parallel processing of CMSY
cl           <- makeCluster(ncores_for_computation)
registerDoParallel(cl, cores = ncores_for_computation)

bd_avail_set <- c("catch", "catch_index", "catch_bsurvey") 

da <- list("Nyears"=20, "Nyears_comp"=20, "comp_sample"=1000) 


bdcmsy_modcombos <- as.matrix(expand.grid("Model"="BD_CMSY", "Data_avail"=bd_avail_set, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "SigmaR"=paste0("SigmaR_",SigmaR_set), "LH"=paste0("LH_", lh_vec)))

## transform model combinations into directories
bdcmsy_dir_vec <- model_paths(modcombos=bdcmsy_modcombos, res_dir=sim_dir)

## run iterations
itervec <- 1:20

##--------------------- setup parallel ----------------------
registerDoParallel(cores=ncores)

## ------------------ simulate data -------------------------

start_datagen <- Sys.time()


## create true population and generated data into directories
foreach(loop=1:length(bdcmsy_dir_vec), .packages=c('LIME','catchMSY')) %dopar% generateData(modpath=bdcmsy_dir_vec[loop], itervec=itervec, spatial=TRUE, Fdynamics=strsplit(bdcmsy_modcombos[loop,"Fdyn"],"_")[[1]][2], Rdynamics=strsplit(bdcmsy_modcombos[loop,"Rdyn"],"_")[[1]][2], LType=1, write=TRUE, lh_list=lh, data_avail_list=da, modname=paste0(bdcmsy_modcombos[loop,"Model"],"_",bdcmsy_modcombos[loop,"Data_avail"]), rewrite=FALSE, param_adjust="SigmaR", val=as.numeric(strsplit(bdcmsy_modcombos[loop,"SigmaR"],"_")[[1]][2]))

end_datagen <- Sys.time() - start_datagen

## -------------- run estimation models -----------------------

start_run <- Sys.time()

foreach(loop=1:length(bdcmsy_dir_vec), .packages=c("catchMSY", "gplots", "R2jags", "coda")) %dopar% run_cmsy_bd(modpath=bdcmsy_dir_vec[loop], itervec=itervec, lh_list=lh, rewrite=FALSE, nsamp=5000, ncores=ncores)

end_run <- Sys.time() - start_run




dir <- bdcmsy_dir_vec[1]
dirx <- cmsy_dir_vec[1]
msy_out <- sapply(1:20, function(x) readMSY(bd=dir, cmsy=dirx, iter=x))


dir <- bdcmsy_dir_vec[2]
dirx <- cmsy_dir_vec[2]
index_out <- sapply(1:20, function(x) readMSY(bd=dir, cmsy=dirx, iter=x))


boxplot(unlist(msy_out["rebd",]), unlist(msy_out["reage",]), unlist(index_out["rebd",]), unlist(index_out["reage",]), col=c(rep(c("gray", "turquoise"),2)))
abline(h=0, lty=2, lwd=4)

dir_vec=cmsy_dir_vec[-c(which(grepl("SigmaR_0.3/", cmsy_dir_vec)), which(grepl("SigmaR_0.9/", cmsy_dir_vec)), which(grepl("ML", cmsy_dir_vec)))]
mod_names=avail_set
Fdyn_vec="Endogenous"
Rdyn_vec=Rdyn_set
SigmaR_vec=c(0,0.6)
lh_num=lh_vec
save=TRUE
fig_name="RE_compare_agebio"
bd_dir_vec <- bdcmsy_dir_vec
itervec <- 1:20

compare_re(dir_vec=cmsy_dir_vec[which(grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("ML", cmsy_dir_vec)==FALSE)], mod_names=avail_set[-which(grepl("ML",avail_set))], Fdyn_vec="Endogenous", Rdyn_vec=Rdyn_set, SigmaR_vec=c(0.6), lh_num=lh_vec, save=TRUE, fig_name="RE_compare_agebio_Fendog_Sigma0.6", bd_dir_vec=bdcmsy_dir_vec[grepl("SigmaR_0.6/",bdcmsy_dir_vec)], itervec=1:20, ylim=c(-2.5,6))



################################################
## Compare to LIME method
################################################
library(LIME)

## data availability scenarios -- LC currently not working, remove for now. 
lc_avail_set <- c("catch_LC", "catch_index_LC", "catch_bsurvey_LC") 
LIME_names <- c("Catch_LC20", "Rich_LC")

da <- data_avail_settings(avail_set=LIME_names, ESS=1000)

## create combos
lime_modcombos <- as.matrix(expand.grid("Model"=c("LIME"), "Data_avail"=LIME_names, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "SigmaR"=paste0("SigmaR_", SigmaR_set), "LH"=paste0("LH_", lh_vec)))

## transform model combinations into directories
lime_dir_vec <- model_paths(modcombos=lime_modcombos, res_dir=sim_dir)

