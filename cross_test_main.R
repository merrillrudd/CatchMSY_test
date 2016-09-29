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
SigmaR_set <- c(0, 0.3, 0.6, 0.9)

################################################
## Run cmsy models
################################################

## data availability scenarios -- LC currently not working, remove for now. 
avail_set <- c("catch", "catch_index", "catch_bsurvey", "catch_ML", "catch_index_ML", "catch_bsurvey_ML") 
avail_set_LC <- c("catch_LC", "catch_index_LC", "catch_bsurvey_LC")

da <- list("Nyears"=20, "Nyears_comp"=20, "comp_sample"=1000) 

## create combos
cmsy_modcombos <- as.matrix(expand.grid("Model"="CMSY", "Data_avail"=avail_set, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "SigmaR"=paste0("SigmaR_",SigmaR_set), "LH"=paste0("LH_", lh_vec)))

## transform model combinations into directories
cmsy_dir_vec <- model_paths(modcombos=cmsy_modcombos, res_dir=sim_dir)

## run iterations
itervec <- 1


##--------------------- setup parallel ----------------------
registerDoParallel(cores=5)

## ------------------ simulate data -------------------------

start_datagen <- Sys.time()

## create true population and generated data into directories
foreach(loop=1:length(cmsy_dir_vec), .packages=c('LIME','catchMSY')) %dopar% generateData(modpath=cmsy_dir_vec[loop], itervec=itervec, spatial=TRUE, Fdynamics=strsplit(cmsy_modcombos[loop,"Fdyn"],"_")[[1]][2], Rdynamics=strsplit(cmsy_modcombos[loop,"Rdyn"],"_")[[1]][2], LType=1, write=TRUE, lh_list=lh, data_avail_list=da, modname=paste0(cmsy_modcombos[loop,"Model"],"_",cmsy_modcombos[loop,"Data_avail"]), rewrite=FALSE, param_adjust="SigmaR", val=as.numeric(strsplit(cmsy_modcombos[loop,"SigmaR"],"_")[[1]][2]))

end_datagen <- Sys.time() - start_datagen

## -------------- run estimation models -----------------------

## catchMSY
start_run <- Sys.time()

foreach(loop=1:length(cmsy_dir_vec), .packages=c('LIME', 'catchMSY')) %dopar% tryCatch(run_cmsy(modpath=cmsy_dir_vec[loop], itervec=itervec, lh_list=lh, data_avail=cmsy_modcombos[loop,"Data_avail"], nyears=20, rewrite=TRUE), error=function(e) print(paste0("issue with ", cmsy_dir_vec[loop])))

end_run <- Sys.time() - start_run

## -------------- figures  -----------------------
compare_re(dir_vec=cmsy_dir_vec, mod_names=avail_set, Fdyn_vec=Fdyn_set, Rdyn_vec=Rdyn_set, lh_num=lh_vec, save=TRUE, fig_name="RE_cmsy")



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

foreach(loop=1:length(bdcmsy_dir_vec), .packages=c("catchMSY", "gplots")) %dopar% run_cmsy_bd(modpath=bdcmsy_dir_vec[loop], itervec=itervec, lh_list=lh, rewrite=FALSE, nsamp=5000, ncores=ncores)

end_run <- Sys.time() - start_run

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

