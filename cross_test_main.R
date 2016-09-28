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
itervec <- 1:20


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

# for(loop in 1:length(cmsy_dir_vec)){
# 	run_cmsy(modpath=cmsy_dir_vec[loop], itervec=itervec, lh_list=lh, data_avail=cmsy_modcombos[loop,"Data_avail"], nyears=20, rewrite=FALSE)
# }

foreach(loop=1:length(cmsy_dir_vec), .packages=c('LIME', 'catchMSY')) %dopar% tryCatch(run_cmsy(modpath=cmsy_dir_vec[loop], itervec=itervec, lh_list=lh, data_avail=cmsy_modcombos[loop,"Data_avail"], nyears=20, rewrite=FALSE), error=function(e) print(paste0("issue with ", cmsy_dir_vec[loop])))

end_run <- Sys.time() - start_run

## -------------- figures  -----------------------
compare_re(dir_vec=cmsy_dir_vec, mod_names=avail_set, Fdyn_vec=Fdyn_set, Rdyn_vec=Rdyn_set, lh_num=lh_vec, save=TRUE, fig_name="RE_cmsy")



################################################
## Sensitivity - dome
################################################

## test across different life histories 
lh_vec <- c("CRSNAP", "SIGSUT", "HAKE")
lh_dome <- lapply(lh_vec, function(x) choose_lh_list(species=x, selex="dome", param_adjust=c("R0"), val=c(1000), start_ages=1))
names(lh_dome) <- lh_vec

## data availability scenarios -- LC currently not working, remove for now. 
nolc_avail_set <- c("catch", "catch_index", "catch_bsurvey", "catch_ML", "catch_index_ML", "catch_bsurvey_ML") 
# avail_set_LC <- c("catch_LC", "catch_index_LC", "catch_bsurvey_LC")

da <- list("Nyears"=20, "Nyears_comp"=20, "comp_sample"=1000) 

## create combos
dome_modcombos <- as.matrix(expand.grid("Model"="CMSY_dome", "Data_avail"=nolc_avail_set, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "LH"=paste0("LH_", lh_vec)))

## transform model combinations into directories
dome_cmsy_dirs <- model_paths(modcombos=dome_modcombos, res_dir=sim_dir)

## run iterations
itervec <- 1:20


##--------------------- setup parallel ----------------------
registerDoParallel(cores=5)

## ------------------ simulate data -------------------------

start_datagen <- Sys.time()

## create true population and generated data into directories
foreach(loop=1:length(dome_cmsy_dirs), .packages=c('LIME','catchMSY')) %dopar% generateData(modpath=dome_cmsy_dirs[loop], itervec=itervec, spatial=TRUE, Fdynamics=strsplit(cmsy_modcombos[loop,"Fdyn"],"_")[[1]][2], Rdynamics=strsplit(cmsy_modcombos[loop,"Rdyn"],"_")[[1]][2], LType=1, write=TRUE, lh_list=lh_dome, data_avail_list=da, modname=paste0(cmsy_modcombos[loop,"Model"],"_",cmsy_modcombos[loop,"Data_avail"]), rewrite=FALSE)

end_datagen <- Sys.time() - start_datagen

## -------------- run estimation models -----------------------

## catchMSY
start_run <- Sys.time()

foreach(loop=1:length(dome_cmsy_dirs), .packages=c('LIME', 'catchMSY')) %dopar% tryCatch(run_cmsy(modpath=dome_cmsy_dirs[loop], itervec=itervec, lh_list=lh, data_avail=cmsy_modcombos[loop,"Data_avail"], nyears=20, rewrite=FALSE), error=function(e) print(paste0("issue with ", dome_cmsy_dirs[loop])))

end_run <- Sys.time() - start_run


## -------------- figures  -----------------------
compare_re(dir_vec=dome_cmsy_dirs, mod_names=nolc_avail_set, Fdyn_vec=Fdyn_set, Rdyn_vec=Rdyn_set, lh_num="CRSNAP", save=FALSE, fig_name="test")

# SPRcover <- interval_coverage(modpath_vec=mixe_dirs, param="SPR", itervec=itervec)
# compare_fits(dir_vec=mixe_dirs, mod_names=paste0("LC", c(1,10)), Ftype="Ramp", Rtype="Pulsed", fig_name="Fits_LC", LHchoose=5, iter=7, save=TRUE, coverage=SPRcover$pcover, type="all")




################################################
## Compare to LIME method
################################################

## data availability scenarios -- LC currently not working, remove for now. 
lc_avail_set <- c("catch_LC", "catch_index_LC", "catch_bsurvey_LC") 

da <- list("Nyears"=20, "Nyears_comp"=20, "comp_sample"=1000)

## create combos
lc_modcombos <- as.matrix(expand.grid("Model"=c("LIME","CMSY"), "Data_avail"=lc_avail_set, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "LH"=paste0("LH_", lh_vec)))

## transform model combinations into directories
lc_dir_vec <- model_paths(modcombos=lc_modcombos, res_dir=sim_dir)
lc_cmsy_dirs <- lc_dir_vec[grep("CMSY", lc_dir_vec)]
lc_lime_dirs <- lc_dir_vec[grep("LIME", lc_dir_vec)]

