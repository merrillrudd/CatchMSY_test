rm(list=ls())

devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
devtools::install_github("merrillrudd/catchMSY", build.vignettes=TRUE, dependencies=TRUE)
library(LIME)
library(catchMSY)
# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper", dep=TRUE)
library(TMBhelper)
library(foreach)
library(doParallel)

# main_dir <- "F:\\Merrill\\Git_Projects\\method_test"
main_dir <- "C:\\Git_Projects\\method_test"
source(file.path(main_dir, "R", "method_test_functions.R"))

## setup results directory
results_dir <- file.path(main_dir, "results")

## remove results (for testing)
# unlink(results_dir, TRUE)

## create directory
dir.create(results_dir, showWarnings=FALSE)

## setup directory for figures
fig_dir <- file.path(main_dir, "figs")
dir.create(fig_dir, showWarnings=FALSE)


#######################
## Settings
#######################

## test across different life histories 
lh_vec <- c("CRSNAP", "SIGSUT", "HAKE")
lh <- lapply(lh_vec, function(x) choose_lh_list(species=x, selex="asymptotic", param_adjust=c("R0"), val=c(1000), start_ages=1))
names(lh) <- lh_vec

## test across different true fishing mortality patterns 
Fdyn_set <- c("Endogenous", "Constant", "Ramp")

## test across different true recruitment patterns
Rdyn_set <- c("Constant", "Pulsed", "BH")


################################################
## Explore length composition data only
## and compare to LBPSR
################################################

model_set <- c("CMSY", "LIME", "LBSPR", "MLmort") 

## data availability scenarios -- LC currently not working, remove for now. 
avail_set <- c("catch", "catch_index", "catch_bsurvey", "catch_LC", "catch_index_LC", "catch_bsurvey_LC", "catch_ML", "catch_index_ML", "catch_bsurvey_ML") 
# avail_set_LC <- c("catch_LC", "catch_index_LC", "catch_bsurvey_LC")

da <- list("Nyears"=20, "Nyears_comp"=20, "comp_sample"=(1000/20)) ## site-specific operating model, 1000 samples total over 20 sites

## create combos
cmsy_modcombos <- as.matrix(expand.grid("Model"="CMSY", "Data_avail"=avail_set, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "LH"=paste0("LH_", lh_vec)))

## transform model combinations into directories
cmsy_dir_vec <- model_paths(modcombos=cmsy_modcombos, res_dir=results_dir)

## run iterations
itervec <- 1:10


##--------------------- setup parallel ----------------------
registerDoParallel(cores=5)

## ------------------ simulate data -------------------------

start_datagen <- Sys.time()

## create true population and generated data into directories
foreach(loop=1:length(cmsy_dir_vec), .packages=c('LIME','catchMSY')) %dopar% generateData(modpath=cmsy_dir_vec[loop], itervec=itervec, spatial=TRUE, Fdynamics=strsplit(cmsy_modcombos[loop,"Fdyn"],"_")[[1]][2], Rdynamics=strsplit(cmsy_modcombos[loop,"Rdyn"],"_")[[1]][2], LType=1, write=TRUE, lh_list=lh, data_avail_list=da, modname=paste0(cmsy_modcombos[loop,"Model"],"_",cmsy_modcombos[loop,"Data_avail"]), rewrite=TRUE)

end_datagen <- Sys.time() - start_datagen

## -------------- run estimation model -----------------------

## catchMSY
start_run <- Sys.time()

foreach(loop=1:length(cmsy_dir_vec), .packages=c('LIME', 'catchMSY')) %dopar% tryCatch(run_cmsy(modpath=cmsy_dir_vec[loop], itervec=itervec, lh_list=lh, data_avail=cmsy_modcombos[loop,"Data_avail"], nyears=20, rewrite=TRUE), error=function(e) print(paste0("issue with ", cmsy_dir_vec[loop])))

end_run <- Sys.time() - start_run


#test
loop <- 1
modpath=cmsy_dir_vec[loop]
itervec=1
lh_list=lh
data_avail=cmsy_modcombos[loop,"Data_avail"]
nyears=20
rewrite=FALSE

modpath=cmsy_dir_vec[loop]
itervec=1
spatial=TRUE
Fdynamics=strsplit(cmsy_modcombos[loop,"Fdyn"],"_")[[1]][2]
Rdynamics=strsplit(cmsy_modcombos[loop,"Rdyn"],"_")[[1]][2]
LType=1
write=TRUE
lh_list=lh
data_avail_list=da
modname=paste0(cmsy_modcombos[loop,"Model"],"_",cmsy_modcombos[loop,"Data_avail"])
rewrite=TRUE

generateData(modpath=cmsy_dir_vec[loop], itervec=1, spatial=TRUE, Fdynamics=strsplit(cmsy_modcombos[loop,"Fdyn"],"_")[[1]][2], Rdynamics=strsplit(cmsy_modcombos[loop,"Rdyn"],"_")[[1]][2], LType=1, write=TRUE, lh_list=lh, data_avail_list=da, modname=paste0(cmsy_modcombos[loop,"Model"],"_",cmsy_modcombos[loop,"Data_avail"]), rewrite=TRUE)

true <- readRDS(file.path(cmsy_dir_vec[loop], 1, "True.rds"))

run_cmsy(modpath=cmsy_dir_vec[loop], itervec=1, lh_list=lh, data_avail=cmsy_modcombos[loop,"Data_avail"], nyears=20, rewrite=TRUE)

stats <- readRDS(file.path(cmsy_dir_vec[loop], "stats.rds"))
sum_stats <- readRDS(file.path(cmsy_dir_vec[loop], "summary_stats.rds"))


# ## lbspr
# # devtools::install_github("AdrianHordyk/LBSPR", dependencies=TRUE)
# library(LBSPR)

# start_run_lbspr <- Sys.time()

# foreach(loop=1:length(lbspr_dirs), .packages=c('LBSPR')) %dopar% try(runLBSPR(modpath=lbspr_dirs[loop], lh_list=lh, itervec=itervec, rewrite=FALSE, simulation=TRUE))

# end_run_lbspr <- Sys.time() - start_run_lbspr

# compare_re(dir_vec=c(mixe_dirs), mod_names=c(paste0("LC", c(1,10))), Rdyn_vec=Rdyn_set, Fdyn_vec=Fdyn_set, lh_num=1:5, save=TRUE, type="all", fig_name="comp_only")

# SPRcover <- interval_coverage(modpath_vec=mixe_dirs, param="SPR", itervec=itervec)
# compare_fits(dir_vec=mixe_dirs, mod_names=paste0("LC", c(1,10)), Ftype="Ramp", Rtype="Pulsed", fig_name="Fits_LC", LHchoose=5, iter=7, save=TRUE, coverage=SPRcover$pcover, type="all")
