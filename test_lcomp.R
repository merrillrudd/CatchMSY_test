rm(list=ls())

##############################################################
##### ------------- install packages ------------------- #####
##############################################################
## use operating model and life histories from LIME package
devtools::install_github("merrillrudd/LIME", build.vignettes=TRUE, dependencies=TRUE)
devtools::install_github("merrillrudd/catchMSY", build.vignettes=TRUE, dependencies=TRUE)
library(catchMSY)
library(LIME)

##############################################################
##### --------------- directories -----------------------#####
##############################################################

main_dir <- "C:\\Git_Projects\\CatchMSY_test"

## setup results directory
test_dir <- file.path(main_dir, "test_lcomp")
## remove results (for testing)
# unlink(test_dir, TRUE)
## create directory
dir.create(test_dir, showWarnings=FALSE)


##############################################################
##### --------- simulation-testing    ------------------ #####
##############################################################

##### --------------- testing options ------------------ #####

## test across different life histories 
lh_vec <- c("CRSNAP")#, "HAKE", "SIGSUT")
lh <- lapply(lh_vec, function(x) choose_lh_list(species=x, selex="asymptotic", param_adjust=c("R0"), val=c(1000), start_ages=1))
names(lh) <- lh_vec

## test across different true fishing mortality patterns 
Fdyn_set <- c("Endogenous")#, "Increasing", "Constant", "Ramp")

## test across different true recruitment patterns
Rdyn_set <- c("Constant")#, "Pulsed", "BH")

## test across different true levels of recruitment variation
SigmaR_set <- c(0)#, 0.6)

################################################
## Run cmsy models
################################################

## data availability scenarios
avail_set_LC <- c("catch_LC")#, "catch_index_LC", "catch_bsurvey_LC")

da <- list("Nyears"=20, "Nyears_comp"=20, "comp_sample"=1000) 

## create combos
lc_modcombos <- as.matrix(expand.grid("Model"="CMSY", "Data_avail"=avail_set_LC, "Fdyn"=paste0("F_",Fdyn_set), "Rdyn"=paste0("R_",Rdyn_set), "SigmaR"=paste0("SigmaR_",SigmaR_set), "LH"=paste0("LH_", lh_vec)))

## transform model combinations into directories
lc_dir_vec <- model_paths(modcombos=lc_modcombos, res_dir=test_dir)

itervec <- 1

loop <- 1
modpath <- lc_dir_vec[loop]
generateData(modpath=lc_dir_vec[loop], itervec=itervec, spatial=TRUE, Fdynamics=strsplit(lc_modcombos[loop,"Fdyn"],"_")[[1]][2], Rdynamics=strsplit(lc_modcombos[loop,"Rdyn"],"_")[[1]][2], LType=1, write=TRUE, lh_list=lh, data_avail_list=da, modname=paste0(lc_modcombos[loop,"Model"],"_",lc_modcombos[loop,"Data_avail"]), rewrite=TRUE, param_adjust=c("SigmaR","SigmaF"), val=c(as.numeric(strsplit(lc_modcombos[loop,"SigmaR"],"_")[[1]][2]), ifelse(as.numeric(strsplit(lc_modcombos[loop,"SigmaR"],"_")[[1]][2])==0,0,0.3)))

lh_name <- lh_vec[1]
lh_choose <- lh[[lh_name]]
species <- new_sID(id=lh_name)
species$linf <- lh_choose$linf
species$vbk <- lh_choose$vbk
species$to <- lh_choose$t0
species$a <- lh_choose$lwa
species$b <- lh_choose$lwb
species$winf <- species$a * species$linf^species$b
species$age <- 1:lh_choose$AgeMax
species$binwidth <- lh_choose$binwidth
species$m <- lh_choose$M
# species$m <- 1.50 * lh_choose$vbk
species$ah <- 1.65 / species$m
species$gh <- 0.10 * species$ah
species$sel1 <- lh_choose$S50
species$sel2 <- species$sel1+1

data_gen <- readRDS(file.path(modpath, loop, "True.rds"))
true_msy <- data_gen$MSY
true_fmsy <- data_gen$Fmsy
true_m <- data_gen$M
true_sel1 <- lh_choose$S50

nyears  <- data_gen$Nyears
data_input <- data.frame("year"=1:nyears)
catch <- data_gen$C_t
lc <- data_gen$LF
	bins <- lh_choose$mids
	colnames(lc) <- paste0("lc.",bins)
lencomp_sd <- rep(0.2, nrow(lc))
data_input$catch <- catch
data_input$lencomp_sd <- lencomp_sd
data_input <- cbind(data_input, lc)
  
species$data <- data_input  

# Set parameter sampling frame
species$dfPriorInfo$dist[1] = "lnorm"
species$dfPriorInfo$par1[1] = log(species$m)
species$dfPriorInfo$par2[1] = 0.05 * species$m
species$dfPriorInfo$dist[2] = "unif"
species$dfPriorInfo$par1[2] = 0.20 * species$m
species$dfPriorInfo$par2[2] = 1.50 * species$m
species$dfPriorInfo$dist[3] = "unif"
species$dfPriorInfo$par1[3] = quantile(species$data$catch,0.05)
species$dfPriorInfo$par2[3] = quantile(species$data$catch,0.95)
## age at 50% and 95% selectivity
selexPriorInfo <- data.frame("id"=4, "dist"="lnorm", "par1"=log(species$sel1), "par2"=0.1*species$sel1, "log"=TRUE, "stringAsFactors"=FALSE)
if(nrow(species$dfPriorInfo)==3) species$dfPriorInfo <- rbind.data.frame(species$dfPriorInfo, selexPriorInfo)
if(nrow(species$dfPriorInfo)==4){
  species$dfPriorInfo$dist[4] <- selexPriorInfo$dist
  species$dfPriorInfo$par1[4] <- selexPriorInfo$par1
  species$dfPriorInfo$par2[4] <- selexPriorInfo$par2
}   

species <- sample.sid(sID=species, selex=TRUE, n=10)

## testing function catchMSYModel
# sID <- species
# sID$m <- species$S[2,1]
# sID$fmsy <- species$S[2,2]
# sID$msy <- species$S[2,3]
# sID$sel1 <- species$S[2,4]
# sID$sel2 <- species$S[2,4] + 1
# attach(sID)

## testing function sir.sid
selex <- TRUE
sID <- species
S <- sID$S
fn <- function(s){
	sID$m    <- s[1]
	sID$fmsy <- s[2]
	sID$msy  <- s[3]
	if(selex==TRUE & sID$smodel=="logistic"){
		sID$sel1 <- s[4]
		sID$sel2 <- s[4]+1
	}
	if(selex==TRUE & sID$smodel=="dome") warning("Not programmed to estimate dome-shaped selectivity parameters")
	return(catchMSYModel(sID))
}
cmsy  <- apply(S,1,fn)		
		sID$code   <- plyr::ldply(cmsy,function(x){c("code"=x[['code']])})
		sID$bo     <- plyr::ldply(cmsy,function(x){c("bo"=x[['bo']])})
		sID$h      <- plyr::ldply(cmsy,function(x){c("h"=x[['h']])})
		sID$nll    <- plyr::ldply(cmsy,function(x){c("nll"=x[['nll']])})
		sID$prior  <- plyr::ldply(cmsy,function(x){c("prior"=x[['prior']])})
		sID$ps.bt  <- plyr::ldply(cmsy,function(x){c("bt"=x[['bt']])})
		sID$ps.dt  <- plyr::ldply(cmsy,function(x){c("dt"=x[['dt']])})
		sID$ps.sbt <- plyr::ldply(cmsy,function(x){c("sbt"=x[['sbt']])})
		sID$ps.ft  <- plyr::ldply(cmsy,function(x){c("ft"=x[['ft']])}) 
		sID$wts    <- exp(-(sID$nll + sID$prior))
		sID$ML <- plyr::ldply(cmsy, function(x){c("ML"=x[['ML']])})
		sID$LF <- sapply(1:length(cmsy), function(x) cmsy[[x]][['LF']])
		sID$spr_msy <- plyr::ldply(cmsy, function(x){c("spr_msy"=x[['spr_msy']])})
		sID$spr_t <- plyr::ldply(cmsy, function(x){c("spr_t"=x[['spr_t']])})
		sID$biomass_resid <- plyr::ldply(cmsy, function(x){c("biomass_resid"=x[['biomass_resid']])})
		sID$index_resid <- plyr::ldply(cmsy, function(x){c("index_resid"=x[['index_resid']])})
		sID$lc_resid <- sapply(1:length(cmsy), function(x) cmsy[[x]][['lc_resid']])
		sID$ml_resid <- plyr::ldply(cmsy, function(x){c("ml_resid"=x[['ml_resid']])})

		## non-statistical criterion - sample combinations that meet criterion
		sID$idx    <- which(sID$code==0)
		## statistical criterion - sample combinations that nll!=0
		## choose samples with the highest probability - best likelihoods
		ic <- which(sID$nll!=0)
		if( length(ic) > 0 ){
			prb <- sID$wts[ic,1]
			sID$idx  <- tryCatch(sample(ic,length(ic),replace=TRUE,prob=prb), error=function(e) NULL)
		}

				
## plot lc residuals
lc_resid <- t(sID$lc_resid[[3]]) ## choose iteration where lc_resid!=NULL
plot(x=1, y=1, type="n", xlim=c(0,21), ylim=c(0,98), xaxs="i", yaxs="i")
for(i in 1:nrow(lc_resid)){
	for(j in 1:ncol(lc_resid)){
		points(x=i,y=j, pch=19, cex=abs(as.numeric(lc_resid[i,j]))/mean(abs(as.numeric(lc_resid))), col=ifelse(lc_resid[i,j]>0,"#AA000050","#0000AA50"))
	}
}


## run function
species_new <- sir.sid(sID=species, selex=TRUE, ncores=1)
