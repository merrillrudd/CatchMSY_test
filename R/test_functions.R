run_cmsy <- function(modpath, itervec, lh_list, data_avail, nyears, selex=FALSE, rewrite=FALSE, nsamp=5000, ncores=1){
	
	lh_num <- ifelse(grepl("LH1", modpath), 1, ifelse(grepl("LH2", modpath), 2, ifelse(grepl("LH3", modpath), 3, ifelse(grepl("LH4", modpath), 4, ifelse(grepl("LH5", modpath), 5, ifelse(grepl("CRSNAP", modpath), "CRSNAP", ifelse(grepl("SIGSUT", modpath), "SIGSUT", ifelse(grepl("HAKE", modpath), "HAKE", stop("No match to life history number")))))))))
  	lh_choose <- lh_list[[lh_num]]

  	species <- new_sID(id=lh_num)
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

  	species$sel50 <- lh_choose$S50
  	species$sel95 <- species$sel50+1

  	stats <- list("msy"=matrix(NA, nrow=length(itervec), ncol=3), "fmsy"=matrix(NA, nrow=length(itervec), ncol=3), "m"=matrix(NA, nrow=length(itervec), ncol=3))
    if(selex==TRUE) stats$sel50 <- matrix(NA, nrow=length(itervec), ncol=3)
  		colnames(stats$msy) <- colnames(stats$fmsy) <- colnames(stats$m) <- c("re", "devs", "cover")
      if(selex==TRUE) colnames(stats$sel50) <- c("re", "devs", "cover")

    if(rewrite==TRUE & file.exists(file.path(modpath, "stats.rds"))) unlink(file.path(modpath, "stats.rds"), TRUE)
    if(rewrite==TRUE & file.exists(file.path(modpath, "summary_stats.rds"))) unlink(file.path(modpath, "summary_stats.rds"), TRUE)

  for(iter in 1:length(itervec)){

      data_gen <- readRDS(file.path(modpath, itervec[iter], "True.rds"))
      true_msy <- data_gen$MSY
      true_fmsy <- data_gen$Fmsy
      true_m <- data_gen$M
      true_sel50 <- lh_choose$S50

  		if(rewrite==FALSE & file.exists(file.path(modpath, itervec[iter], "cmsy_output.rds"))){
        species <- readRDS(file.path(modpath, itervec[iter], "cmsy_output.rds"))

        stats$msy[iter,"re"] <- (species$msy.stats["Median"] - true_msy)/true_msy
        stats$fmsy[iter,"re"] <- (species$fmsy.stats["Median"] - true_fmsy)/true_fmsy
        stats$m[iter,"re"] <- (species$m.stats["Median"] - true_m)/true_m
        if(selex==TRUE) stats$sel50[iter,"re"] <- (species$sel50.stats["Median"] - true_sel50)/true_sel50   

        stats$msy[iter,"devs"] <- (species$msy.stats["Median"] - true_msy)^2
        stats$fmsy[iter,"devs"] <- (species$fmsy.stats["Median"] - true_fmsy)^2
        stats$m[iter,"devs"] <- (species$m.stats["Median"] - true_m)^2
        if(selex==TRUE) stats$sel50[iter,"devs"] <- (species$sel50.stats["Median"] - true_sel50)^2    

        stats$msy[iter,"cover"] <- ifelse(species$msy.stats["Min."]<= true_msy & species$msy.stats["Max."]>= true_msy, 1, 0)
        stats$fmsy[iter,"cover"] <- ifelse(species$fmsy.stats["Min."]<= true_fmsy & species$fmsy.stats["Max."]>= true_fmsy, 1, 0)
        stats$m[iter,"cover"] <- ifelse(species$m.stats["Min."]<= true_m & species$m.stats["Max."]>= true_m, 1, 0)
        if(selex==TRUE) stats$sel50[iter,"cover"] <- ifelse(species$sel50.stats["Min."]<= true_sel50 & species$sel50.stats["Max."]>=true_sel50, 1, 0)
        next
      }

      if(rewrite==TRUE | file.exists(file.path(modpath, itervec[iter], "cmsy_output.rds"))==FALSE){

      		catch <- data_gen$C_t
      		index <- data_gen$I_t
      		index.lse <- rep(0.2, length(index))
      		biomass <- index/lh_choose$qcoef
      		biomass.lse <- rep(0.2, length(catch))
      		lc <- data_gen$LF
      			bins <- lh_choose$mids
      			colnames(lc) <- paste0("lc.",bins)
      		ml <- data_gen$ML_t
      		ml.lse <- rep(0.6, length(ml))    

      		data_input <- data.frame("year"=1:nyears)
      		if(grepl("catch", data_avail)) data_input$catch <- catch
      		if(grepl("index", data_avail)){
      			data_input$index <- index
      			data_input$index.lse <- index.lse
      		}
      		if(grepl("bsurvey", data_avail)){
      			data_input$biomass <- biomass
      			data_input$biomass.lse <- biomass.lse
      		}
      		if(grepl("LC", data_avail)){
            data_input <- cbind(data_input, lc)
            data_input$ess <- 1
          }
      		if(grepl("ML", data_avail)){
      			data_input$meanlength <- ml
      			data_input$meanlength.lse <- ml.lse
      		}   

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
        selexPriorInfo <- data.frame("id"=4, "dist"="lnorm", "par1"=log(species$sel50), "par2"=0.1*species$sel50, "log"=TRUE, "stringAsFactors"=FALSE)
        if(nrow(species$dfPriorInfo)==3) species$dfPriorInfo <- rbind.data.frame(species$dfPriorInfo, selexPriorInfo)
        if(nrow(species$dfPriorInfo)==4){
          species$dfPriorInfo$dist[4] <- selexPriorInfo$dist
          species$dfPriorInfo$par1[4] <- selexPriorInfo$par1
          species$dfPriorInfo$par2[4] <- selexPriorInfo$par2
        }   

    		species <- sample.sid(sID=species, selex=selex, n=nsamp)
         		# species$m <- species$S[1,1]
         		# species$fmsy <- species$S[1,2]
         		# species$msy <- species$S[1,3]
    		species_new <- tryCatch(sir.sid(sID=species, selex=selex, ncores=ncores), error=function(e) NA)
        if(all(is.na(species_new))) write("error in model run", file.path(modpath, itervec[iter], "error.txt"))
        if(all(is.na(species_new))==FALSE){
          species <- species_new
          species$msy.stats <- summary(species$S[species$idx,3])
          species$fmsy.stats <- summary(species$S[species$idx,2])
          species$m.stats <- summary(species$S[species$idx,1])
          if(selex==TRUE) species$sel50.stats <- summary(species$S[species$idx,4])         
          saveRDS(species, file.path(modpath, itervec[iter], "cmsy_output.rds"))
        }
      }

    if(all(is.na(species))==FALSE){
      stats$msy[iter,"re"] <- (species$msy.stats["Median"] - true_msy)/true_msy
      stats$fmsy[iter,"re"] <- (species$fmsy.stats["Median"] - true_fmsy)/true_fmsy
      stats$m[iter,"re"] <- (species$m.stats["Median"] - true_m)/true_m
      if(selex==TRUE) stats$sel50[iter,"re"] <- (species$sel50.stats["Median"] - true_sel50)/true_sel50 

      stats$msy[iter,"devs"] <- (species$msy.stats["Median"] - true_msy)^2
      stats$fmsy[iter,"devs"] <- (species$fmsy.stats["Median"] - true_fmsy)^2
      stats$m[iter,"devs"] <- (species$m.stats["Median"] - true_m)^2
      if(selex==TRUE) stats$sel50[iter,"devs"] <- (species$sel50.stats["Median"] - true_sel50)^2  

      stats$msy[iter,"cover"] <- ifelse(species$msy.stats["Min."]<= true_msy & species$msy.stats["Max."]>= true_msy, 1, 0)
      stats$fmsy[iter,"cover"] <- ifelse(species$fmsy.stats["Min."]<= true_fmsy & species$fmsy.stats["Max."]>= true_fmsy, 1, 0)
      stats$m[iter,"cover"] <- ifelse(species$m.stats["Min."]<= true_m & species$m.stats["Max."]>= true_m, 1, 0)
      if(selex==TRUE) stats$sel50[iter,"cover"] <- ifelse(species$sel50.stats["Min."]<= true_sel50 & species$sel50.stats["Max."]>=true_sel50, 1, 0)      
    }
    if(all(is.na(species))) next
  }

  saveRDS(stats, file.path(modpath, "stats.rds"))
	summary_stats <- list("msy"=NULL, "fmsy"=NULL, "m"=NULL)
  if(selex==TRUE) summary_stats$sel50 <- NULL
	summary_stats$msy$mre <- median(stats$msy[,"re"])
	summary_stats$msy$rmse <- sqrt(sum(stats$msy[,"devs"])/length(stats$msy[,"devs"]))
	summary_stats$msy$pcover <- sum(stats$msy[,"cover"])/length(stats$msy[,"cover"])
	summary_stats$fmsy$mre <- median(stats$fmsy[,"re"])
	summary_stats$fmsy$rmse <- sqrt(sum(stats$fmsy[,"devs"])/length(stats$fmsy[,"devs"]))
	summary_stats$fmsy$pcover <- sum(stats$fmsy[,"cover"])/length(stats$fmsy[,"cover"])
	summary_stats$m$mre <- median(stats$m[,"re"])
	summary_stats$m$rmse <- sqrt(sum(stats$m[,"devs"])/length(stats$m[,"devs"]))
	summary_stats$m$pcover <- sum(stats$m[,"cover"])/length(stats$m[,"cover"])
  if(selex==TRUE){
      summary_stats$sel50$mre <- median(stats$sel50[,"re"])
      summary_stats$sel50$rmse <- sqrt(sum(stats$sel50[,"devs"])/length(stats$sel50[,"devs"]))
      summary_stats$sel50$pcover <- sum(stats$sel50[,"cover"])/length(stats$sel50[,"cover"])
  }
	saveRDS(summary_stats, file.path(modpath, "summary_stats.rds"))

	return(paste0("ran ", length(itervec), " iters in ", modpath))
}

compare_re <- function(dir_vec, mod_names, Fdyn_vec, Rdyn_vec, lh_num, save, fig_name, itervec){
  for(ll in 1:length(lh_num)){
    if(save==TRUE){
      write_fig <- paste0(fig_name, "_LH", lh_num[ll], ".png")
      png(file.path(fig_dir, write_fig), res=200, units="in", width=29, height=14)
    }

      if(length(Fdyn_vec)==3 & length(Rdyn_vec)==3) par(mfcol=c(3,3), mar=c(0,0,0,0), omi=c(1,1.5,1,1), mgp=c(4,1,0))
      if(length(Fdyn_vec)==1 & length(Rdyn_vec)==1) par(mfcol=c(1,1), mar=c(0,0,0,0), omi=c(1,1.5,1,1), mgp=c(4,1,0))

      if(grepl("F_", dir_vec[1]) & grepl("R_", dir_vec[1])){
      for(ff in 1:length(Fdyn_vec)){
        for(rr in 1:length(Rdyn_vec)){
          dirs <- dir_vec[which(grepl(paste0("F_", Fdyn_vec[ff]), dir_vec) & grepl(paste0("R_", Rdyn_vec[rr]), dir_vec))]
          subdirs <- dirs[which(sapply(1:length(dirs), function(x) grepl(paste0("LH_",lh_num[ll]), dirs[x])))]
          stats <- sapply(1:length(subdirs), function(x) readRDS(file.path(subdirs[x], "stats.rds"))$msy[,"re"])
          sumstats <- sapply(1:length(subdirs), function(x) readRDS(file.path(subdirs[x], "summary_stats.rds"))$msy)
          boxplot(stats, ylim=c(-3, 3), col="turquoise", xaxt="n", yaxt="n", lwd=2)
          abline(h=0, lwd=5, lty=2)
          for(i in 1:ncol(sumstats)){
            text(x=i, y=-2.5, paste0("(", round(as.numeric(sumstats["mre",i]),2), ")\n", round(as.numeric(sumstats["pcover",i]),2)), cex=3)
          }
        if(ff==1){
          if(length(Rdyn_vec)>1) mtext(paste0( Rdyn_vec[rr]), side=2, line=5, font=2, cex=1.5)
          axis(side=2, las=2, cex.axis=1.2)
        }
        if(rr==1) if(length(Rdyn_vec)>1) mtext(paste0(Fdyn_vec[ff]), side=3, line=1, font=2, cex=1.5)
          mod_names_plot <- gsub("_", "\n", mod_names)
        if(rr==length(Rdyn_vec)) axis(side=1, at=1:length(mod_names), mod_names_plot, cex.axis=1.3, font=2)
        }
      }
    }
    if(grepl("F_", dir_vec[1])==FALSE & grepl("R_", dir_vec[1])==FALSE){
      dirs <- dir_vec
      subdirs <- dirs[which(sapply(1:length(dirs), function(x) grepl(paste0("LH_",lh_num[ll]), dirs[x])))]
          stats <- sapply(1:length(subdirs), function(x) readRDS(file.path(subdirs[x], "stats.rds"))$msy[,"re"])
          sumstats <- sapply(1:length(subdirs), function(x) readRDS(file.path(subdirs[x], "summary_stats.rds"))$msy)
          boxplot(stats, ylim=c(-3, 3), col="turquoise", xaxt="n", yaxt="n", lwd=2)
          abline(h=0, lwd=5, lty=2)
          for(i in 1:ncol(sumstats)){
            text(x=i, y=-2.5, paste0("(", round(as.numeric(sumstats["mre",i]),2), ")\n", round(as.numeric(sumstats["pcover",i]),2)), cex=3)
          } 
    }
    mtext("Model", side=1, cex=1.5, outer=TRUE, line=4)
    mtext("Relative error", side=2, cex=1.5, outer=TRUE, line=3)
    mtext("Fishing mortality pattern", side=3, line=3, cex=2, font=2, outer=TRUE)
    mtext("Recruitment pattern", side=2, line=7, cex=2, font=2, outer=TRUE)
    if(save==TRUE) dev.off()
  }
}

