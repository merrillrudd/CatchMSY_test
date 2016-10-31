
catch_dirs_nv <- cmsy_dir_vec[which(grepl("catch/", cmsy_dir_vec) & grepl("SigmaR_0/", cmsy_dir_vec) & grepl("CRSNAP", cmsy_dir_vec))]
catch_dirs_v <- cmsy_dir_vec[which(grepl("catch/", cmsy_dir_vec) & grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("CRSNAP", cmsy_dir_vec))]

png(file.path(fig_dir, "FR_scenarios.png"), height=8, width=12, units="in", res=200)
par(mfrow=c(2,4), mar=c(0,0,0,0), omi=c(1,1,1,1))
## loop over 4 fishing mortality scenarios
for(i in 1:4){
  plot(x=1, y=1, type="n", xlim=c(1,20), ylim=c(0,1.5), xaxt="n", yaxt="n", xaxs="i", yaxs="i", cex.axis=2)
  if(i==1) axis(2, at=pretty(c(0,1)), las=2, cex.axis=2)
  if(i==4) axis(1, at=pretty(c(0,20)), cex.axis=2)
  if(i==1) mtext(side=2, "Fishing mortality", cex=1.5, line=4)

  for(ii in itervec){
      true_v <- readRDS(file.path(catch_dirs_v[i], ii, "True.rds"))
      lines(true_v$F_t, col="#AA000050", lwd=2, xaxt="n", yaxt="n", xaxs="i", yaxs="i")
      if(ii==length(itervec)){
        true_nv <- readRDS(file.path(catch_dirs_nv[i], ii, "True.rds"))
        lines(true_nv$F_t, col="black", lwd=4, xaxt="n", yaxt="n", xaxs="i", yaxs="i")
      }
  }
  if(i==1) mtext(side=3, "Endogenous", cex=1.5, line=-2)
  if(i==2) mtext(side=3, "Increasing", cex=1.5, line=-2)
  if(i==3) mtext(side=3, "Constant", cex=1.5, line=-2)
  if(i==4) mtext(side=3, "Ramped", cex=1.5, line=-2)
}
## loop over 3 recruitment scenarios
for(i in c(1,5,9)){
  plot(x=1, y=1, type="n", xlim=c(1,20), ylim=c(0,4), yaxt="n", xaxs="i", yaxs="i", cex.axis=2)
  if(i==1) axis(2, at=pretty(c(0,3)), las=2, cex.axis=2)
  if(i==1) mtext(side=2, "Relative recruitment", cex=1.5, line=4)
  for(ii in itervec){
      true_v <- readRDS(file.path(catch_dirs_v[i], ii, "True.rds"))
      lines(true_v$R_t/mean(true_v$R_t), col="#0000AA50", lwd=2, xaxt="n", yaxt="n", xaxs="i", yaxs="i")
      if(ii==length(itervec)){
        true_nv <- readRDS(file.path(catch_dirs_nv[i], ii, "True.rds"))
        lines(true_nv$R_t/mean(true_nv$R_t), col="black", lwd=4, xaxt="n", yaxt="n", xaxs="i", yaxs="i")
      }
  }
  if(i==1) mtext(side=3, "Constant", cex=1.5, line=-2)
  if(i==5) mtext(side=3, "Pulsed", cex=1.5, line=-2)
  if(i==9) mtext(side=3, "Beverton-Holt", cex=1.5, line=-2)
}
mtext(side=1, "Year", outer=TRUE, cex=1.5, line=3)
dev.off()


## -------------- figures  -----------------------
compare_re(dir_vec=cmsy_dir_vec[grepl("SigmaR_0/", cmsy_dir_vec)], mod_names=avail_set, Fdyn_vec=Fdyn_set, Rdyn_vec=Rdyn_set, SigmaR_vec=c(0), lh_num=lh_vec, save=TRUE, fig_name="RE_cmsy")
compare_re(dir_vec=cmsy_dir_vec[grepl("SigmaR_0.6/", cmsy_dir_vec)], mod_names=avail_set, Fdyn_vec=Fdyn_set, Rdyn_vec=Rdyn_set, SigmaR_vec=c(0.6), lh_num=lh_vec, save=TRUE, fig_name="RE_cmsy_RecVar")

dir <- cmsy_dir_vec[which(grepl("catch/", cmsy_dir_vec) & grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("HAKE", cmsy_dir_vec) & grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec))]
reout <- sapply(itervec, function(x) re_calc(dir=dir, iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=dir, iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fendog_RecVar.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

dir <- cmsy_dir_vec[which(grepl("index/", cmsy_dir_vec) & grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("HAKE", cmsy_dir_vec) & grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec))]
reout <- sapply(itervec, function(x) re_calc(dir=dir, iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=dir, iter=x))
png(file.path(fig_dir, "MSY_iters_index_Fendog_RecVar.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()



dir <- cmsy_dir_vec[which(grepl("catch/", cmsy_dir_vec) & grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("HAKE", cmsy_dir_vec) & grepl("F_Constant", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec))]
reout <- sapply(itervec, function(x) re_calc(dir=dir, iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=dir, iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fconstant_RecVar.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()


dir <- cmsy_dir_vec[which(grepl("index/", cmsy_dir_vec) & grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("HAKE", cmsy_dir_vec) & grepl("F_Constant", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec))]
reout <- sapply(itervec, function(x) re_calc(dir=dir, iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=dir, iter=x))
png(file.path(fig_dir, "MSY_iters_index_Fconstant_RecVar.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()


png(file.path(fig_dir, "MSY_iters_catchBD_Fendogenous_RecVar.png"), height=10, width=15, res=200, units="in")
par(mfrow=c(2,1), mar=c(4.5,5,1,1))
dir <- bdcmsy_dir_vec[which(grepl("catch/", bdcmsy_dir_vec) & grepl("SigmaR_0.6/", bdcmsy_dir_vec) & grepl("HAKE", bdcmsy_dir_vec) & grepl("F_Endogenous", bdcmsy_dir_vec) & grepl("R_Constant", bdcmsy_dir_vec))]
reoutbd <- sapply(itervec, function(x) re_calc_bd(bd=dir, iter=x))
covoutbd <- sapply(itervec, function(x) cover_calc_bd(bd=dir, iter=x))
plot(x=1,y=1,type="n",ylim=c(-1,10),xlim=c(0,max(itervec)+1), cex.lab=2, xlab="Iteration", ylab="Relative Error", xaxt="n")
axis(1, at=pretty(c(1,length(itervec))))
for(i in 1:ncol(reoutbd)){
  segments(x0=i, y0=as.numeric(reoutbd["relcl",i]), y1=as.numeric(reoutbd["reucl",i]), col="gray", lwd=4)
}
points(x=1:ncol(reoutbd), y=reoutbd["re",], col=covoutbd, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4), pch=19, cex=1.5)
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reoutbd)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covoutbd)=="gray"))/length(covoutbd),2)), xy=c(0.925,0.885), cex=1.5)

dir <- cmsy_dir_vec[which(grepl("catch/", cmsy_dir_vec) & grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("HAKE", cmsy_dir_vec) & grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec))]
reout <- sapply(itervec, function(x) re_calc(dir=dir, iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=dir, iter=x))
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, xlim=c(1,length(itervec)+1), ylim=c(-1,10), xaxt="n")
abline(h=0, lwd=4, col="red")
axis(1, at=pretty(c(1,length(itervec))))
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

png(file.path(fig_dir, "MSY_iters_indexBD_Fendogenous_RecVar.png"), height=10, width=15, res=200, units="in")
par(mfrow=c(2,1), mar=c(4.5,5,1,1))
dir <- bdcmsy_dir_vec[which(grepl("index/", bdcmsy_dir_vec) & grepl("SigmaR_0.6/", bdcmsy_dir_vec) & grepl("HAKE", bdcmsy_dir_vec) & grepl("F_Endogenous", bdcmsy_dir_vec) & grepl("R_Constant", bdcmsy_dir_vec))]
reoutbd <- sapply(itervec, function(x) re_calc_bd(bd=dir, iter=x))
covoutbd <- sapply(itervec, function(x) cover_calc_bd(bd=dir, iter=x))
plot(x=1,y=1,type="n",ylim=c(-1,2.5),xlim=c(0,max(itervec)+1), cex.lab=2, xlab="Iteration", ylab="Relative Error", xaxt="n")
axis(1, at=pretty(c(1,length(itervec))))
for(i in 1:ncol(reoutbd)){
  segments(x0=i, y0=as.numeric(reoutbd["relcl",i]), y1=as.numeric(reoutbd["reucl",i]), col="gray", lwd=4)
}
points(x=1:ncol(reoutbd), y=reoutbd["re",], col=covoutbd, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4), pch=19, cex=1.5)
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reoutbd)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covoutbd)=="gray"))/length(covoutbd),2)), xy=c(0.925,0.885), cex=1.5)

dir <- cmsy_dir_vec[which(grepl("index/", cmsy_dir_vec) & grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("HAKE", cmsy_dir_vec) & grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec))]
reout <- sapply(itervec, function(x) re_calc(dir=dir, iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=dir, iter=x))
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-1,2.5), xaxt="n")
axis(1, at=pretty(c(1,length(itervec))))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()








plot(x=1:20, y=true$F_t, type="l", lwd=4, ylim=c(0,1))
med <- sapply(1:ncol(out$ps.ft), function(x) median(out$ps.ft[out$idx,x], na.rm=TRUE))
ucl <- sapply(1:ncol(out$ps.ft), function(x) quantile(out$ps.ft[out$idx,x], 0.975, na.rm=TRUE))
lcl <- sapply(1:ncol(out$ps.ft), function(x) quantile(out$ps.ft[out$idx,x], 0.025, na.rm=TRUE))
polygon(x=c(1:20, 20:1), y=c(lcl, rev(ucl)), col="#AAAAAA50", border=NA)
lines(x=1:20, y=med, col="red", lwd=4)




model <- 2
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_index_Fendog_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 3
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_ml_Fendog_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()


model <- 4
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fconstant_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 6
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_ml_Fconstant_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 28
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fendog_sigR0.6.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 31
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fconstant_sigR0.6.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- 33
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_ml_Fconstant_sigR0.6.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()


model <- 19
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fendog_Rpulse_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- which(grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec) & grepl("CRSNAP", cmsy_dir_vec) & grepl("SigmaR_0/", cmsy_dir_vec) & grepl("catch/", cmsy_dir_vec))
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fendog_Rconstant_sigR0_SNAP.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- which(grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec) & grepl("SIGSUT", cmsy_dir_vec) & grepl("SigmaR_0/", cmsy_dir_vec) & grepl("catch/", cmsy_dir_vec))
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fendog_Rconstant_sigR0_SIG.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- which(grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec) & grepl("SIGSUT", cmsy_dir_vec) & grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("catch/", cmsy_dir_vec))
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_catch_Fendog_Rconstant_sigR0.6_SIG.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- which(grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Pulsed", cmsy_dir_vec) & grepl("HAKE", cmsy_dir_vec) & grepl("SigmaR_0/", cmsy_dir_vec) & grepl("catch_index/", cmsy_dir_vec))
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_index_Fendog_Rpulse_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()


model <- which(grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec) & grepl("CRSNAP", cmsy_dir_vec) & grepl("SigmaR_0/", cmsy_dir_vec) & grepl("catch_index/", cmsy_dir_vec))
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_index_Fendog_Rconstant_sigR0_SNAP.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- which(grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec) & grepl("SIGSUT", cmsy_dir_vec) & grepl("SigmaR_0/", cmsy_dir_vec) & grepl("catch_index/", cmsy_dir_vec))
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_index_Fendog_Rconstant_sigR0_SIG.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- which(grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Pulsed", cmsy_dir_vec) & grepl("HAKE", cmsy_dir_vec) & grepl("SigmaR_0/", cmsy_dir_vec) & grepl("catch_ML/", cmsy_dir_vec))
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_ml_Fendog_Rpulsed_sigR0.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- which(grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec) & grepl("CRSNAP", cmsy_dir_vec) & grepl("SigmaR_0/", cmsy_dir_vec) & grepl("catch_ML/", cmsy_dir_vec))
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_ml_Fendog_Rconstant_sigR0_SNAP.png"), height=10, width=15, res=200, units="in")
boxplot(reout, col=covout, xlab="Iteration", ylab="Relative Error", cex.lab=2, ylim=c(-2,4))
abline(h=0, lwd=4, col="red")
print.letter(paste0("RE = ", round(median(unlist(reout)),2)), xy=c(0.925,0.925), cex=1.5)
print.letter(paste0("Cover = ", round(length(which(unlist(covout)=="gray"))/length(covout),2)), xy=c(0.925,0.885), cex=1.5)
dev.off()

model <- which(grepl("F_Endogenous", cmsy_dir_vec) & grepl("R_Constant", cmsy_dir_vec) & grepl("SIGSUT", cmsy_dir_vec) & grepl("SigmaR_0/", cmsy_dir_vec) & grepl("catch_ML/", cmsy_dir_vec))
reout <- sapply(itervec, function(x) re_calc(dir=cmsy_dir_vec[model], iter=x))
covout <- sapply(itervec, function(x) cover_calc(dir=cmsy_dir_vec[model], iter=x))
png(file.path(fig_dir, "MSY_iters_ml_Fendog_Rconstant_sigR0_SIG.png"), height=10, width=15, res=200, units="in")
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



compare_re(dir_vec=cmsy_dir_vec[which(grepl("SigmaR_0.6/", cmsy_dir_vec) & grepl("ML", cmsy_dir_vec)==FALSE)], mod_names=avail_set[-which(grepl("ML",avail_set))], Fdyn_vec="Endogenous", Rdyn_vec=Rdyn_set, SigmaR_vec=c(0.6), lh_num=lh_vec, save=TRUE, fig_name="RE_compare_agebio_Fendog_Sigma0.6", bd_dir_vec=bdcmsy_dir_vec[grepl("SigmaR_0.6/",bdcmsy_dir_vec)], itervec=1:20, ylim=c(-2.5,6))




#### no other data types included
### increasing fishing mortality, zero recruitment variation
model=1
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
png(file.path(fig_dir, "resid_catch_Fendog_SigmaR0.png"), width=10, height=8, res=200, units="in")
resid_plot(out=out, true=true, ylim=c(-4,4))
dev.off()

model=1
  out <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
  true <- readRDS(file.path(cmsy_dir_vec[model], 1, "True.rds"))
par(mfrow=c(1,1), mar=c(5,5,3,3), omi=c(0.2,0.2,0.2,0.2))
hist(out$S[,"msy"], col="gray", xlim=c(0,800), ylim=c(0,600), xaxs="i", yaxs="i", border=NA, main="", xlab="MSY", cex.lab=2, cex.axis=2)
par(new=TRUE)
hist(out$S[out$idx,"msy"], col="red", xlim=c(0,800),ylim=c(0,600), border=NA, main="", xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i", )

model <- 2
  out2 <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
par(new=TRUE)
hist(out2$S[out2$idx,"msy"], col="blue", xlim=c(0,800),ylim=c(0,600), border=NA, main="", xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i", )

model <- 4
  out3 <- readRDS(file.path(cmsy_dir_vec[model], 1, "cmsy_output.rds"))
par(new=TRUE)
hist(out3$S[out3$idx,"msy"], col="yellow", xlim=c(0,800),ylim=c(0,600), border=NA, main="", xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i", )



png(file.path(fig_dir, "example_datagen_SigmaR0.6.png"), height=10, width=8, res=200, units="in")
par(mfrow=c(3,2), mar=c(0,0,0,0), omi=c(1,1.5,1,1))
  # true1 <- readRDS(file.path(cmsy_dir_vec[1], 1, "True.rds"))
  # true2 <- readRDS(file.path(cmsy_dir_vec[7], 1, "True.rds"))
  # true3 <- readRDS(file.path(cmsy_dir_vec[19], 1, "True.rds"))
  # true4 <- readRDS(file.path(cmsy_dir_vec[25], 1, "True.rds"))
  true1 <- readRDS(file.path(cmsy_dir_vec[28], 1, "True.rds"))
  true2 <- readRDS(file.path(cmsy_dir_vec[31], 1, "True.rds"))
  true3 <- readRDS(file.path(cmsy_dir_vec[37], 1, "True.rds"))
  true4 <- readRDS(file.path(cmsy_dir_vec[41], 1, "True.rds"))
plot(true1$F_t, ylim=c(0,0.5), type="l", lwd=4, xaxt="n", yaxt="n")
lines(true2$F_t, lwd=4, col="red", lty=2)
axis(2, cex.axis=2, las=2)
mtext(side=2, "Fishing\nmortality", cex=2, line=4)

plot(true3$F_t, ylim=c(0,0.5), type="l", lwd=4, xaxt="n", yaxt="n")
lines(true4$F_t, lwd=4, col="red", lty=2)

plot(true1$R_t, ylim=c(0,2000), type="l", lwd=4, xaxt="n", yaxt="n")
lines(true2$R_t, lwd=4, col="red", lty=2)
axis(2, cex.axis=2, las=2)
mtext(side=2, "Recruitment", cex=2, line=6)

plot(true3$R_t, ylim=c(0,2000), type="l", lwd=4, xaxt="n", yaxt="n")
lines(true4$R_t, lwd=4, col="red", lty=2)

plot(true1$D_t, ylim=c(0,1), type="l", lwd=4, cex.axis=2, yaxt="n")
lines(true2$D_t, lwd=4, col="red", lty=2)
axis(2, cex.axis=2, las=2)
mtext(side=2, "Depletion\n(B/K)", cex=2, line=4)

plot(true3$D_t, ylim=c(0,1), type="l", lwd=4, cex.axis=2, yaxt="n")
lines(true4$D_t, lwd=4, col="red", lty=2)

# plot(true1$ML_t, ylim=c(0,100), type="l", lwd=4, cex.axis=2, yaxt="n")
# lines(true2$ML_t, lwd=4, col="red", lty=2)
# axis(2, cex.axis=2, las=2)
# mtext(side=2, "Mean length", cex=2, line=4)

# plot(true3$ML_t, ylim=c(0,100), type="l", lwd=4, cex.axis=2, yaxt="n")
# lines(true4$ML_t, lwd=4, col="red", lty=2)

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



### scatterplot with priors
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

## scatterplot with priors and posteriors
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


## depletion trajectories
## catch only
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

## catch + index
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



