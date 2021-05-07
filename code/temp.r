#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

parentdir = getwd()
library(methods)
library(TMB)
library(mgcv)
dyn.load(dynlib(paste0(parentdir,"/code/binom_gamm_nsmooth")))
dyn.load(dynlib(paste0(parentdir,"/code/betabin_gamm_nsmooth")))
source(paste0(parentdir,"/code/estimate.binom.gamm.nsmooth.r"))
source(paste0(parentdir,"/code/estimate.betabin.gamm.nsmooth.r"))
source(paste0(parentdir,"/code/bootstrap.cl.less.r"))
source(paste0(parentdir,"/code/fit_tmb.R"))

sp.info = readRDS("data/sp.info.RDS")
args = c("023456", "100", "0","bb7","2")
args = c("023456", "3", "0","bb7","2")
set.seed(as.integer(args[1]))
n = as.integer(args[2])
boot.set = as.integer(args[3])
model = args[4]
sp.i = as.integer(args[5])
sp = sp.info$sp.names[sp.i]
x = readRDS(paste0(parentdir, "/results/big_results/", sp, "_model_fits.RDS"))
combined.data = readRDS(paste0(parentdir, "/data/", sp, "_data.RDS"))

nc = length(x$model.fits[[model]]$model.res$rep$mean_pred_eta)
plen = combined.data$predict.length
new.dat = data.frame(dn = factor(rep(c("day","night"), each = length(plen))), length = rep(plen,2))
do_pred = TRUE
boot.pred.eta = matrix(NA, n, nc)

for(i in 1:n)
{
  #print(paste0("i: ", i))
  dat = bootstrap.cl.less(combined.data$cl.less)
}
  dat$dn = dat$day.night
  out = try(eval(x$model.fits[[model]]$call))
  out$sdrep = try(summary(sdreport(out$model.res)))
  if(!is.character(out$sdrep)) boot.pred.eta[i,] = out$sdrep[which(rownames(out$sdrep) == "mean_pred_eta"),1]
  #saveRDS(boot.pred.eta, file = paste0(parentdir, "/results/", sp, "_boot_pred_eta_", boot.set, ".RDS"))
  #print(paste0("bootset: ", boot.set, ", i: ", i, " done"))  
}

