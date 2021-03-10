#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(methods)
library(TMB)
library(mgcv)
dyn.load(dynlib("../binom_gamm_nsmooth"))
dyn.load(dynlib("../betabin_gamm_nsmooth"))
source("../bootstrap.cl.less.r")
source("../fit_tmb.R")
load("combined_data.RData")
x = readRDS("model_fits.RDS")
set.seed(as.integer(args[1]))
n = as.integer(args[2])
boot.set = as.integer(args[3])
model = args[4]
nc = length(x$model.fits[[model]]$model.res$rep$mean_pred_eta)
plen = combined.data$predict.length
new.dat = data.frame(dn = factor(rep(c("day","night"), each = length(plen))), length = rep(plen,2))
do_pred = TRUE
boot.pred.eta = matrix(NA, n, nc)

for(i in 1:n)
{
  #print(paste0("i: ", i))
  dat = bootstrap.cl.less(combined.data$cl.less)
  dat$dn = dat$day.night
  out = try(eval(x$model.fits[[model]]$call))
  out$sdrep = try(summary(sdreport(out$model.res)))
  if(!is.character(out$sdrep)) boot.pred.eta[i,] = out$sdrep[which(rownames(out$sdrep) == "mean_pred_eta"),1]
  save(boot.pred.eta, file = paste0("boot_pred_eta_", boot.set, ".RData"))
  print(paste0("bootset: ", boot.set, ", i: ", i, " done"))  
}

