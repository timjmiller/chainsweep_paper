#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(methods)
library(TMB)
library(mgcv)
dyn.load(dynlib("../binom_gamm"))
dyn.load(dynlib("../betabin_gamm"))
load("combined.dn.data.for.boot.RData")
set.seed(as.integer(args[1]))
n = as.integer(args[2])
#do.day = as.logical(as.integer(args[3]))
#do.night = as.logical(as.integer(args[4]))
#boot.set = as.integer(args[5])
#model = args[6]
boot.set = as.integer(args[3])
model = args[4]

cl.less = combined.dn.data$cl.less
cl.aug = combined.dn.data$cl.aug
pred.length = combined.dn.data$predict.length

boot.pred.eta = matrix(NA, n, length(pred.length))

for(i in 1:n)
{
  print(paste0("i: ", i))
  out = final.estimate.efficiency.fn(cl.less, cl.aug, fit.models = model, bootstrap = TRUE,
    predict.length = pred.length)
  out$sdrep = try(summary(sdreport(out$model.fits[[model]]$model.res)))
  if(!is.character(out$sdrep)) boot.pred.eta[i,] = out$sdrep[which(rownames(out$sdrep) == "mean_pred_eta"),1]
  save(boot.pred.eta, file = paste0("boot_pred_eta_", boot.set, ".RData"))
  print(paste0("bootset: ", boot.set, ", i: ", i, " done"))  
}

