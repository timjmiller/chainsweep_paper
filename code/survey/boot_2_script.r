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
do.day = as.logical(as.integer(args[3]))
do.night = as.logical(as.integer(args[4]))
boot.set = as.integer(args[5])
model = args[6]

cl.less.day = combined.dn.data$cl.less.day
cl.less.night = combined.dn.data$cl.less.night
cl.aug.day = combined.dn.data$cl.less.day
cl.aug.night = combined.dn.data$cl.less.night
pred.length = combined.dn.data$predict.length

if(do.night) boot.pred.eta.night = matrix(NA, n, length(pred.length))
if(do.day) boot.pred.eta.day = matrix(NA, n, length(pred.length))

for(i in 1:n)
{
  print(paste0("i: ", i))
  if(do.night)
  {
    night = final.estimate.efficiency.fn(cl.less.night, cl.aug.night, fit.models = model, bootstrap = TRUE,
      predict.length = pred.length)
    night$sdrep = try(summary(sdreport(night$model.fits[[model]]$model.res)))
    if(!is.character(night$sdrep)) boot.pred.eta.night[i,] = night$sdrep[which(rownames(night$sdrep) == "mean_pred_eta"),1]
    save(boot.pred.eta.night, file = paste0("boot_pred_eta_night_", boot.set, ".RData"))
  }
  if(do.day)
  {
    day = final.estimate.efficiency.fn(cl.less.day, cl.aug.day, fit.models = model, bootstrap = TRUE,
      predict.length = pred.length)
    day$sdrep = try(summary(sdreport(day$model.fits[[model]]$model.res)))
    if(!is.character(day$sdrep)) boot.pred.eta.day[i,] = day$sdrep[which(rownames(day$sdrep) == "mean_pred_eta"),1]
    save(boot.pred.eta.day, file = paste0("boot_pred_eta_day_", boot.set, ".RData"))
  }
  print(paste0("day: ", do.day, ", night: ", do.night, ", bootset: ", boot.set, ", i: ", i, " done"))  
}

