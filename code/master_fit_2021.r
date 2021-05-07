library(TMB)
library(Hmisc)

#datasets by species are already compiled by code in .gitignored code/gather_data

parentdir = getwd()
setwd("code")
compile("binom_gamm_nsmooth.cpp")
compile("betabin_gamm_nsmooth.cpp")
dyn.load(dynlib("binom_gamm_nsmooth"))
dyn.load(dynlib("betabin_gamm_nsmooth"))
source("estimate.betabin.gamm.nsmooth.r")
source("estimate.binom.gamm.nsmooth.r")
source("estimate.efficiency.models.r")
source("fit_tmb.R")
source("get.best.r")
setwd(parentdir)

x = read.csv(paste0(parentdir, "/data/spp_list.csv"))
x = x[which(x$NSTATIONS>30 & !(x$COMNAME %in% c("FOURSPOT FLOUNDER", "GULF STREAM FLOUNDER"))),]
spps = c(103, 102, 108, 106, 105, 107, 197, 22, 28, 77, 73)
x = x[match(spps,x$SVSPP),]
x$NSTOCKS = c(1,1,2,3,3,1,2,1,1,2,2)
x$sp.names = c("fluke", "plaice", "windowpane", "winter_flounder", "yellowtail_flounder", "witch_flounder", "goosefish", "barndoor", "thorny", "red_hake", "cod")
x$sp.pretty.names = capitalize(tolower(x$COMNAME))
sp.info = x
saveRDS(sp.info, file = paste0("data/sp.info.RDS"))

system(paste0("mkdir ", parentdir, "/results/big_results"))

#first fit all the models and figure out the best models for each species
for(i in 1:NROW(sp.info))
{
  sp = sp.info$sp.names[i]
  
  spp = sp.info$SVSPP[i]
  print(spp)
  print(sp)
  combined.data = readRDS(paste0(parentdir, "/data/",sp, "_data.RDS"))
  model_fits = estimate.efficiency.models(combined.data$cl.less, do_pred = TRUE, bootstrap = FALSE, predict.length = combined.data$predict.length)
  saveRDS(model_fits, file = paste0(parentdir,"/results/big_results/",sp, "_model_fits.RDS"))
  rm(model_fits)
  rm(combined.data)
  #rm(list = ls()[!(ls() %in% c("stocks","parentdir", "sp.info"))])
  gc()
}

first.names = paste0(rep(c("bi", "bb"), c(5,8)), c(0:4,0:7))


for(i in 1:NROW(sp.info))
{
  print(i)
  sp = sp.info$sp.names[i]
  print(sp)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))
  ind = which(names(x$model.fits) %in% first.names) #remove dn models if already fit.
  x$model.fits = x$model.fits[ind]
  saveRDS(x, file = paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))
  spp = sp.info$SVSPP[i]
  y = get.best(sp, parentdir)
  print(y$best.model) #bi3

  if(!("sdrep" %in% names(x$model.fits[[y$best.model]]$model.res))) {
    x$model.fits[[y$best.model]]$model.res$sdrep = sdreport(x$model.fits[[y$best.model]]$model.res)
    saveRDS(x, file = paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))
  }
  
  combined.data = readRDS(paste0(parentdir, "/data/", sp, "_data.RDS"))
  if(!("dn" %in% names(combined.data$new.dat))){
    combined.data$new.dat = data.frame(dn = factor(rep(c("day","night"), each = length(combined.data$predict.length))), length = rep(combined.data$predict.length,2))
    saveRDS(combined.data, file = paste0(parentdir, "/data/", sp, "_data.RDS"))
  }
  temp = x$model.fits[[y$best.model]]
  best.pars = as.list(x$model.fits[[y$best.model]]$model.res$sdrep, "Estimate")
  do.length.pop = sum(grepl("length", as.character(temp$call[["mu.mean.form"]])))>0
  do.linear.pop = ifelse(do.length.pop, max(abs(best.pars$mean_smooth_re)) < 1e-5, FALSE)
  do.length.sta = sum(grepl("length", as.character(temp$call[["mu.station.form"]])))>0
  do.linear.sta = ifelse(do.length.sta, max(abs(best.pars$station_smooth_re)) < 1e-5, FALSE)

  dat = combined.data$cl.less
  dat$dn = dat$day.night
  is.binom = grepl("bi", y$best.model)
  next.model = ifelse(is.binom, paste0("bi", sum(grepl("bi", names(x$model.fits)))), paste0("bb", sum(grepl("bb", names(x$model.fits))))) 
  #model names start at 0
  # print(next.model)
  # if(do.linear.pop) {
  #   temp$call[["mu.mean.form"]] = cbind(recnumlen.ch, recnumlen.rh) ~ dn + I(length - mean(length))    
  #   temp$call[["use_mean_smooth_re"]] = 0
  # } else {
    if(do.length.pop) {
      temp$call[["mu.mean.form"]] = cbind(recnumlen.ch, recnumlen.rh) ~ dn + s(length, bs = "cr", k = 10)
    } else temp$call[["mu.mean.form"]] = cbind(recnumlen.ch, recnumlen.rh) ~ dn
  # }
  # if(do.linear.sta) {
  #   temp$call[["mu.station.form"]] = cbind(recnumlen.ch, recnumlen.rh) ~ I(length - mean(length))
  #   temp$call[["use_station_smooth_re"]] = 0
  # } #don't have to do anything else at station level. 
  temp$call[["do_pred"]] = TRUE
  temp$call[["new_dat"]] = quote(combined.data$new.dat)
  temp$call[["fit"]] = TRUE
  x$model.fits[[next.model]] = eval(temp$call)
  if(is.character(x$model.fits[[next.model]]$model.res$opt)){ #start at values from next best model if model does not converge. like for plaice
    next.in = x$model.fits[[next.model]]$input
    last.par = as.list(x$model.fits[[y$best.model]]$model.res$sdrep, "Estimate")
    if(length(last.par$beta)<2) {
      last.par$beta = c(last.par$beta,0)
    } else last.par$beta = c(last.par$beta[1], 0, last.par$beta[2:length(last.par$beta)])
    if(length(next.in$par$mean_smooth_re) == 0) last.par$mean_smooth_re = numeric()
    if(length(next.in$par$lambda_par) == 0) last.par$lambda_par = numeric()
    if(length(next.in$par$station_smooth_re) == 0) last.par$station_smooth_re = numeric()
    if(length(next.in$par$lambda_par_station) == 0) last.par$lambda_par_station = numeric()
    next.in$par = last.par
    temp$call[["fit"]] = TRUE
    temp$call[["input"]] = quote(next.in)
    x$model.fits[[next.model]] = eval(temp$call)  
  } else if(x$model.fits[[next.model]]$model.res$opt$conv == 1){ #don't estimate correlation of random effects if it causes lack of convergence
    temp$call[["use_beta_re_cor"]] = 0
    x$model.fits[[next.model]] = eval(temp$call)
  }  
  saveRDS(x, file = paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))
    
  next.model = ifelse(is.binom, paste0("bi", sum(grepl("bi", names(x$model.fits)))), paste0("bb", sum(grepl("bb", names(x$model.fits))))) 
  #model names start at 0
  print(next.model)
  if(do.length.pop) { #if not, then we're all done (like for cod)
    if(do.linear.pop) {
      temp$call[["mu.mean.form"]] = cbind(recnumlen.ch, recnumlen.rh) ~ dn * I(length - mean(length))
      temp$call[["use_mean_smooth_re"]] = 0
    } else temp$call[["mu.mean.form"]] = cbind(recnumlen.ch, recnumlen.rh) ~ dn + s(length, by = factor(dn), bs = "cr", k = 10)
    if("input" %in% names(temp$call)) temp$call = temp$call[-which(names(temp$call) == "input")]
    x$model.fits[[next.model]] = eval(temp$call)
    if(is.character(x$model.fits[[next.model]]$model.res$opt)){ #start at values from next best model if model does not converge
      next.in = x$model.fits[[next.model]]$input
      last.par = as.list(x$model.fits[[y$best.model]]$model.res$sdrep, "Estimate")
      if(length(last.par$beta)<3) last.par$beta = c(last.par$beta[1],0,last.par$beta[2],0)
      else last.par$beta = c(last.par$beta[1], 0, last.par$beta[2],0, last.par$beta[3:length(last.par$beta)])
      last.par$mean_smooth_re = rep(last.par$mean_smooth_re,2)
      last.par$lambda_par = rep(last.par$lambda_par,2)
      next.in$par = last.par
      temp$call[["fit"]] = TRUE
      temp$call[["input"]] = quote(next.in)
      x$model.fits[[next.model]] = eval(temp$call)
    } else if(x$model.fits[[next.model]]$model.res$opt$conv == 1){ #don't estimate correlation of random effects if it causes lack of convergence
      temp$call[["use_beta_re_cor"]] = 0
      x$model.fits[[next.model]] = eval(temp$call)
    }  
  }
  saveRDS(x, file = paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))
  y = get.best(sp,parentdir)
  print(y$best.model)

  if(!("sdrep" %in% names(x$model.fits[[y$best.model]]$model.res))) {
    x$model.fits[[y$best.model]]$model.res$sdrep = sdreport(x$model.fits[[y$best.model]]$model.res)
  }
  x$model.fits = x$model.fits[c(sort(grep("bi", names(x$model.fits), value = TRUE)), sort(grep("bb", names(x$model.fits), value = TRUE)))]
  saveRDS(x, file = paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))
  remove(x)
  remove(combined.data)
}

fname = paste0(parentdir, "/code/twintrawl_bootstrap_commands.txt")
write("#commands to run on venus in corresponding directory.", file = fname, append = FALSE)
write("#bootstrap estimates of relative efficiency for from twin trawl data.", file = fname, append = TRUE)
for(i in 1:NROW(sp.info))
{
  print(i)
  sp = sp.info$sp.names[i]
  print(sp)

  #setwd(sp.info$sp.names[i])
  #setwd(parentdir)
  #file.copy("boot_script.r", sp)
  x = readRDS(paste0(parentdir, "/results/big_results/", sp,"_model_fits.RDS"))
  #convergence issues for the models?
  #aic for the models
  best.model = get.best(sp, parentdir)$best.model 
  for(j in 0:9) write(paste0("Rscript --vanilla code/boot_script.r ", j, "23456 100 ",  j, " ", best.model, " ", i, " &"), file = fname, append = TRUE)
}


#push all info needed for bootstrapping to server
#rsync -arvxu ~/work/paired_tow_studies/R/2020/winter_flounder saturn.nefsc.noaa.gov:/home7/tmiller2/work/paired_tow_studies/R/2020

#pull all bootstrap results back to local machine
#rsync -arvxu saturn.nefsc.noaa.gov:/home7/tmiller2/work/paired_tow_studies/R/2020/winter_flounder ~/work/paired_tow_studies/R/2020

# setwd(sp)
# load("boot_pred_eta_0.RData")
# temp = boot.pred.eta
# for(j in 1:9)
# {
#   load(paste0("boot_pred_eta_", j, ".RData"))
#   temp = rbind(temp, boot.pred.eta) 
# }
# boot.pred.eta = temp
# saveRDS(boot.pred.eta, file = "boot.pred.eta.RDS")
# setwd(parentdir)
  
# source("plot.results.r")
# plot.results("red_hake", ymax = c(40,20), xlim = c(10,50)) #red hake



#make plots of results for each species including bootstrap-based confidence intervals.
plot.results.fn(1) 
for(i in 2:6) plot.results.fn(i) 
plot.results.fn(7) #goosefish
plot.results.fn(3, 10, sp.info) #windowpane, need wider y-axis
plot.results.fn(8, ymax = 8) #barndoor skate
plot.results.fn(9) #thorny skate
plot.results.fn(5, ymax = 10) #yellowtail flounder
plot.results.fn(6, ymax = 10) #witch flounder
plot.results.fn(10, ymax = 10) #red hake
plot.results.fn(11,6) #cod


#make summary table of AIC for each model/species
setwd(parentdir)
source("get.best.r")
source("make.all.AIC.table.fn.r")
aic.table = make.all.AIC.table.fn(sp.info$sp.names, sp.info$sp.pretty.names, latex.table.file = paste0(parentdir, "/paper/model_compare.tex"))
models = paste0(rep(c("BI_","BB_"), c(5,8)), c(0:4, 0:7))
rownames(aic.table) = c(models, paste0(models, "_DN"))
write.csv(aic.table, file = 'model_compare.csv')

#plots for working paper to fluke assessment
#setwd(parentdir)
#source("get.best.r")
#source("plot.all.rho.results.fn.r")
#plot.all.rho.results.fn(c(6,6,10,10,10,10))
#setwd(parentdir)
#source("get.best.r")
#source("plot.all.rho.dn.results.fn.r")
#plot.all.rho.results.fn(c(6,6,15,6,10,10))

