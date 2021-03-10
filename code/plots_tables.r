library(TMB)
library(Hmisc)

#datasets by species are already compiled by code in .gitignored code/gather_data

parentdir = getwd()
setwd("code")
source("get.best.r")
setwd(parentdir)

x = read.csv(paste0(parentdir,"/data/spp_list.csv"))
x = x[which(x$NSTATIONS>30 & !(x$COMNAME %in% c("FOURSPOT FLOUNDER", "GULF STREAM FLOUNDER"))),]
spps = c(103, 102, 108, 106, 105, 107, 197, 22, 28, 77, 73)
x = x[match(spps,x$SVSPP),]
x$NSTOCKS = c(1,1,2,3,3,1,2,1,1,2,2)
x$sp.names = c("fluke", "plaice", "windowpane", "winter_flounder", "yellowtail_flounder", "witch_flounder", "goosefish", "barndoor", "thorny", "red_hake", "cod")
x$sp.pretty.names = capitalize(tolower(x$COMNAME))
sp.info = x

first.names = paste0(rep(c("bi", "bb"), c(5,8)), c(0:4,0:7))


for(i in 1:NROW(sp.info))
{
  sp = sp.info$sp.names[i]
  combinded.data = readRDS(paste0(parentdir, "data/", sp, "_combined.data.RDS")
  dat = combined.dn.data$catch.length
  
  x = aggregate(cbind(expnumlen.rh,expnumlen.ch) ~ length * day.night, data = combined.dn.data$catch.length, FUN = sum) 
  colnames(x)[3:4] = c("expnumlen.rockhopper", "expnumlen.chainsweep")
  setwd(parentdir)
  write.csv(x, file = paste0(parentdir,"/results/", sp, "/twin_trawl_length_frequencies.csv"))
}
#make summary table of AIC for each model/species
setwd(parentdir)
source("get.best.r")
source("make.all.AIC.table.fn.r")
aic.table = make.all.AIC.table.fn(sp.info$sp.names, sp.info$sp.pretty.names, latex.table.file = paste0(parentdir,"/paper/model_compare.tex"))
models = paste0(rep(c("BI_","BB_"), c(5,8)), c(0:4, 0:7))
rownames(aic.table) = c(models, paste0(models, "_DN"))
colnames(aic.table)[1] = "No. parameters"
write.csv(aic.table, file = paste0(parentdir,"/paper/model_compare.csv"))

setwd(parentdir)
source("make.nobs.table.fn.r")
nobs.table = make.nobs.table.fn()

#Then do bootstrap fits of the models
#set up bootstrap scripts to run on unix servers
fname = "twintrawl_bootstrap_commands.txt"
for(i in 1:NROW(sp.info))
{
  setwd(sp.info$sp.names[i])
  file.copy(paste0(parentdir,"/code/boot_1_script.r"), paste0(parentdir, "/results/",sp))
  file.copy(paste0(parentdir,"/code/boot_2_script.r"), paste0(parentdir, "/results/",sp))
  load(paste0(parentdir,"/results/",sp,"/twintrawl.dn.res.RData")
  x = twintrawl.dn.res
  setwd(parentdir)
  source("get.best.r")
  #convergence issues for the models?
  #aic for the models
  y = get.best(sp.info$sp.names[i])
  setwd(sp.info$sp.names[i])
  best.model = y$best.model 
  best.dn.model = y$best.dn.model
  write("#commands to run on venus in corresponding directory.", file = fname, append = FALSE)
  write("#bootstrap estimates of relative efficiency for from twin trawl data.", file = fname, append = TRUE)
  for(j in 0:9) write(paste0("Rscript --vanilla boot_1_script.r ", j, "23456 100 ",  j, " ", best.model, " &"), file = fname, append = TRUE)
  write("#bootstrap estimates of relative efficiency for day and night from twin trawl data.", file = fname, append = TRUE)
  write("#day", file = fname, append = TRUE)
  for(j in 0:9) write(paste0("Rscript --vanilla boot_2_script.r ", j, "23456 100 1 0 ",  j, " ", best.dn.model, " &"), file = fname, append = TRUE)
  write("#night", file = fname, append = TRUE)
  for(j in 0:9) write(paste0("Rscript --vanilla boot_2_script.r ", j, "23456 100 0 1 ",  j, " ", best.dn.model, " &"), file = fname, append = TRUE)
  setwd(parentdir)
}

#run bootstrap scripts to run on unix servers
#then copy results files back to local directories
#then combine bootstrap results into a single file
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
for(i in 1:NROW(x)) # i = 7 #goosefish
{

  setwd(sp.info$sp.names[i])
  load("boot_pred_eta_0.RData")
  temp = boot.pred.eta
  for(j in 1:9)
  {
    load(paste0("boot_pred_eta_", j, ".RData"))
    temp = rbind(temp, boot.pred.eta) 
  }
  boot.pred.eta = temp
  save(boot.pred.eta, file = "boot.pred.eta.RData")

  load("boot_pred_eta_day_0.RData")
  temp = boot.pred.eta.day
  for(j in 1:9)
  {
    load(paste0("boot_pred_eta_day_", j, ".RData"))
    temp = rbind(temp, boot.pred.eta.day) 
  }
  boot.pred.eta.day = temp
  save(boot.pred.eta.day, file = "boot.pred.eta.day.RData")

  load("boot_pred_eta_night_0.RData")
  temp = boot.pred.eta.night
  for(j in 1:9)
  {
    load(paste0("boot_pred_eta_night_", j, ".RData"))
    temp = rbind(temp, boot.pred.eta.night) 
  }
  boot.pred.eta.night = temp
  save(boot.pred.eta.night, file = "boot.pred.eta.night.RData")
  
  setwd(parentdir)
}

#make plots of results for each species including bootstrap-based confidence intervals.
source("plot.results.fn.r")
source("get.best.r")
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

