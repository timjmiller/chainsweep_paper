#pull all ADIOS survey strata info
#home0/pdy/pub/STOCKEFF/ADIOS/ADIOS_SV/website/webfiles
#rsync -arvxu saturn.nefsc.noaa.gov:/home0/pdy/pub/STOCKEFF/ADIOS/ADIOS_SV/website/webfiles ~/work/paired_tow_studies/R

source("code/get.best.r")
source("code/survey/estimate_biomass.r")
source("code/survey/trac.biomass.dn.fn.r")
source("code/survey/get.survey.Nal.hat.fn.r")
load("data/survey/stock.strata.RData")

parentdir = getwd() #main repo directory, exists in stock.strata.RData for the server so need to make sure to use the right one.
setwd("code/survey")
library(TMB)
compile("wal.cpp")
dyn.load(dynlib("wal"))
setwd(parentdir)

file.name.prefix = paste0(parentdir, "/results/")
survey.years = 2009:2019
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)

#make estimates of biomass and numbers at length by converting from rockhopper to chainsweep gear.
for(i in 1:length(x))
{
  print(stocks[i])
  #30+ biomass for GOM winter flounder
  twintrawl.res = readRDS(paste0("results/big_results/", sp.info$sp.names[x[i]], "_model_fits.RDS"))
  if(stocks[i] == "gom_winter_flounder") estimate_biomass(sp.i = x[i], stock.i = i, sp.info, twintrawl.res, 
    years = survey.years, filename.mod = file.name.prefix, min.length.calibrate = 30, pdir = parentdir)
	 #goosefish uses different TOGA definition
  if(stocks[i] %in% c("north_goose", "south_goose")) estimate_biomass(sp.i = x[i], stock.i = i, sp.info, twintrawl.res, 
    years = survey.years, filename.mod = file.name.prefix,  TOGA.operation = 2, pdir = parentdir) 
  if(!(stocks[i] %in% c("gom_winter_flounder", "north_goose", "south_goose"))) estimate_biomass(sp.i = x[i], stock.i = i, sp.info, twintrawl.res, 
    years = survey.years, filename.mod = file.name.prefix, pdir = parentdir)
}


##############################################
#do bootstraps for each stock
#on venus server  using commands in unix_boot_biomass_commands.txt
#move results back to local machine
##############################################

#make the unix commands to run on venus. copy and paste from this file.
fname = paste0(parentdir"bootstrap_biomass_commands.txt"
write("#commands for bootstrapping biomass estimates on venus in corresponding directory.", file = fname, append = FALSE)
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
for(i in 1:length(x))
{
  for(d in 1:2) {
    for(y in survey.years) {
      write(paste0("Rscript --vanilla boot_biomass_2019_script.r 1000 ", x[i], " ", i, " ", y, " ", d, " 1 1 &"), file = fname, append = TRUE)
    }
    #fall not available yet for 2019
    #write(paste0("Rscript --vanilla boot_biomass_2019_script.r 1000 ", x[i], " ", i, " ", 2019, " ", d, " 1 0 &"), file = fname, append = TRUE)
  }
  write("\n", file = fname, append = TRUE)
}

#push all info needed for bootstrapping to serveri
#rsync -arvxu ~/work/paired_tow_studies saturn.nefsc.noaa.gov:/home7/tmiller2/work

#pull all bootstrap results back to local machine
#rsync -arvxu saturn.nefsc.noaa.gov:/home7/tmiller2/work/paired_tow_studies/R/2019 ~/work/paired_tow_studies/R

#make summary files
source("boot_biomass_2019.r")
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
for(i in 1:length(x)) {
  boot_biomass(sp.i = x[i], stock.i = i, do.boot.biomass.1 = FALSE, do.boot.biomass.2 = FALSE, 
    do.boot.1.res = TRUE, do.boot.2.res = TRUE, years = 2009:2019, n.boot = 1000, do.spring = TRUE,
    do.fall = TRUE, summ.file.mod = "2019_")
}
#parentdir = getwd()

setwd(parentdir)
load("stock.strata.RData")
source("plot.biomass.fn.r")
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
for( i in 1:length(x)) plot.biomass.fn(i=x[i], sp.info, stock = stocks[i], use.stock.names[i], file.loc = "~/work/paired_tow_studies/R/2019")
plot.biomass.fn(i=1, sp.info, stock = stocks[1], use.stock.names[1], file.loc = "~/work/paired_tow_studies/R/2019")

load("stock.strata.RData")
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
y = read.csv("model_compare.csv", as.is = TRUE)
for( i in 1:NROW(sp.info)) {
  system(paste0("mkdir ../assessments/",sp.info$sp.names[i]))
  dn = which(y[,2+i] == 0)>13
  setwd(sp.info$sp.names[i])
  if(dn) ftypes <- c("\\.2\\.", "2\\.csv", "2\\.png", "dn\\_rho") else ftypes <- c("\\.1\\.", "1\\.csv", "1\\.png", "^(?=.*rho)(?!.*dn)")
  fnames = unique(unlist(sapply(ftypes, grep, x = dir(), perl = TRUE, value = TRUE)))
  system(paste0("cp ", paste(fnames, collapse = " "), " ", "../../assessments/", sp.info$sp.names[i]))
  setwd(parentdir)
}

load("stock.strata.RData")
setwd(parentdir)
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
for( i in 1:length(stocks)) {
  setwd(sp.info$sp.names[x[i]])
  load(paste0(stocks[i], "_N.W.RData"))
  s.ind = match(as.character(2009:2019),colnames(all.spring.N.W))
  f.ind = match(as.character(2009:2019),colnames(all.fall.N.W))
  y = cbind(spring = all.spring.N.W[2,s.ind], fall = all.fall.N.W[2,f.ind])/1000
  write.csv(y, file = paste0("../../assessments/", sp.info$sp.names[x[i]], "/", stocks[i], "_bigelow_biomass.csv"))
  y = cbind(spring = all.spring.N.W[3,s.ind], fall = all.fall.N.W[3,f.ind])
  write.csv(y, file = paste0("../../assessments/", sp.info$sp.names[x[i]], "/", stocks[i], "_bigelow_n_per_tow.csv"))
  y = cbind(spring = all.spring.N.W[4,s.ind], fall = all.fall.N.W[4,f.ind])
  write.csv(y, file = paste0("../../assessments/", sp.info$sp.names[x[i]], "/", stocks[i], "_bigelow_kg_per_tow.csv"))
  setwd(parentdir)
}

#push all results to assessment folderss on server
#rsync -arvxu ~/work/paired_tow_studies/R/assessments saturn.nefsc.noaa.gov:/home7/tmiller2/work/paired_tow_studies/R
