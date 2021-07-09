#pull all ADIOS survey strata info
#home0/pdy/pub/STOCKEFF/ADIOS/ADIOS_SV/website/webfiles
#rsync -arvxu saturn.nefsc.noaa.gov:/home0/pdy/pub/STOCKEFF/ADIOS/ADIOS_SV/website/webfiles ~/work/paired_tow_studies/R

source("code/get.best.r")
source("code/survey/estimate_biomass.r")
source("code/survey/trac.biomass.dn.fn.r")
source("code/survey/get.survey.Nal.hat.fn.r")
source("code/survey/boot_biomass.r")
source("code/survey/plot.biomass.fn.r")

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
fname = "code/survey/bootstrap_biomass_commands_shell.txt"
write("#commands for bootstrapping biomass estimates on venus in corresponding directory.", file = fname, append = FALSE)
for(i in 1:length(x))
{
  for(y in survey.years) {
    write(paste0("Rscript --vanilla code/survey/boot_biomass_script.r 1000 ", x[i], " ", i, " ", y, " 1 1 &"), file = fname, append = TRUE)
  }
}

#push all info needed for bootstrapping to server

#pull all bootstrap results back to local machine

