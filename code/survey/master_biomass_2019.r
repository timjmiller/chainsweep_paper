#pull all ADIOS survey strata info
#home0/pdy/pub/STOCKEFF/ADIOS/ADIOS_SV/website/webfiles
#rsync -arvxu saturn.nefsc.noaa.gov:/home0/pdy/pub/STOCKEFF/ADIOS/ADIOS_SV/website/webfiles ~/work/paired_tow_studies/R

#make estimates of biomass and numbers at length by converting from rockhopper to chainsweep gear.
library(RODBC)
creds = readLines("/home/tmiller2/work/oracle_uid_pwd.txt")
sole <- odbcConnect(dsn="sole", uid=creds[1], pwd=creds[2], believeNRows=FALSE)
source("define.stock.strata.r")
remove(sole)
save(spring.strata.sizes,fall.strata.sizes, fall.strata, spring.strata, big.swept.area, stocks, use.stock.names, sp.info, parentdir, file = "stock.strata.RData")

  estimate_biomass(sp.i = x[i], stock.i = i, sp.info, sole, do.lendat = TRUE, do.lw.dat = TRUE, do.N.W = TRUE, years = 2009, filename.mod = "test_") 
  source("estimate_biomass_2019.r")
sole <- odbcConnect(dsn="sole", uid=creds[1], pwd=creds[2], believeNRows=FALSE)
estimate_biomass(sp.i = x[i], stock.i = i, sp.info, sole, do.lendat = TRUE, do.lw.dat = TRUE, do.N.W = TRUE, 
	years = 2009:2019, filename.mod = "2019_", min.length.calibrate = 30)
  setwd(parentdir)
  odbcCloseAll()
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
for(i in 1:length(x))
{
  setwd(parentdir)
  load("stock.strata.RData")
  print(stocks[i])
  source("get.best.r")
  source("estimate_biomass_2019.r")
  sole <- odbcConnect(dsn="sole", uid=creds[1], pwd=creds[2], believeNRows=FALSE) 
  #30+ biomass for GOM winter flounder
  if(stocks[i] == "gom_winter_flounder") estimate_biomass(sp.i = x[i], stock.i = i, sp.info, sole, do.lendat = TRUE, do.lw.dat = TRUE, do.N.W = TRUE, 
	years = 2009:2019, filename.mod = "2019_", min.length.calibrate = 30)
	 #goosefish uses different TOGA definition
  if(stocks[i] %in% c("north_goose", "south_goose")) estimate_biomass(sp.i = x[i], stock.i = i, sp.info, sole, do.lendat = TRUE, do.lw.dat = TRUE, do.N.W = TRUE, 
	years = 2009:2019, filename.mod = "2019_",  TOGA.operation = 2) 
  if(!(stocks[i] %in% c("gom_winter_flounder", "north_goose", "south_goose"))) estimate_biomass(sp.i = x[i], stock.i = i, sp.info, sole, do.lendat = TRUE, do.lw.dat = TRUE, do.N.W = TRUE, 
	years = 2009:2019, filename.mod = "2019_")
  setwd(parentdir)
  odbcCloseAll()
}

#after fall comes in
y = 2019
for(i in 1:length(x))
{
  setwd(parentdir)
  load("stock.strata.RData")
  print(stocks[i])
  source("get.best.r")
  source("estimate_biomass_2019.r")
  sole <- odbcConnect(dsn="sole", uid=creds[1], pwd=creds[2], believeNRows=FALSE) #for some reason does not work in Rstudio on linux.
  if(stocks[i] == "gom_winter_flounder") estimate_biomass(sp.i = x[i], stock.i = i, sp.info, sole, do.lendat = TRUE, do.lw.dat = TRUE, do.N.W = TRUE, 
	years = y, filename.mod = paste0(y,"_"), min.length.calibrate = 30) else 
  estimate_biomass(sp.i = x[i], stock.i = i, sp.info, sole, do.lendat = TRUE, do.lw.dat = TRUE, do.N.W = TRUE, years = y, filename.mod = paste0(y,"_"))
  setwd(parentdir)
  odbcCloseAll()
}

##############################################
#do bootstraps for each stock
#on venus server  using commands in unix_boot_biomass_commands.txt
#move results back to local machine
##############################################

#make the unix commands to run on venus. copy and paste from this file.
fname = "bootstrap_biomass_commands_2019.txt"
write("#commands for bootstrapping biomass estimates on venus in corresponding directory.", file = fname, append = FALSE)
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
for(i in 1:length(x))
{
  for(d in 1:2) {
    for(y in 2009:2018) {
      write(paste0("Rscript --vanilla boot_biomass_2019_script.r 1000 ", x[i], " ", i, " ", y, " ", d, " 1 1 &"), file = fname, append = TRUE)
    }
    #fall not available yet for 2019
    write(paste0("Rscript --vanilla boot_biomass_2019_script.r 1000 ", x[i], " ", i, " ", 2019, " ", d, " 1 0 &"), file = fname, append = TRUE)
  }
  write("\n", file = fname, append = TRUE)
}

#after fall comes in
years = 2019
seasons = c(1,1) #1 = spring, 2 = fall
fname = paste0("bootstrap_biomass_commands_", years, ".txt")
write("#commands for bootstrapping biomass estimates on venus in corresponding directory.", file = fname, append = FALSE)
x = rep(1:NROW(sp.info), sp.info$NSTOCKS)
for(y in years) {
  for(d in 1:2) {
    for(i in 1:length(x))
    {
      write(paste("Rscript --vanilla boot_biomass_2019_script.r 1000", x[i], i, y, d, paste(seasons, collapse = " "), collapse = " "), file = fname, append = TRUE)
    }
    write("\n", file = fname, append = TRUE)
  }
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
