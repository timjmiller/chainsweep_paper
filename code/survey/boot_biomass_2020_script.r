#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(methods)
library(TMB)
library(RODBC)


#args[1] is the number of bootstraps
#args[2] is the species index
#args[3] is the stock index
#args[4] is the year
#args[5] is whether to do spring (1)
#args[6] is whether to do fall (1)

source("code/survey/boot.trac.biomass.dn.fn.r")
source("code/survey/boot.lendat.fn.r")
source("code/survey/get.survey.Nal.hat.fn.r")
source("code/survey/boot.dfo.lendat.fn.r")
source("code/survey/get.dfo.survey.Nal.hat.fn.r")  
source("code/survey/boot_biomass_2020.r")
dyn.load(dynlib("code/survey/wal"))
load("data/stock.strata.RData")

parentdir = getwd()

#setwd(sp.info$sp.names[as.integer(args[2])])
boot_biomass(sp.i = as.integer(args[2]), stock.i = as.integer(args[3]), do.boot.biomass = TRUE, 
  do.boot.res = FALSE, years = as.integer(args[4]), n.boot = as.integer(args[1]), do.spring = as.integer(args[5]) == 1,
  do.fall = as.integer(args[6]) == 1, summ.file.mod = paste0(parentdir,"/results/"))
setwd(parentdir)
#dyn.unload(dynlib("wal"))
