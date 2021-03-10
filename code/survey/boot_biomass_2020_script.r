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

load("stock.strata.RData")
parentdir = getwd()
dyn.load(dynlib("wal"))
source("boot_biomass_2020.r")
setwd(sp.info$sp.names[as.integer(args[2])])
boot_biomass(sp.i = as.integer(args[2]), stock.i = as.integer(args[3]), do.boot.biomass = TRUE, 
  do.boot.res = FALSE, years = as.integer(args[4]), n.boot = as.integer(args[1]), do.spring = as.integer(args[5]) == 1,
  do.fall = as.integer(args[6]) == 1, summ.file.mod = "2020_")
setwd(parentdir)
#dyn.unload(dynlib("wal"))
