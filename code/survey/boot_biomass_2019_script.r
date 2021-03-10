#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(methods)
library(TMB)
library(RODBC)

#args[1] is the number of bootstraps
#args[2] is the species index
#args[3] is the stock index
#args[4] is the year
#args[5] is whether to use day/night (2) model or not (1)
#args[6] is whether to do spring (1)
#args[7] is whether to do fall (1)

parentdir = getwd()
#load("~/work/paired_tow_studies/R/2017/.RData")
load("stock.strata.RData")
source("boot_biomass_2019.r")
print(as.integer(args[1:7]))
print(as.integer(args[5]) == 2)
if(as.integer(args[5]) == 2) {
  boot_biomass(sp.i = as.integer(args[2]), stock.i = as.integer(args[3]), do.boot.biomass.1 = FALSE, do.boot.biomass.2 = TRUE, 
    do.boot.1.res = FALSE, do.boot.2.res = FALSE, years = as.integer(args[4]), n.boot = as.integer(args[1]), do.spring = as.integer(args[6]) == 1,
    do.fall = as.integer(args[7]) == 1, summ.file.mod = "2019_")
} else {
  boot_biomass(sp.i = as.integer(args[2]), stock.i = as.integer(args[3]), do.boot.biomass.1 = TRUE, do.boot.biomass.2 = FALSE, 
    do.boot.1.res = FALSE, do.boot.2.res = FALSE, years = as.integer(args[4]), n.boot = as.integer(args[1]), do.spring = as.integer(args[6]) == 1,
    do.fall = as.integer(args[7]) == 1, summ.file.mod = "2019_")
}
