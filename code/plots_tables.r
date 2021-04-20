library(TMB)
library(Hmisc)
load("data/survey/stock.strata.RData")
source("code/get.best.r")
source("code/make.all.AIC.table.fn.r")
parentdir = getwd()

#datasets by species are already compiled by code in .gitignored code/gather_data


for(i in 1:NROW(sp.info))
{
  sp = sp.info$sp.names[i]
  combined.data = readRDS(paste0(parentdir, "/data/", sp, "_data.RDS"))
  dat = combined.data$catch.length  
  x = aggregate(cbind(expnumlen.rh,expnumlen.ch) ~ length * day.night, data = combined.data$catch.length, FUN = sum) 
  colnames(x)[3:4] = c("expnumlen.rockhopper", "expnumlen.chainsweep")
  write.csv(x, file = paste0(parentdir,"/results/", sp, "_twin_trawl_length_frequencies.csv"))
}


source("make.nobs.table.fn.r")
nobs.table = make.nobs.table.fn()


#make summary table of AIC for each model/species
species.get.best = lapply(sp.info$sp.names, get.best, pdir = parentdir)
all.models = unique(unlist(sapply(x, function(y) names(y$np))))

n.rep = c(5,8)
first.names = paste0(rep(c("bi", "bb"), n.rep), c(0:4,0:7))

out = cbind.data.frame(
  Model = paste0(rep(c("BI", "BB"), n.rep), "$_{", c(0:(n.rep[1]-1),0:(n.rep[2]-1)), "}$"),
  "$\\log\\left(\\rho\\right)$" = paste0("$\\sim$ ", c(
  1, 
  "1 + 1$|$pair", 
  "s(length)", 
  "s(length) + 1$|$pair",
  "s(length) + s(length)$|$pair",
  1,
  "1 + 1$|$pair", 
  "s(length)", 
  "s(length)", 
  "s(length) + 1$|$pair",
  "s(length) + 1$|$pair",
  "s(length) + s(length)$|$pair",
  "s(length) + s(length)$|$pair")),
  "$\\log\\left(\\phi\\right)$" = c(
  "--", 
  "--", 
  "--", 
  "--", 
  "--", 
  paste0("$\\sim$ ", c(
  1, 
  1, 
  1, 
  "s(length)", 
  1, 
  "s(length)", 
  1, 
  "s(length)"))),  
#  "log-likelihood" = round(nll.rep,2),
  "$n_p$" = species.get.best[[1]]$np[first.names]#,
  #aic = round(aic.rep,2)
#  "log-likelihood" = round(-sapply(x$model.fits, function(x) x$model.res$opt$obj),2),
#  "$P$" = sapply(x$model.fits, function(x) length(x$model.res$opt$par)),
#  aic = round(sapply(x$model.fits, function(x) 2*(x$model.res$opt$obj + length(x$model.res$opt$par))),2)
)
#out[["$\\Delta$AIC"]] = round(out$aic - min(out$aic, na.rm = TRUE), 2)
out$Description = c(
  "population-level mean for all observations",
  "population- and random station-level $\\rho$",
  "population-level smooth size effect on $\\rho$",
  "population-level smooth size effect and random station-level intercept for $\\rho$",
  "population-level and random station-level smooth size effects for $\\rho$",
  "population-level $\\rho$ and $\\phi$",
  "population-level and random station-level intercept for $\\rho$ and population-level $\\phi$",
  "population-level smooth size effect on $\\rho$ and population-level $\\phi$",
  "population-level smooth size effect on $\\rho$ and $\\phi$",
  "population-level smooth size effect and random station-level intercept for $\\rho$ and population-level $\\phi$",
  "population-level smooth size effect on $\\rho$ and $\\phi$ and random station-level intercepts for $\\rho$",
  "population-level and random station-level smooth size effects on $\\rho$ and population-level $\\phi$",
  "population-level and random station-level smooth size effects on $\\rho$ and population-level smooth size effects on $\\phi$")#,
#  "population-level smooth size and day/night effects and random station-level smooth size effects on $\\rho$ and population-level smooth size effects on $\\phi$",
#  "BB$_8$ with no correlation of station-specific random effects",
#  "BB$_9$ with interaction of population-level day/night and smooth size effects")

x = latex(out, file = paste0(parentdir, "/paper/base_model_description.tex"), 
  table.env = FALSE, col.just = c("l","l","l","r","p{0.5\\textwidth}"), rowname = NULL)

aic.table = make.all.AIC.table.fn(sp.info$sp.names, sp.info$sp.pretty.names, pdir = parentdir)
aic.table = t(aic.table[,-12])
aic.table = aic.table[c(2:7,11,8:10),]
aic.table[11,] = round(aic.table[11,], 0)
cell.format = matrix("", NROW(aic.table), NCOL(aic.table))
cell.format[cbind(1:10, apply(aic.table[1:10,], 1, function(x) which(x == 0)))] = "bfseries"
x = aic.table
x[] = as.character(x)
x = latex(x, file = paste0(parentdir,"/paper/base_model_compare.tex"), 
  #rgroup = c("", ""), first.hline.double = FALSE,
  #n.rgroup = c(length(sp.info$sp.names)-1,1), 
  col.just = rep("r",NCOL(x)),
  numeric.dollar = FALSE, cellTexCmds = cell.format,
  rowlabel = '', table.env = FALSE)#, rowlabel.just = "c")


source("code/make_best_model_table.r") #makes a data.frame called finalout
ind = c(1:6,10,7:9)
temp = finalout[,-1]
cell.format = matrix("", NROW(temp), NCOL(temp))
cell.format[temp[[5]] == "0",5] = "bfseries"
x = latex(temp, file = paste0(parentdir,"/paper/best_model_compare.tex"), 
  rgroup = sp.info$sp.pretty.names[ind],
  n.rgroup = rep(3,length(ind)), 
  col.just = c("l", "l", "l", rep("r",2)),
  numeric.dollar = FALSE, cellTexCmds = cell.format,
  rowlabel = '', rowname = rep("",NROW(temp)), table.env = FALSE)#, rowlabel.just = "c")


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
