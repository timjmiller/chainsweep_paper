library(TMB)
library(Hmisc)
load("data/survey/stock.strata.RData")
source("code/get.best.r")
source("code/make.all.AIC.table.fn.r")
parentdir = getwd()

#datasets by species are already compiled by code in .gitignored code/gather_data

species.order = c(1:6,10,7:9)
sp.info$sp.pretty.names[species.order]

for(i in 1:NROW(sp.info))
{
  sp = sp.info$sp.names[i]
  combined.data = readRDS(paste0(parentdir, "/data/", sp, "_data.RDS"))
  dat = combined.data$catch.length  
  x = aggregate(cbind(expnumlen.rh,expnumlen.ch) ~ length * day.night, data = combined.data$catch.length, FUN = sum) 
  colnames(x)[3:4] = c("expnumlen.rockhopper", "expnumlen.chainsweep")
  write.csv(x, file = paste0(parentdir,"/results/", sp, "_twin_trawl_length_frequencies.csv"))
}


source("code/make.nobs.table.fn.r")
nobs.table = make.nobs.table.fn()
nobs.table = nobs.table[species.order,]
nobs.table = format.df(nobs.table, big.mark = ",", numeric.dollar=F)
x = latex(nobs.table, file = paste0(parentdir,"/paper/data_table.tex"), 
  cgroup = c("Paired Tows", "Captured", "\\thead{Both Gears \\\\ Measured}", "\\thead{Chainsweep \\\\ Measured}", "\\thead{Rockhopper \\\\ Measured}"), first.hline.double = FALSE,
  n.cgroup = c(3,1,3,3,3), 
  col.just = rep("r",NCOL(nobs.table)),
  cgroupTexCmd="normalfont",#numeric.dollar = FALSE, #cellTexCmds = cell.format,
  rowlabel = 'Species', table.env = FALSE)#, rowlabel.just = "c")

#make summary table of AIC for each model/species
species.get.best = lapply(sp.info$sp.names, get.best, pdir = parentdir)
all.models = unique(unlist(sapply(species.get.best, function(y) names(y$np))))

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

# x = readRDS(paste0(parentdir,"/results/big_results/",sp.info$sp.names[i],"_model_fits.RDS"))
# best.model = species.get.best[[i]]$best.model
# ind = length(x$pred.length)
# #png(filename = paste0("~/work/paired_tow_studies/R/2020/",sp,"/plot_for_Dave.png"), width = 12*144, height = 12*144, res = 144, pointsize = 12, family = "Times")#,
# plot(x$pred.length,x$model.fits$bi2$model.res$rep$mean_pred_eta[1:ind], type = 'n', ylim = c(-1,log(20)), xlab = "Length", ylab = "Ln Relative Efficiency")
# grid()
# lines(x$pred.length,x$model.fits[[best.model]]$model.res$rep$mean_pred_eta[1:ind], col = 'red')
# lines(x$pred.length,x$model.fits[[best.model]]$model.res$rep$mean_pred_eta[1:ind + ind], col = 'blue')

# x = readRDS(paste0("data/", sp.info$sp.names[i], "_data.RDS"))$catch.length
# x$ratio = x$expnumlen.ch/x$expnumlen.rh
# x$ratio[is.infinite(x$ratio)] = NA
# x$ratio[x$ratio == 0] = NA
# y = aggregate(ratio ~ length * day.night, data = x, FUN = mean, na.rm = TRUE)
# points(y$length[y$day.night == 'day'], log(y[y$day.night == 'day',3]), col = 'red', pch = 19)
# points(y$length[y$day.night == 'night'], log(y[y$day.night == 'night',3]), col = 'blue', pch = 19)
# y = aggregate(cbind(expnumlen.ch, expnumlen.rh) ~ length * day.night, data = x, FUN = sum)
# points(y$length[y$day.night == 'day'], log(y[,3]/y[,4])[y$day.night == 'day'], col = 'red', pch = 2)
# points(y$length[y$day.night == 'night'], log(y[,3]/y[,4])[y$day.night == 'night'], col = 'blue', pch = 2)
# #dev.off()

#done to here

i = 2
sp.name = sp.info$sp.names[i]
x = readRDS(paste0(parentdir,"/results/big_results/",sp.name,"_model_fits.RDS"))
best.model = species.get.best[[i]]$best.model
x = summary(x$model.fits[[best.model]]$model.res$sdrep)
x = x[which(rownames(x) == "mean_pred_eta"),]
temp = readRDS(paste0("/home/tmiller2/work/paired_tow_studies/R/2020/",sp.name,"/model_fits.RDS"))
temp = summary(temp$model.fits[[best.model]]$model.res$sdrep)
temp = temp[which(rownames(temp) == "mean_pred_eta"),]

#make plots of results for each species including bootstrap-based confidence intervals.
source("code/plot.results.r")
#source("get.best.r")
n.good.boot = c()
for(i in species.order) 
{
  boot.pred.eta = readRDS(paste0("results/", sp.info$sp.names[i], "_boot_pred_eta_0.RDS"))
  for(j in 1:9) boot.pred.eta = rbind(boot.pred.eta, readRDS(paste0("results/", sp.info$sp.names[i], "_boot_pred_eta_",j,".RDS")))
  n.good.boot[i] = sum(!is.na(boot.pred.eta[,1]))
}

boot.pred.eta = readRDS(paste0("results/", sp.name, "_boot_pred_eta_0.RDS"))
#temp = boot.pred.eta
for(j in 1:9) boot.pred.eta = rbind(boot.pred.eta, readRDS(paste0("results/", sp.name, "_boot_pred_eta_",j,".RDS")))
temp = readRDS(paste0("/home/tmiller2/work/paired_tow_studies/R/2020/",sp.name,"/boot.pred.eta.RDS"))
apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.025, na.rm = TRUE)
apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.5, na.rm = TRUE)
apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.975, na.rm = TRUE)
    if(boot) lines(plen, exp(apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.975, na.rm = TRUE)), col = 'black', lty = 2)

ymax.sp = c(6,6,10,6,10,10,6,10,6,10)
png(filename = paste0("paper/sp_rho_plot.png"), width = 8*144, height = 12*144, res = 144, pointsize = 12, family = "Times")
par(oma = c(3,3,0,0), mar = c(2,2,3,1), mfcol = c(5,2))
for(i in species.order) 
{
  plot.results(sp.info$sp.names[i], pdir = parentdir, ymax = ymax.sp[i]) 
  mtext(side = 3, sp.info$sp.pretty.names[i], line = 0, cex = 1.5)
}
mtext(side = 1, 'Length (cm)', line = 1, cex = 1.5, outer = TRUE)
mtext(side = 2, "Relative Catch Efficiency (Chain:Rockhopper)", line = 1, cex = 1.5, outer = TRUE)
dev.off()

plot.results(sp.info$sp.names[1], pdir = parentdir) 
mtext(side = 3, sp.info$sp.pretty.names[1], line = 0, cex = 1.5)
plot.results(sp.info$sp.names[2], pdir = parentdir) 
mtext(side = 3, sp.info$sp.pretty.names[2], line = 0, cex = 1.5)
plot.results(sp.info$sp.names[3], pdir = parentdir, ymax = 10) 
mtext(side = 3, sp.info$sp.pretty.names[3], line = 0, cex = 1.5)
plot.results(sp.info$sp.names[4], pdir = parentdir) 
mtext(side = 3, sp.info$sp.pretty.names[4], line = 0, cex = 1.5)
plot.results(sp.info$sp.names[5], pdir = parentdir) 
mtext(side = 3, sp.info$sp.pretty.names[5], line = 0, cex = 1.5)
for(i in species.order) plot.results.fn(i) 
plot.results.fn(7) #goosefish
plot.results.fn(3, 10, sp.info) #windowpane, need wider y-axis
plot.results.fn(8, ymax = 8) #barndoor skate
plot.results.fn(9) #thorny skate
plot.results.fn(5, ymax = 10) #yellowtail flounder
plot.results.fn(6, ymax = 10) #witch flounder
plot.results.fn(10, ymax = 10) #red hake
plot.results.fn(11,6) #cod

