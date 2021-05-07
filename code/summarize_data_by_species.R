# This will give the numbers in the text of working papers that summarinze the amounts of data available in various ways for each species. 
library(Hmisc)
x = read.csv("spp_list.csv")
x = x[which(x$NSTATIONS>30 & !(x$COMNAME %in% c("FOURSPOT FLOUNDER", "GULF STREAM FLOUNDER"))),]
spps = c(103, 102, 108, 106, 105, 107, 197, 22, 28, 77, 73)
x = x[match(spps,x$SVSPP),]
x$NSTOCKS = c(1,1,2,3,3,1,2,1,1,2,2)
x$sp.names = c("fluke", "plaice", "windowpane", "winter_flounder", "yellowtail_flounder", "witch_flounder", "goosefish", "barndoor", "thorny", "red_hake", "cod")
x$sp.pretty.names = capitalize(tolower(x$COMNAME))
sp.info = x
parentdir = getwd()


x = matrix(NA, NROW(sp.info), 13)
colnames(x) = c("n.sta", "n.sta.d", "n.sta.n", "n.len", "n.len.d", "n.len.n", "exp.n.len", "chain.n.len", "rock.n.len", "chain.n.len.d", "chain.n.len.n", "rock.n.len.d", "rock.n.len.n")
for(i in 1:NROW(sp.info))
{
  sp = sp.info$sp.names[i]
  spp = sp.info$SVSPP[i]
  combined.data = readRDS(paste0(parentdir,"/data/",sp,"_data.RDS"))

  n.sta = aggregate(id ~ 1, data = combined.data$catch.length, FUN = function(x) length(unique(x)))
  
  n.sta.dn = aggregate(id ~ day.night, data = combined.data$catch.length, FUN = function(x) length(unique(x)))

  n.len = aggregate(recnumlen.ch + recnumlen.rh ~ 1, data = combined.data$catch.length, FUN = sum)
  n.len.dn = aggregate(recnumlen.ch + recnumlen.rh ~ day.night, data = combined.data$catch.length, FUN = sum)

  expn.len = aggregate(expnumlen.ch + expnumlen.rh ~ 1, data = combined.data$catch.length, FUN = sum)

  chain.n.len = aggregate(recnumlen.ch ~ 1, data = combined.data$catch.length, FUN = sum)
  rock.n.len = aggregate(recnumlen.rh ~ 1, data = combined.data$catch.length, FUN = sum)

  chain.n.len.dn = aggregate(recnumlen.ch ~ day.night, data = combined.data$catch.length, FUN = sum)
  rock.n.len.dn = aggregate(recnumlen.rh ~ day.night, data = combined.data$catch.length, FUN = sum)
  temp = unlist(c(n.sta, n.sta.dn, n.len, n.len.dn, expn.len, chain.n.len, rock.n.len, chain.n.len.dn, rock.n.len.dn))
  print(length(temp))
  print(temp)
  x[i,] = unlist(c(n.sta, n.sta.dn[,2], n.len, n.len.dn[,2], expn.len, chain.n.len, rock.n.len, chain.n.len.dn[,2], rock.n.len.dn[,2]))
}

#make some plots of this (length frequencies
for(i in 1:NROW(sp.info))
{
  sp = sp.info$sp.names[i]
  setwd(sp)
  load("combined.dn.data.for.boot.RData")
  dat = combined.dn.data$catch.length
  
  x = aggregate(cbind(expnumlen.rh,expnumlen.ch) ~ length * day.night, data = combined.dn.data$catch.length, FUN = sum) 
  colnames(x)[3:4] = c("expnumlen.rockhopper", "expnumlen.chainsweep")
  setwd(parentdir)
  write.csv(x, file = paste0(parentdir,"/results/", sp, "/twin_trawl_length_frequencies.csv"))
}

#make summary table of AIC for each model/species
setwd(parentdir)
source(paste0(parentdir,"/code/get.best.r")
source(paste0(parentdir,"/code/make.all.AIC.table.fn.r"))
aic.table = make.all.AIC.table.fn(sp.info$sp.names, sp.info$sp.pretty.names, latex.table.file = paste0(parentdir,"/paper/model_compare.tex"))
models = paste0(rep(c("BI_","BB_"), c(5,8)), c(0:4, 0:7))
rownames(aic.table) = c(models, paste0(models, "_DN"))
colnames(aic.table)[1] = "No. parameters"
write.csv(aic.table, file = paste0(parentdir,"/paper/model_compare.csv"))


#done to here

setwd(parentdir)
source("make.nobs.table.fn.r")
nobs.table = make.nobs.table.fn()
