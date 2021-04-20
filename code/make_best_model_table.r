ind = c(1:6,10,7:9)
i = ind[1] #fluke
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + s(length) + 1$|$pair",
  "$\\sim$ dn * s(length) + 1$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic)



i = ind[2] #plaice
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + s(length) + s(length)$|$pair",
  "$\\sim$ dn * s(length) + s(length)$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = rbind(finalout,
    cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic))

i = ind[3] #windowpane
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + length + s(length)$|$pair",
  "$\\sim$ dn * length + s(length)$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = rbind(finalout,
    cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic))

i = ind[4] #winter flounder
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + s(length) + length$|$pair",
  "$\\sim$ dn * s(length) + length$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = rbind(finalout,
    cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic))

i = ind[5] #ytf
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + s(length) + s(length)$|$pair",
  "$\\sim$ dn * s(length) + s(length)$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = rbind(finalout,
    cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic))

i = ind[6] #witch
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + length + s(length)$|$pair",
  "$\\sim$ dn * length + s(length)$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = rbind(finalout,
    cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic))

i = ind[7] #red hake
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + s(length) + s(length)$|$pair",
  "$\\sim$ dn * s(length) + s(length)$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = rbind(finalout,
    cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic))

i = ind[8] #goose
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + s(length) + 1$|$pair",
  "$\\sim$ dn * s(length) + 1$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = rbind(finalout,
    cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic))

i = ind[9] #barndoor
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + length + length$|$pair",
  "$\\sim$ dn * length + length$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = rbind(finalout,
    cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic))

i = ind[10] #thorny
  sp = sp.info$sp.names[i]
  print(sp)
  x = species.get.best[[i]]
  aic = x$aic.converged
  aic = aic[names(aic) %in% first.names]
  first.best = names(sort(aic - min(aic,na.rm=T)))[1] 
  dist = substr(first.best,1,2)
  if(dist == "bi") temp = c(first.best, "bi5","bi6") else temp = c(first.best, "bb8","bb9")
  aic = x$aic.converged[temp]  
  daic = aic - min(aic)
  x = readRDS(paste0(parentdir,"/results/big_results/", sp, "_model_fits.RDS"))$model.fits
  x = x[temp]
  sapply(x, function(y) y$call)
  np = sapply(x, function(y) length(y$model.res$opt$par))
  rho.mod = c(out[temp[1]==first.names,2],
  "$\\sim$ dn + length + length$|$pair",
  "$\\sim$ dn * length + length$|$pair")
  phi.mod = rep(out[temp[1]==first.names,3],3)
  modnum = c(0:4,0:7)[temp[1] == first.names]
  dist = c("BI","BB")[match(dist,c("bi","bb"))]
  models = c(out$Model[temp[1]==first.names], paste0(dist,"$_{", modnum, c("a","b"), "}$"))
  finalout = rbind(finalout,
    cbind.data.frame(sp.info$sp.pretty.names[i], models, rho.mod, phi.mod, np, daic))

colnames(finalout)[2:5] = colnames(out)[1:4]
colnames(finalout)[6] = "$\\Delta$AIC"
finalout[[6]] = as.character(round(finalout[[6]],2))
