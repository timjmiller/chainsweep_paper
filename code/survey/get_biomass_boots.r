get_biomass_boots = function(sp.i = 1, stock.i = 1, sp.dat=sp.info, stock.vec=stocks, years)
{
  svspp = sp.dat$SVSPP[sp.i]
  stock = stock.vec[stock.i]
  fall.str.size = fall.strata.sizes[[stock]]
  spring.str.size = spring.strata.sizes[[stock]]

  boot.pred.eta = readRDS(paste0("results/", sp.dat$sp.names[sp.i], "_boot_pred_eta_0.RDS"))
  #temp = boot.pred.eta
  for(j in 1:9) boot.pred.eta = rbind(boot.pred.eta, readRDS(paste0("results/", sp.dat$sp.names[sp.i], "_boot_pred_eta_",j,".RDS")))
  nsim = NROW(boot.pred.eta)
  ind = which(!is.na(boot.pred.eta[,1]))
  nsim = length(ind)
  print(1)
  fall.calib.boots = t(sapply(years, function(x)
  {
    paste0(stock,"_boot_",x, "_")
    f = grep(paste0(stock,"_boot_",x, "_"), dir("results"), value=TRUE)
    print(f)
    boot.res = readRDS(paste0("results/",f))
    print(names(boot.res))
    out = rep(NA,nsim)
    if(length(grep("_F", f)) & length(boot.res$fall)) {
      out = sapply(boot.res$fall[ind], function(x) x$fall.calib.biomass)
    }
    return(out)
  }))
  spring.calib.boots = t(sapply(years, function(x)
  {
    paste0(stock,"_boot_",x, "_")
    f = grep(paste0(stock,"_boot_",x, "_"), dir("results"), value=TRUE)
    print(f)
    boot.res = readRDS(paste0("results/",f))
    print(names(boot.res))
    out = rep(NA,nsim)
    if(length(grep("_S", f)) & length(boot.res$spring)) {
      out = sapply(boot.res$spring[ind], function(x) x$spring.calib.biomass)
    } 
    return(out)
  }))
  fall.uncalib.boots = t(sapply(years, function(x)
  {
    paste0(stock,"_boot_",x, "_")
    f = grep(paste0(stock,"_boot_",x, "_"), dir("results"), value=TRUE)
    print(f)
    boot.res = readRDS(paste0("results/",f))
    print(names(boot.res))
    out = rep(NA,nsim)
    if(length(grep("_F", f)) & length(boot.res$fall)) {
      out = sapply(boot.res$fall[ind], function(x) x$boot.fall.biomass)
    }
    return(out)
  }))
  spring.uncalib.boots = t(sapply(years, function(x)
  {
    paste0(stock,"_boot_",x, "_")
    f = grep(paste0(stock,"_boot_",x, "_"), dir("results"), value=TRUE)
    print(f)
    boot.res = readRDS(paste0("results/",f))
    print(names(boot.res))
    out = rep(NA,nsim)
    if(length(grep("_S", f)) & length(boot.res$spring)) {
      out = sapply(boot.res$spring[ind], function(x) x$boot.spring.biomass)
    } 
    return(out)
  }))
  out = list(fall.calib.boots, fall.uncalib.boots, spring.calib.boots,spring.uncalib.boots)
  names(out) = c("fall.calib.boots", "fall.uncalib.boots", "spring.calib.boots","spring.uncalib.boots")
  return(out)
}