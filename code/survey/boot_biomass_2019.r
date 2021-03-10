
boot_biomass = function(sp.i = 1, stock.i = 1, sp.dat=sp.info, stock.vec=stocks,
  do.boot.biomass.1 = TRUE, do.boot.biomass.2 = TRUE, do.boot.1.res = TRUE, do.boot.2.res = TRUE, 
  years, pred.length.integer = NULL, n.boot = 1000, do.spring = TRUE, do.fall = TRUE, summ.file.mod = "", min.length.calibrate)
{
  svspp = sp.dat$SVSPP[sp.i]
  sp = sp.dat$sp.names[sp.i]
  stock = stock.vec[stock.i]
  fall.str.size = fall.strata.sizes[[stock]]
  spring.str.size = spring.strata.sizes[[stock]]
  source("boot.trac.biomass.1.fn.r")
  source("boot.trac.biomass.2.fn.r")
  source("boot.lendat.fn.r")
  source("get.survey.Nal.hat.fn.r")
  source("boot.dfo.lendat.fn.r")
  source("get.dfo.survey.Nal.hat.fn.r")  
  setwd(sp)
  load("twintrawl.dn.res.RData")
  if(is.null(pred.length.integer)) pred.length.integer = as.integer(read.csv(paste0(summ.file.mod, stock,".spring.Nal.1.csv"))[[1]])
  if(!missing(min.length.calibrate)) pred.length.integer = pred.length.integer[which(pred.length.integer>= min.length.calibrate)]
  #stop()

  library(TMB)
  dyn.load(dynlib("../wal"))

  temp = pred.length.integer
  temp = sapply(temp, function(x) 
  { 
    y = which(twintrawl.dn.res$pred.length == x)
    if(length(y)) return(y)
    else 
    {
      if(x< min(twintrawl.dn.res$pred.length)) return(min(which(twintrawl.dn.res$pred.length %in% temp)))
      else return(max(which(twintrawl.dn.res$pred.length %in% temp)))
    }
  })
  pred.eta.index.for.survey = temp

  ####################################################
  ######### make bootstrap biomass estimates #########
  ####################################################

  load("boot.pred.eta.RData")
  small.boot.eta = boot.pred.eta[,pred.eta.index.for.survey]

  load("boot.pred.eta.night.RData")
  small.boot.eta.night = boot.pred.eta.night[,pred.eta.index.for.survey]

  load("boot.pred.eta.day.RData")
  small.boot.eta.day = boot.pred.eta.day[,pred.eta.index.for.survey]
  
  load(paste0(stock, "_lendat.RData"))
  load(paste0(stock, "_lw_dat.RData"))
  
  if(do.boot.biomass.1)
  {
    for(y in years)
    {
      trac.year = as.character(y)
      #fall.year = as.character(i-1)
      fall.year = as.character(y)
      boot.res = list()
      if(do.fall) boot.res$fall = list()
      if(do.spring) boot.res$spring = list()
      set.seed(060817)
      for(i in 1:n.boot)
      {
      
        if(do.fall) if(fall.year %in% names(all.fall.lendat)) {
          boot.res$fall[[i]] = boot.trac.biomass.1.fn(fall.lendat = all.fall.lendat[[fall.year]], spring.lendat = all.spring.lendat[[trac.year]],
            fall.lw.dat = all.fall.lw.dat[[fall.year]], spring.lw.dat = all.spring.lw.dat[[trac.year]], 
            boot.pred.eta = small.boot.eta[i,], nefsc.str.size = fall.str.size,
            include.fall = TRUE, include.spring = FALSE, include.dfo = FALSE, lengths = pred.length.integer)
        }
        if(do.spring) if(trac.year %in% names(all.spring.lendat)) {
          boot.res$spring[[i]] = boot.trac.biomass.1.fn(fall.lendat = all.fall.lendat[[fall.year]], spring.lendat = all.spring.lendat[[trac.year]],
            fall.lw.dat = all.fall.lw.dat[[fall.year]], spring.lw.dat = all.spring.lw.dat[[trac.year]], 
            boot.pred.eta = small.boot.eta[i,], nefsc.str.size = spring.str.size,
            include.fall = FALSE, include.spring = TRUE, include.dfo = FALSE, lengths = pred.length.integer)
        }
        print(paste0("trac.year: ", trac.year, ", i: ", i, " done"))  
      }
      fn = paste0(stock, "_boot_1_", trac.year)
      if(do.spring) fn = paste0(fn, "_S")
      if(do.fall) fn = paste0(fn, "_F")
      fn = paste0(fn, ".RData")
      save(boot.res, file = fn)
      remove(boot.res)
    }
  }

  if(do.boot.biomass.2)
  {
    for(y in years)
    {
      trac.year = as.character(y)
      fall.year = as.character(y)
      #fall.year = as.character(i-1)
      boot.res = list()
      if(do.fall) boot.res$fall = list()
      if(do.spring) boot.res$spring = list()
      boot.res = list(fall = list(), spring = list())
      set.seed(060816)
      for(i in 1:n.boot)
      {
        if(do.fall) if(fall.year %in% names(all.fall.lendat)) {
          boot.res$fall[[i]] = boot.trac.biomass.2.fn(fall.lendat = all.fall.lendat[[fall.year]], spring.lendat = all.spring.lendat[[trac.year]],
            fall.lw.dat = all.fall.lw.dat[[fall.year]], spring.lw.dat = all.spring.lw.dat[[trac.year]], 
            boot.pred.eta.night = small.boot.eta.night[i,], boot.pred.eta.day = small.boot.eta.day[i,], nefsc.str.size = fall.str.size, 
            include.fall = TRUE, include.spring = FALSE, include.dfo = FALSE, lengths = pred.length.integer)
        }
        if(do.spring)  if(trac.year %in% names(all.spring.lendat)) {
          boot.res$spring[[i]] = boot.trac.biomass.2.fn(fall.lendat = all.fall.lendat[[fall.year]], spring.lendat = all.spring.lendat[[trac.year]],
            fall.lw.dat = all.fall.lw.dat[[fall.year]], spring.lw.dat = all.spring.lw.dat[[trac.year]], 
            boot.pred.eta.night = small.boot.eta.night[i,], boot.pred.eta.day = small.boot.eta.day[i,], nefsc.str.size = spring.str.size, 
            include.fall = FALSE, include.spring = TRUE, include.dfo = FALSE, lengths = pred.length.integer)
        }
        print(paste0("trac.year: ", trac.year, ", i: ", i, " done"))  
      }
      fn = paste0(stock, "_boot_2_", trac.year)
      if(do.spring) fn = paste0(fn, "_S")
      if(do.fall) fn = paste0(fn, "_F")
      fn = paste0(fn, ".RData")
      save(boot.res, file = fn)
      remove(boot.res)
    }
  }


  ####################################################
  ######### generate bootstrap biomass results #######
  ####################################################
  if(do.boot.1.res)
  {
    #need to remove bootstraps when NAs occur for rockhopper efficiency because na.rm in get.survey.Nal.fn does sum(NA,na.rm = TRUE) = 0
    load("boot.pred.eta.RData")
    ind = which(!is.na(boot.pred.eta[,1]))
    boot.1.res = t(sapply(years, function(x)
    {
      f = grep(paste0(stock,"_boot_1_",x, "_"), dir(), value=TRUE)
      load(f)
      out = c()
      if(do.fall & length(grep("_F", f)) & length(boot.res$fall)) {
        out = c(out,sd(sapply(boot.res$fall[ind], function(x) x$fall.calib.biomass))/mean(sapply(boot.res$fall[ind], function(x) x$fall.calib.biomass)))
        out = c(out,quantile(sapply(boot.res$fall[ind], function(x) x$fall.calib.biomass), probs = c(0.025,0.975))/1000)
      } else out = c(out, rep(NA,3))
      if(do.spring & length(grep("_S", f)) & length(boot.res$spring)) {
        out = c(out,sd(sapply(boot.res$spring[ind], function(x) x$spring.calib.biomass))/mean(sapply(boot.res$spring[ind], function(x) x$spring.calib.biomass)))
        out = c(out,quantile(sapply(boot.res$spring[ind], function(x) x$spring.calib.biomass), probs = c(0.025,0.975))/1000)
      } else out = c(out, rep(NA,3))
      return(out)
    }))
    
    boot.1.ratio.res = t(sapply(years, function(x)
    {
      f = grep(paste0(stock,"_boot_1_",x, "_"), dir(), value=TRUE)
      load(f)
      out = c()#
      if(do.fall & length(grep("_F", f)) & length(boot.res$fall)) {
        out = c(out, sd(sapply(boot.res$fall[ind], function(x) x$boot.fall.biomass/x$fall.calib.biomass),na.rm=TRUE)/
          mean(sapply(boot.res$fall[ind], function(x) x$boot.fall.biomass/x$fall.calib.biomass),na.rm=TRUE))
        out = c(out, quantile(sapply(boot.res$fall[ind], function(x) x$boot.fall.biomass/x$fall.calib.biomass), probs = c(0.025,0.975),na.rm=TRUE))
      } else out = c(out, rep(NA,3))
      if(do.spring & length(grep("_S", f)) & length(boot.res$spring)) {
        out = c(out, sd(sapply(boot.res$spring[ind], function(x) x$boot.spring.biomass/x$spring.calib.biomass),na.rm=TRUE)/
          mean(sapply(boot.res$spring[ind], function(x) x$boot.spring.biomass/x$spring.calib.biomass),na.rm=TRUE))
        out = c(out, quantile(sapply(boot.res$spring[ind], function(x) x$boot.spring.biomass/x$spring.calib.biomass), probs = c(0.025,0.975),na.rm=TRUE))
      } else out = c(out, rep(NA,3))
      return(out)
    }))
    colnames(boot.1.res) = colnames(boot.1.res) = c("cv.fall", "lo.fall", "hi.fall", "cv.spring", "lo.spring", "hi.spring")
    if(do.spring | do.fall) {
      rownames(boot.1.ratio.res) = rownames(boot.1.res) = years
      write.csv(round(boot.1.res,2), file = paste0(summ.file.mod, stock, ".cv.biomass.1.csv"))
      write.csv(round(boot.1.ratio.res,2), file = paste0(summ.file.mod, stock, ".cv.biomass.1.ratio.csv"))
    }
  }
  
  if(do.boot.2.res)
  {
    #need to remove bootstraps when NAs occur for rockhopper efficiency because na.rm in get.survey.Nal.fn does sum(NA,na.rm = TRUE) = 0
    load("boot.pred.eta.day.RData")
    load("boot.pred.eta.night.RData")
    ind = which(!is.na(boot.pred.eta.night[,1]) & !is.na(boot.pred.eta.day[,1]))
    boot.2.res = t(sapply(years, function(x)
    {
      f = grep(paste0(stock,"_boot_2_",x, "_"), dir(), value=TRUE)
      load(f)
      out = c()
      if(do.fall & length(grep("_F", f)) & length(boot.res$fall)) {
        out = c(out, sd(sapply(boot.res$fall[ind], function(x) x$fall.calib.biomass))/mean(sapply(boot.res$fall[ind], function(x) x$fall.calib.biomass)))
        out = c(out, quantile(sapply(boot.res$fall[ind], function(x) x$fall.calib.biomass), probs = c(0.025,0.975))/1000)
      } else out = c(out, rep(NA,3))
      if(do.spring & length(grep("_S", f)) & length(boot.res$spring)) {
        out = c(out, sd(sapply(boot.res$spring[ind], function(x) x$spring.calib.biomass))/mean(sapply(boot.res$spring[ind], function(x) x$spring.calib.biomass)))
        out = c(out, quantile(sapply(boot.res$spring[ind], function(x) x$spring.calib.biomass), probs = c(0.025,0.975))/1000)
      } else out = c(out, rep(NA,3))
      return(out)
    }))
    
    boot.2.ratio.res = t(sapply(years, function(x)
    {
      f = grep(paste0(stock,"_boot_2_",x, "_"), dir(), value=TRUE)
      load(f)
      out = c()#
      if(do.fall & length(grep("_F", f)) & length(boot.res$fall)) {
        out = c(out, sd(sapply(boot.res$fall[ind], function(x) x$boot.fall.biomass/x$fall.calib.biomass),na.rm=TRUE)/
          mean(sapply(boot.res$fall[ind], function(x) x$boot.fall.biomass/x$fall.calib.biomass),na.rm=TRUE))
        out = c(out, quantile(sapply(boot.res$fall[ind], function(x) x$boot.fall.biomass/x$fall.calib.biomass), probs = c(0.025,0.975),na.rm=TRUE))
      } else out = c(out, rep(NA,3))
      if(do.spring & length(grep("_S", f)) & length(boot.res$spring)) {
        out = c(out, sd(sapply(boot.res$spring[ind], function(x) x$boot.spring.biomass/x$spring.calib.biomass),na.rm=TRUE)/
          mean(sapply(boot.res$spring[ind], function(x) x$boot.spring.biomass/x$spring.calib.biomass),na.rm=TRUE))
        out = c(out, quantile(sapply(boot.res$spring[ind], function(x) x$boot.spring.biomass/x$spring.calib.biomass), probs = c(0.025,0.975),na.rm=TRUE))
      } else out = c(out, rep(NA,3))
      return(out)
    }))
    colnames(boot.2.res) = colnames(boot.2.res) = c("cv.fall", "lo.fall", "hi.fall", "cv.spring", "lo.spring", "hi.spring")
    if(do.spring | do.fall) {
      rownames(boot.2.ratio.res) = rownames(boot.2.res) = years
      write.csv(round(boot.2.res,2), file = paste0(summ.file.mod, stock, ".cv.biomass.2.csv"))
      write.csv(round(boot.2.ratio.res,2), file = paste0(summ.file.mod, stock, ".cv.biomass.2.ratio.csv"))
    }
  }

  dyn.unload(dynlib("../wal"))
  setwd(parentdir)
  return(pred.length.integer)
}
