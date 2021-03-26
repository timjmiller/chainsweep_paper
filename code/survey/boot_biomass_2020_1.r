
boot_biomass = function(sp.i = 1, stock.i = 1, sp.dat=sp.info, stock.vec=stocks,
  do.boot.biomass = TRUE, do.boot.res = TRUE, 
  years, pred.length.integer = NULL, n.boot = 1000, do.spring = TRUE, do.fall = TRUE, summ.file.mod = "", min.length.calibrate, max.length.calibrate)
{
  #added in max.length.calibrate argument
  svspp = sp.dat$SVSPP[sp.i]
  stock = stock.vec[stock.i]
  fall.str.size = fall.strata.sizes[[stock]]
  spring.str.size = spring.strata.sizes[[stock]]
  
  twintrawl.res = readRDS(paste0("results/big_results/", sp.dat$sp.names[sp.i], "_model_fits.RDS"))
  load(paste0("data/survey/", stock, "_lendat.RData"))
  load(paste0("data/survey/", stock, "_lw_dat.RData"))
  boot.pred.eta = readRDS("boot.pred.eta.RDS")
  
  if(is.null(pred.length.integer)) pred.length.integer = as.integer(read.csv(paste0(summ.file.mod, stock,".spring.Nal.csv"))[[1]])

  if(!missing(min.length.calibrate)) pred.length.integer = pred.length.integer[which(pred.length.integer>= min.length.calibrate)]
  if(!missing(max.length.calibrate)) pred.length.integer = pred.length.integer[which(pred.length.integer<= max.length.calibrate)]
  if(!length(pred.length.integer)) stop("pred.length.integer does not work with supplied min.length.calibrate and/or max.length.calibrate")
  #stop()

  temp = pred.length.integer
  temp = sapply(temp, function(x) 
  { 
    y = which(twintrawl.res$pred.length == x)
    if(length(y)) return(y)
    else 
    {
      if(x< min(twintrawl.res$pred.length)) return(min(which(twintrawl.res$pred.length %in% temp)))
      else return(max(which(twintrawl.res$pred.length %in% temp)))
    }
  })
  pred.eta.index.for.survey = temp

  ####################################################
  ######### make bootstrap biomass estimates #########
  ####################################################

  
  small.boot.eta.day = boot.pred.eta[,pred.eta.index.for.survey]
  if(NCOL(boot.pred.eta)> length(twintrawl.res$pred.length)) #day first, night second
  {
    small.boot.eta.night = boot.pred.eta[,length(twintrawl.res$pred.length) + pred.eta.index.for.survey]
  }
  else small.boot.eta.night = boot.pred.eta.day
    
  
  if(do.boot.biomass)
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
          boot.res$fall[[i]] = boot.trac.biomass.dn.fn(fall.lendat = all.fall.lendat[[fall.year]], spring.lendat = all.spring.lendat[[trac.year]],
            fall.lw.dat = all.fall.lw.dat[[fall.year]], spring.lw.dat = all.spring.lw.dat[[trac.year]], 
            boot.pred.eta.night = small.boot.eta.night[i,], boot.pred.eta.day = small.boot.eta.day[i,], nefsc.str.size = fall.str.size, 
            include.fall = TRUE, include.spring = FALSE, include.dfo = FALSE, lengths = pred.length.integer)
        }
        if(do.spring)  if(trac.year %in% names(all.spring.lendat)) {
          boot.res$spring[[i]] = boot.trac.biomass.dn.fn(fall.lendat = all.fall.lendat[[fall.year]], spring.lendat = all.spring.lendat[[trac.year]],
            fall.lw.dat = all.fall.lw.dat[[fall.year]], spring.lw.dat = all.spring.lw.dat[[trac.year]], 
            boot.pred.eta.night = small.boot.eta.night[i,], boot.pred.eta.day = small.boot.eta.day[i,], nefsc.str.size = spring.str.size, 
            include.fall = FALSE, include.spring = TRUE, include.dfo = FALSE, lengths = pred.length.integer)
        }
        print(paste0("trac.year: ", trac.year, ", i: ", i, " done"))  
      }
      fn = paste0("results/", stock, "_boot_", trac.year)
      if(do.spring) fn = paste0(fn, "_S")
      if(do.fall) fn = paste0(fn, "_F")
      fn = paste0(fn, ".RDS")
      saveRDS(boot.res, file = fn)
      remove(boot.res)
    }
  }


  ####################################################
  ######### generate bootstrap biomass results #######
  ####################################################  
  if(do.boot.res)
  {
    #need to remove bootstraps when NAs occur for rockhopper efficiency because na.rm in get.survey.Nal.fn does sum(NA,na.rm = TRUE) = 0
    ind = which(!is.na(boot.pred.eta[,1]))
    boot.res = t(sapply(years, function(x)
    {
      f = grep(paste0("results/", stock,"_boot_",x, "_"), dir(), value=TRUE)
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
    
    boot.ratio.res = t(sapply(years, function(x)
    {
      f = grep(paste0(stock,"_boot_",x, "_"), dir(), value=TRUE)
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
    colnames(boot.res) = colnames(boot.res) = c("cv.fall", "lo.fall", "hi.fall", "cv.spring", "lo.spring", "hi.spring")
    if(do.spring | do.fall) {
      rownames(boot.ratio.res) = rownames(boot.res) = years
      write.csv(round(boot.res,2), file = paste0("results/", stock, ".cv.biomass.csv"))
      write.csv(round(boot.ratio.res,2), file = paste0("results/", stock, ".cv.biomass.ratio.csv"))
    }
  }

  #dyn.unload(dynlib("../wal"))
  print(parentdir)
  #setwd(parentdir)
  #return(pred.length.integer)
}
