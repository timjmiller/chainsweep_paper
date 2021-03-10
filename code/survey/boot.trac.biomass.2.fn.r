boot.trac.biomass.2.fn = function(fall.lendat = all.fall.lendat[[7]], spring.lendat = all.spring.lendat[[8]], dfo.lendat = all.dfo.lendat[[7]],
  fall.lw.dat, spring.lw.dat, dfo.lw.dat, boot.pred.eta.night, boot.pred.eta.day, nefsc.str.size, dfo.str.size, include.fall=TRUE, include.spring = TRUE, include.dfo = TRUE, fit.models = TRUE, lengths = 1:50)
{
  #source("boot.lendat.fn.r")
  #source("get.survey.Nal.hat.fn.r")
  #source("boot.dfo.lendat.fn.r")
  #source("get.dfo.survey.Nal.hat.fn.r")  
  #fall
  out = list(mean.calib.biomass = 0)
  if(include.fall)
  {
    boot.fall.lendat = boot.lendat.fn(fall.lendat)
    boot.fall.lw.dat = fall.lw.dat[sample(1:NROW(fall.lw.dat), size = NROW(fall.lw.dat), replace = TRUE),]
    if(fit.models)
    {
      out$boot.fall.Nal = get.survey.Nal.hat.fn(len.data = boot.fall.lendat, str.size = nefsc.str.size, lengths = lengths)$Nal.hat
      calib.boot.fall.lendat = boot.fall.lendat
      for(i in lengths) 
      {
        ind = which(calib.boot.fall.lendat$LENGTH == i & calib.boot.fall.lendat$day.night == "night")
        calib.boot.fall.lendat$EXPNUMLEN[ind] = exp(boot.pred.eta.night[which(lengths== i)])*calib.boot.fall.lendat$EXPNUMLEN[ind] #change rockhopper night to chainsweep night
        ind = which(calib.boot.fall.lendat$LENGTH == i & calib.boot.fall.lendat$day.night == "day")
        calib.boot.fall.lendat$EXPNUMLEN[ind] = exp(boot.pred.eta.day[which(lengths== i)])*calib.boot.fall.lendat$EXPNUMLEN[ind] #change rockhopper day to chainsweep day
      }
      out$cal.boot.fall.Nal = get.survey.Nal.hat.fn(len.data = calib.boot.fall.lendat, str.size = nefsc.str.size, lengths = lengths)$Nal.hat
      dat = list(wal = boot.fall.lw.dat$INDWT, lens = boot.fall.lw.dat$LENGTH, pred_lens = lengths)
      lwmod = MakeADFun(dat, par = list(lw_pars = c(0,0,0)), DLL="wal", silent = TRUE)
      lwmod.opt = nlminb(lwmod$par,lwmod$fn,lwmod$gr)
      out$boot.fall.lw.pred = exp(lwmod$report()$pred_log_wal)
      out$fall.calib.biomass = sum(out$cal.boot.fall.Nal * out$boot.fall.lw.pred)
      out$boot.fall.biomass = sum(out$boot.fall.Nal * out$boot.fall.lw.pred)
      out$mean.calib.biomass = out$mean.calib.biomass + out$fall.calib.biomass
    }
  }  
  
  #spring
  if(include.spring)
  {
    boot.spring.lendat = boot.lendat.fn(spring.lendat)
    boot.spring.lw.dat = spring.lw.dat[sample(1:NROW(spring.lw.dat), size = NROW(spring.lw.dat), replace = TRUE),]
    if(fit.models)
    {
      out$boot.spring.Nal = get.survey.Nal.hat.fn(len.data = boot.spring.lendat, str.size = nefsc.str.size, lengths = lengths)$Nal.hat
      calib.boot.spring.lendat = boot.spring.lendat
      for(i in lengths) 
      {
        ind = which(calib.boot.spring.lendat$LENGTH == i & calib.boot.spring.lendat$day.night == "night")
        calib.boot.spring.lendat$EXPNUMLEN[ind] = exp(boot.pred.eta.night[which(lengths== i)])*calib.boot.spring.lendat$EXPNUMLEN[ind] #change rockhopper night to chainsweep night
        ind = which(calib.boot.spring.lendat$LENGTH == i & calib.boot.spring.lendat$day.night == "day")
        calib.boot.spring.lendat$EXPNUMLEN[ind] = exp(boot.pred.eta.day[which(lengths== i)])*calib.boot.spring.lendat$EXPNUMLEN[ind] #change rockhopper day to chainsweep day
      }
      out$cal.boot.spring.Nal = get.survey.Nal.hat.fn(len.data = calib.boot.spring.lendat, str.size = nefsc.str.size, lengths = lengths)$Nal.hat
      dat = list(wal = boot.spring.lw.dat$INDWT, lens = boot.spring.lw.dat$LENGTH, pred_lens = lengths)
      lwmod = MakeADFun(dat, par = list(lw_pars = c(0,0,0)), DLL="wal", silent = TRUE)
      lwmod.opt = nlminb(lwmod$par,lwmod$fn,lwmod$gr)
      out$boot.spring.lw.pred = exp(lwmod$report()$pred_log_wal)
      out$spring.calib.biomass = sum(out$cal.boot.spring.Nal * out$boot.spring.lw.pred)
      out$boot.spring.biomass = sum(out$boot.spring.Nal * out$boot.spring.lw.pred)
      out$mean.calib.biomass = out$mean.calib.biomass + out$spring.calib.biomass
    }
  }  
  
  get.dfo.moddat.fn = function(cal.dat)
  {
    length.cols = c(6:(NCOL(cal.dat)-2))
    lengths = as.integer(substr(colnames(cal.dat)[length.cols],2,nchar(colnames(cal.dat)[length.cols])))
    ids = cal.dat$SET
    nal.tow <- sapply(lengths, function(x) sapply(ids, function(y) {
      if(x %in% lengths) return(cal.dat[[which(lengths == x)+5]][cal.dat$SET == y])
      else return(0)
    }))
    x = cbind.data.frame(count = as.vector(t(nal.tow[,1:length(lengths)])), length = rep(lengths, NROW(nal.tow)), dn = rep(cal.dat$day.night, each = length(lengths)),
      set = rep(cal.dat$SET, each = length(lengths)), stratum = factor(rep(cal.dat$STRATA, each = length(lengths))))
    return(x)
  }
  if(include.dfo)
  {
    boot.dfo.lendat = boot.dfo.lendat.fn(dfo.lendat)
    boot.dfo.lw.dat = dfo.lw.dat[sample(1:NROW(dfo.lw.dat), size = NROW(dfo.lw.dat), replace = TRUE),]
    #must do this twice to make sure original results are repeatable with a given random seed. 
    boot.dfo.lw.dat = dfo.lw.dat[sample(1:NROW(dfo.lw.dat), size = NROW(dfo.lw.dat), replace = TRUE),]
    if(fit.models)
    {
      out$boot.dfo.Nal = get.dfo.survey.Nal.hat.fn(lengths, str.size = dfo.str.size, cal.dat = boot.dfo.lendat)$Nal.hat
      calib.boot.dfo.lendat = boot.dfo.lendat
	    length.cols = which(!(colnames(calib.boot.dfo.lendat) %in% c("STRATA","SLAT","SLONG","UNITAREA","SET","TOTAL","day.night")))
	    lengths.present = as.numeric(substr(colnames(calib.boot.dfo.lendat),2,nchar(colnames(calib.boot.dfo.lendat)))[length.cols])
      for(i in lengths) 
      {
        ind = length.cols[which(lengths.present == i)]
        ind2 = which(calib.boot.dfo.lendat$day.night == 'night')
        if(length(ind)) calib.boot.dfo.lendat[ind2,ind] = exp(boot.pred.eta.night[which(lengths== i)])*calib.boot.dfo.lendat[ind2,ind] #change rockhopper night to chainsweep night
        ind2 = which(calib.boot.dfo.lendat$day.night == 'day')
        if(length(ind)) calib.boot.dfo.lendat[ind2,ind] = exp(boot.pred.eta.day[which(lengths== i)])*calib.boot.dfo.lendat[ind2,ind] #change rockhopper night to chainsweep night
      }
      out$cal.boot.dfo.Nal = get.dfo.survey.Nal.hat.fn(lengths, str.size = dfo.str.size, cal.dat = calib.boot.dfo.lendat)$Nal.hat
      dat = list(wal = boot.dfo.lw.dat$INDWT, lens = boot.dfo.lw.dat$LENGTH, pred_lens = lengths)
      lwmod = MakeADFun(dat, par = list(lw_pars = c(0,0,0)), DLL="wal", silent = TRUE)
      lwmod.opt = nlminb(lwmod$par,lwmod$fn,lwmod$gr)
      out$boot.dfo.lw.pred = exp(lwmod$report()$pred_log_wal)
      out$dfo.calib.biomass = sum(out$cal.boot.dfo.Nal * out$boot.dfo.lw.pred)
      out$boot.dfo.biomass = sum(out$boot.dfo.Nal * out$boot.dfo.lw.pred)
      out$mean.calib.biomass = out$mean.calib.biomass + out$dfo.calib.biomass
    }
  }
  
  out$mean.calib.biomass = out$mean.calib.biomass/sum(include.fall,include.spring,include.dfo)
  
  return(out)
  
  
}
