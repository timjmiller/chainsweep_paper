trac.biomass.2.fn = function(fall.lendat = all.fall.lendat[[7]], spring.lendat = all.spring.lendat[[8]], dfo.lendat = all.dfo.lendat[[7]],
  fall.lw.dat, spring.lw.dat, dfo.lw.dat, pred.eta.night, pred.eta.day, nefsc.str.size, dfo.str.size, include.fall=TRUE, include.spring = TRUE, include.dfo = TRUE, lengths = 1:50)
{
  #source("get.survey.Nal.hat.fn.r")
  #source("get.dfo.survey.Nal.hat.fn.r")  
  #fall
  out = list(mean.calib.biomass = 0)
  if(include.fall)
  {
    calib.fall.lendat = fall.lendat
    for(i in lengths) 
    {
      ind = which(calib.fall.lendat$LENGTH == i & calib.fall.lendat$day.night == "night")
      calib.fall.lendat$EXPNUMLEN[ind] = exp(pred.eta.night[which(lengths== i)])*calib.fall.lendat$EXPNUMLEN[ind] #change rockhopper night to chainsweep night
      ind = which(calib.fall.lendat$LENGTH == i & calib.fall.lendat$day.night == "day")
      calib.fall.lendat$EXPNUMLEN[ind] = exp(pred.eta.day[which(lengths== i)])*calib.fall.lendat$EXPNUMLEN[ind] #change rockhopper day to chainsweep day
    }
    out$cal.fall.Nal = get.survey.Nal.hat.fn(len.data = calib.fall.lendat, str.size = nefsc.str.size, lengths = lengths)$Nal.hat
    out$fall.Nal = get.survey.Nal.hat.fn(len.data = fall.lendat, str.size = nefsc.str.size, lengths = lengths)$Nal.hat
    dat = list(wal = fall.lw.dat$INDWT, lens = fall.lw.dat$LENGTH, pred_lens = lengths)
    lwmod = MakeADFun(dat, par = list(lw_pars = c(0,0,0)), DLL="wal", silent = TRUE)
    lwmod.opt = nlminb(lwmod$par,lwmod$fn,lwmod$gr)
    out$fall.lw.pred = exp(lwmod$report()$pred_log_wal)
    out$fall.calib.biomass = sum(out$cal.fall.Nal * out$fall.lw.pred)
    out$fall.biomass = sum(out$fall.Nal * out$fall.lw.pred)
    out$mean.calib.biomass = out$mean.calib.biomass + out$fall.calib.biomass
  }  
  
  #spring
  if(include.spring)
  {
    calib.spring.lendat = spring.lendat
    for(i in lengths) 
    {
      ind = which(calib.spring.lendat$LENGTH == i & calib.spring.lendat$day.night == "night")
      calib.spring.lendat$EXPNUMLEN[ind] = exp(pred.eta.night[which(lengths== i)])*calib.spring.lendat$EXPNUMLEN[ind] #change rockhopper night to chainsweep night
      ind = which(calib.spring.lendat$LENGTH == i & calib.spring.lendat$day.night == "day")
      calib.spring.lendat$EXPNUMLEN[ind] = exp(pred.eta.day[which(lengths== i)])*calib.spring.lendat$EXPNUMLEN[ind] #change rockhopper day to chainsweep day
    }
    out$cal.spring.Nal = get.survey.Nal.hat.fn(len.data = calib.spring.lendat, str.size = nefsc.str.size, lengths = lengths)$Nal.hat
    out$spring.Nal = get.survey.Nal.hat.fn(len.data = spring.lendat, str.size = nefsc.str.size, lengths = lengths)$Nal.hat
    dat = list(wal = spring.lw.dat$INDWT, lens = spring.lw.dat$LENGTH, pred_lens = lengths)
    lwmod = MakeADFun(dat, par = list(lw_pars = c(0,0,0)), DLL="wal", silent = TRUE)
    lwmod.opt = nlminb(lwmod$par,lwmod$fn,lwmod$gr)
    out$spring.lw.pred = exp(lwmod$report()$pred_log_wal)
    out$spring.calib.biomass = sum(out$cal.spring.Nal * out$spring.lw.pred)
    out$spring.biomass = sum(out$spring.Nal * out$spring.lw.pred)
    out$mean.calib.biomass = out$mean.calib.biomass + out$spring.calib.biomass
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
    calib.dfo.lendat = dfo.lendat
	  length.cols = which(!(colnames(calib.dfo.lendat) %in% c("STRATA","SLAT","SLONG","UNITAREA","SET","TOTAL","day.night")))
	  lengths.present = as.numeric(substr(colnames(calib.dfo.lendat),2,nchar(colnames(calib.dfo.lendat)))[length.cols])
    for(i in lengths) 
    {
      ind = length.cols[which(lengths.present == i)]
      ind2 = which(calib.dfo.lendat$day.night == 'night')
      if(length(ind)) calib.dfo.lendat[ind2,ind] = exp(pred.eta.night[which(lengths== i)])*calib.dfo.lendat[ind2,ind] #change rockhopper night to chainsweep night
      ind2 = which(calib.dfo.lendat$day.night == 'day')
      if(length(ind)) calib.dfo.lendat[ind2,ind] = exp(pred.eta.day[which(lengths== i)])*calib.dfo.lendat[ind2,ind] #change rockhopper night to chainsweep night
    }
    out$cal.dfo.Nal = get.dfo.survey.Nal.hat.fn(lengths, str.size = dfo.str.size, cal.dat = calib.dfo.lendat)$Nal.hat
    out$dfo.Nal = get.dfo.survey.Nal.hat.fn(lengths, str.size = dfo.str.size, cal.dat = dfo.lendat)$Nal.hat
    dat = list(wal = dfo.lw.dat$INDWT, lens = dfo.lw.dat$LENGTH, pred_lens = lengths)
    lwmod = MakeADFun(dat, par = list(lw_pars = c(0,0,0)), DLL="wal", silent = TRUE)
    lwmod.opt = nlminb(lwmod$par,lwmod$fn,lwmod$gr)
    out$dfo.lw.pred = exp(lwmod$report()$pred_log_wal)
    out$dfo.calib.biomass = sum(out$cal.dfo.Nal * out$dfo.lw.pred)
    out$dfo.biomass = sum(out$dfo.Nal * out$dfo.lw.pred)
    out$mean.calib.biomass = out$mean.calib.biomass + out$dfo.calib.biomass
  }
  
  out$mean.calib.biomass = out$mean.calib.biomass/sum(include.fall,include.spring,include.dfo)
  
  return(out)
  
  
}
