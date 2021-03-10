
estimate_biomass = function(sp.i = 1, stock.i = 1, sp.info, sole, do.lendat = TRUE, do.lw.dat = TRUE, do.N.W = TRUE, years = NULL, pred.length.integer = NULL,
min.length.calibrate, max.length.calibrate, filename.mod = "", TOGA.type = 1, TOGA.operation = 3, TOGA.gear = 2)
{
  
  #added in max.length.calibrate argument

  ####################################################
  ############### get survey data ####################
  ####################################################
  svspp = sp.info$SVSPP[sp.i]
  sp = sp.info$sp.names[sp.i]
  stock = stocks[stock.i]
  x = get.best(sp)
  best.model = x$best.model #"bi3"
  setwd(sp)
 
  #fluke
  library(TMB)
  dyn.load(dynlib("../wal"))

  get.survey.data.fn <- function(oc = sole, purpose.code = 10){
    q.surveys <- "select cruise6,purpose_code,season,year, startdate, enddate from mstr_cruise where year > 2008 and purpose_code= 10 and cruise6 is not null"
    survey.view <- sqlQuery(oc,q.surveys)
    return(survey.view)
  }
  survey.cruises <- get.survey.data.fn()
  fall.survey.cruises <- survey.cruises[survey.cruises$SEASON == 'FALL',]
  fall.surveys <- sort(unique(fall.survey.cruises$CRUISE6))
  
  spring.survey.cruises <- survey.cruises[survey.cruises$SEASON == 'SPRING',]
  spring.surveys <- sort(unique(spring.survey.cruises$CRUISE6))
  if(is.null(years)) years = as.numeric(substr(fall.surveys,1,4))
  else 
  {
    fall.surveys = fall.surveys[which(as.numeric(substr(fall.surveys,1,4)) %in% years)]
    spring.surveys = spring.surveys[which(as.numeric(substr(spring.surveys,1,4)) %in% years)]
  }
  twintrawl.dn.res = readRDS("model_fits.RDS")
  if(is.null(pred.length.integer)) pred.length.integer = floor(min(twintrawl.dn.res$pred.length)):ceiling(max(twintrawl.dn.res$pred.length))
  source("../sunPosition.fn.r")
  source("../get.lendat.dn.TOGA.fn.r")
  source("../get.survey.stratum.estimates.TOGA.fn.r")
  source("../get.survey.Nal.hat.fn.r")
  source("../get.dfo.survey.Nal.hat.fn.r")
  source("../trac.biomass.1.fn.r")
  source("../trac.biomass.2.fn.r")

  get.length.range.fn = function(survey,svspp=105, lens = 1:50, strata = svspp.strata)
  {
    x = try(get.survey.stratum.estimates.TOGA.fn(spp = svspp, survey = survey, strata = strata, do.length = TRUE, lengths = lens, do.age = FALSE,
    typecode = TOGA.type, operationcode = TOGA.operation, gearcode = TOGA.gear))
    if(!is.character(x)) lengths = x$lengths
    else lengths = NA
    return(lengths)
  }
  fall.lengths = lapply(fall.surveys, get.length.range.fn, svspp = svspp, lens = pred.length.integer, strata = fall.strata[[stock]])
  no.fall.data = sapply(fall.lengths, function(x) length(x) == 1)
  if(any(no.fall.data)) {
    fall.lengths = fall.lengths[which(!no.fall.data)]
    fall.surveys = fall.surveys[which(!no.fall.data)]
  }
  spring.lengths = lapply(spring.surveys, get.length.range.fn, svspp = svspp, lens = pred.length.integer, strata = spring.strata[[stock]])
  no.spring.data = sapply(fall.lengths, function(x) length(x) == 1)
  if(any(no.spring.data)) {
    spring.lengths = spring.lengths[which(!no.spring.data)]
    spring.surveys = spring.surveys[which(!no.spring.data)]
  }

  min.lengths = min(sapply(fall.lengths, min),sapply(spring.lengths, min))
  max.lengths = max(sapply(fall.lengths, max),sapply(spring.lengths, max))
  if(min(pred.length.integer)>min.lengths || max(pred.length.integer) < max.lengths)
  {
    warning("adjusting pred.length.integer to deal with ranges of survey data lengths")
    if(min(pred.length.integer)>min.lengths) pred.length.integer = sort(c(min.lengths:(min(pred.length.integer)-1), pred.length.integer))
    if(max(pred.length.integer) < max.lengths) pred.length.integer = sort(c(pred.length.integer, (max(pred.length.integer)+1):max.lengths))
  }

  if(do.lendat)
  {
    all.fall.lendat = lapply(fall.surveys, get.lendat.dn.TOGA.fn, svspp = svspp, lens = pred.length.integer, strata = fall.strata[[stock]], 
		typec = TOGA.type, operc = TOGA.operation, gearc = TOGA.gear)
    all.spring.lendat = lapply(spring.surveys, get.lendat.dn.TOGA.fn, svspp = svspp, lens = pred.length.integer, strata = spring.strata[[stock]],
		typec = TOGA.type, operc = TOGA.operation, gearc = TOGA.gear)
    names(all.fall.lendat) = as.numeric(substr(fall.surveys,1,4))
    names(all.spring.lendat) = as.numeric(substr(spring.surveys,1,4))
    save(all.fall.lendat,all.spring.lendat, file = paste0(stock, "_lendat.RData"))
  }
  else load(paste0(stock, "_lendat.RData"))
  fall.str.size = get.survey.stratum.estimates.TOGA.fn(spp = svspp, survey = fall.surveys[1], lengths = pred.length.integer,tow_swept_area = 0.007, 
	strata = fall.strata[[stock]])$stratum.size
  spring.str.size = get.survey.stratum.estimates.TOGA.fn(spp = svspp, survey = spring.surveys[1], lengths = pred.length.integer,tow_swept_area = 0.007, 
	strata = spring.strata[[stock]])$stratum.size

  get.lw.data.fn = function(survey,lendat=NULL,svspp=105, strata = svspp.strata, oc = sole)
  {
    q.lw <- paste("select  cruise6, stratum, tow, station, sex, length, age, indwt, maturity from union_fscs_svbio ",
      "where cruise6 = ", survey, " and STRATUM IN('", paste(strata, collapse = "','"), "')",
      "and svspp = ", svspp, " order by cruise6, stratum, tow, station", sep = '')
    lw.view <- sqlQuery(oc,q.lw,stringsAsFactors= FALSE)
    return(lw.view)
  }
  
  if(do.lw.dat)
  {
    all.fall.lw.dat = lapply(fall.surveys, function(x)
    {
      lw.dat <- get.lw.data.fn(svspp = svspp, survey = x, strata = fall.strata[[stock]])
      lw.dat = na.omit(lw.dat[c("LENGTH","INDWT")])
      lw.dat = lw.dat[which(lw.dat$INDWT>0 & lw.dat$LENGTH>0),]
    })
    all.spring.lw.dat = lapply(spring.surveys, function(x)
    {
      lw.dat <- get.lw.data.fn(svspp = svspp, survey = x, strata = spring.strata[[stock]])
      lw.dat = na.omit(lw.dat[c("LENGTH","INDWT")])
      lw.dat = lw.dat[which(lw.dat$INDWT>0 & lw.dat$LENGTH>0),]
    })
    names(all.fall.lw.dat) = as.numeric(substr(fall.surveys,1,4))
    names(all.spring.lw.dat) = as.numeric(substr(spring.surveys,1,4))
    save(all.spring.lw.dat,all.fall.lw.dat, file = paste0(stock, "_lw_dat.RData"))
  }
  else load(paste0(stock, "_lw_dat.RData"))

  if(do.N.W)
  {
    all.fall.N.W = sapply(fall.surveys, function(x)
    {
      print(x)
      dat <- get.survey.stratum.estimates.TOGA.fn(spp = svspp, survey = x, do.length = TRUE, lengths = pred.length.integer, do.age = TRUE, 
		strata = fall.strata[[stock]], typecode = TOGA.type, operationcode = TOGA.operation, gearcode = TOGA.gear)
      if(any(dat$out[,3] ==0)) warning(paste0("Some strata for survey ", x, " have zero tows, will extrapolate data to these areas"))
      N = sum(dat$out[,2] * dat$out[,4], na.rm = TRUE) * sum(dat$out[,2])/sum(dat$out[which(dat$out[,3]>0),2])
      N.pertow = N/ sum(dat$out[,2])
      W = sum(dat$out[,2] * dat$out[,5], na.rm = TRUE) * sum(dat$out[,2])/sum(dat$out[which(dat$out[,3]>0),2])
      W.pertow = W/sum(dat$out[,2])
      return(c(N = N, W = W, n = N.pertow, w = W.pertow, min.len = min(dat$lengths), max.len = max(dat$lengths)))
    })
    #some strata each year have zero tows for fall survey

    all.spring.N.W = sapply(spring.surveys, function(x)
    {
      print(x)
      dat <- get.survey.stratum.estimates.TOGA.fn(spp = svspp, survey = x, do.length = TRUE, lengths = pred.length.integer, do.age = TRUE, 
      strata = spring.strata[[stock]], typecode = TOGA.type, operationcode = TOGA.operation, gearcode = TOGA.gear)
      if(any(dat$out[,3] ==0)) warning(paste0("Some strata for survey ", x, " have zero tows, will extrapolate data to these areas"))
      N = sum(dat$out[,2] * dat$out[,4], na.rm = TRUE) * sum(dat$out[,2])/sum(dat$out[which(dat$out[,3]>0),2])
      N.pertow = N/ sum(dat$out[,2])
      W = sum(dat$out[,2] * dat$out[,5], na.rm = TRUE) * sum(dat$out[,2])/sum(dat$out[which(dat$out[,3]>0),2])
      W.pertow = W/sum(dat$out[,2])
      return(c(N = N, W = W, n = N.pertow, w = W.pertow, min.len = min(dat$lengths), max.len = max(dat$lengths)))
    })
    colnames(all.spring.N.W) = as.numeric(substr(spring.surveys,1,4))
    colnames(all.fall.N.W) = as.numeric(substr(fall.surveys,1,4))
    
    save(all.fall.N.W, all.spring.N.W, file = paste0(stock, "_N.W.RData"))
  }
  else load(paste0(stock, "_N.W.RData"))
  
  calibrate.length.integer = pred.length.integer
  if(!missing(min.length.calibrate)) calibrate.length.integer = calibrate.length.integer[which(calibrate.length.integer>= min.length.calibrate)]
  if(!missing(max.length.calibrate)) calibrate.length.integer = calibrate.length.integer[which(calibrate.length.integer<= max.length.calibrate)]
  if(!length(calibrate.length.integer)) stop("pred.length.integer does not work with supplied min.length.calibrate and/or max.length.calibrate")
  #print("here")
  #temp = pred.length.integer
  temp = calibrate.length.integer
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
  #print("here")

  #scaling method 2 (if no day/night effects just use the same mean_pred_eta
  sdrep = summary(twintrawl.dn.res$model.fits[[best.model]]$model.res$sdrep)
  temp = list()
  temp[[1]] = temp[[2]] = sdrep[which(rownames(sdrep) == "mean_pred_eta"),1]
  if(length(temp[[1]] > length(twintrawl.dn.res$pred.length))) # mean_pred_eta should be 1 or 2 times length(pred.length)
  {
    temp[[2]] = temp[[1]][length(twintrawl.dn.res$pred.length) + 1:length(twintrawl.dn.res$pred.length)] #night
    temp[[1]] = temp[[1]][1:length(twintrawl.dn.res$pred.length)] #day
  }
  #print("here")
  temp = cbind(temp[[1]],temp[[2]]) #day,night
  #print(exp(head(temp)))
  temp = temp[pred.eta.index.for.survey,]
  biomass.estimates = list(fall = list(), spring = list())
  for(y in names(all.fall.lendat))
  {
    biomass.estimates$fall[[y]] = trac.biomass.2.fn(fall.lendat = all.fall.lendat[[y]], 
      fall.lw.dat = all.fall.lw.dat[[y]], pred.eta.night = temp[,2], pred.eta.day = temp[,1],
      nefsc.str.size = fall.str.size, include.fall = TRUE, include.spring = FALSE, include.dfo = FALSE, lengths = calibrate.length.integer)
  }
  #print("here")
  for(y in names(all.spring.lendat))
  {
    biomass.estimates$spring[[y]] = trac.biomass.2.fn(spring.lendat = all.spring.lendat[[y]], 
      spring.lw.dat = all.spring.lw.dat[[y]], pred.eta.night = temp[,2], pred.eta.day = temp[,1],
      nefsc.str.size = spring.str.size, include.fall = FALSE, include.spring = TRUE, include.dfo = FALSE, lengths = calibrate.length.integer)
    #print(biomass.estimates$spring[[y]])
  }
  #print("here")
  
  temp.a.fn = function(yrs=years,var= "fall.calib.biomass",biomass = biomass.estimates$fall)
  {
    sapply(as.character(yrs), function(y) 
    ifelse(y %in% names(biomass), biomass[[y]][[var]], NA))
  }
  temp.b.fn = function(yrs=years,var= "fall.Nal",biomass = biomass.estimates$fall)
  {
    sapply(as.character(yrs), function(y)
    {
      #print(y)
      ind =  y %in% names(biomass)
      #print(var)
      #print(names(biomass[[y]]))
      #print(biomass[[y]][[var]])
      if(ind) return(biomass[[y]][[var]])
      else return(rep(NA,length(pred.length.integer)))
    })
  }
  #print(names(biomass.estimates$fall))
  #print(names(biomass.estimates$spring))
  
  biomass.res = temp.a.fn(years, "fall.calib.biomass", biomass.estimates$fall)
  biomass.res = cbind(biomass.res, temp.a.fn(years, "spring.calib.biomass", biomass.estimates$spring))
  colnames(biomass.res) = c("fall","spring")
  rownames(biomass.res) = years
  write.csv(round(biomass.res/1000,0), file = paste0(paste0(filename.mod,stock),".biomass.csv"))

  biomass.unc = temp.a.fn(years, "fall.biomass", biomass.estimates$fall)
  biomass.unc = cbind(biomass.unc, temp.a.fn(years, "spring.biomass", biomass.estimates$spring))
  temp = round(biomass.unc/biomass.res,2)
  colnames(temp) = c("fall", "spring")
  rownames(temp) = years
  write.csv(temp, file = paste0(paste0(filename.mod,stock),".biomass.efficiency.csv"))
  
  temp =cbind(all.fall.N.W[2,match(as.character(years), colnames(all.fall.N.W))],all.spring.N.W[2,match(as.character(years), colnames(all.spring.N.W))])
  temp = round(temp/biomass.unc,2)
  colnames(temp) = c("fall", "spring")
  rownames(temp) = years
  write.csv(temp, file = paste0(paste0(filename.mod,stock),".wal.biomass.ratios.csv"))

  #print(biomass.estimates)
  fall.Nal = temp.b.fn(years, var = "fall.Nal", biomass = biomass.estimates$fall)
  rownames(fall.Nal) = calibrate.length.integer
  colnames(fall.Nal) = years
  write.csv(fall.Nal, file = paste0(paste0(filename.mod,stock),".fall.Nal.csv"))
  fall.calib.Nal = temp.b.fn(years, var = "cal.fall.Nal", biomass = biomass.estimates$fall)
  rownames(fall.calib.Nal) = calibrate.length.integer
  colnames(fall.calib.Nal) = years
  write.csv(fall.calib.Nal, file = paste0(paste0(filename.mod,stock),".fall.calib.Nal.csv"))

  #print(biomass.estimates$spring)
  spring.Nal = temp.b.fn(years, var = "spring.Nal", biomass = biomass.estimates$spring)
  #stop()
  rownames(spring.Nal) = calibrate.length.integer
  colnames(spring.Nal) = years
  write.csv(spring.Nal, file = paste0(paste0(filename.mod,stock),".spring.Nal.csv"))
  spring.calib.Nal = temp.b.fn(years, var = "cal.spring.Nal", biomass = biomass.estimates$spring)
  rownames(spring.calib.Nal) = calibrate.length.integer
  colnames(spring.calib.Nal) = years
  write.csv(spring.calib.Nal, file = paste0(paste0(filename.mod,stock),".spring.calib.Nal.csv"))

  dyn.unload(dynlib("../wal"))
  return(biomass.estimates)
}
