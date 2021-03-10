
estimate_biomass = function(sp.i = 1, stock.i = 1, sp.info, twintrawl.res, do.N.W = TRUE, years = NULL, pred.length.integer = NULL,
min.length.calibrate, max.length.calibrate, filename.mod = "", TOGA.type = 1, TOGA.operation = 3, TOGA.gear = 2, pdir)
{
  
  #added in max.length.calibrate argument
  #now length obs and len-weight obs data are generated elsewhere and read in here.

  ####################################################
  ############### get survey data ####################
  ####################################################
  svspp = sp.info$SVSPP[sp.i]
  sp = sp.info$sp.names[sp.i]
  stock = stocks[stock.i]
  x = get.best(sp)
  best.model = x$best.model #"bi3"
  setwd(sp)
 
  #read in survey data generated by gather_survey_data.r
  load(paste0(pdir, "/data/survey/", stock, "_lendat.RData"))
  load(paste0(pdir, "/data/survey/", stock, "_lw_dat.RData"))  
  load(paste0(pdir, "/data/survey/", stock, "_N.W.RData"))
  
  if(is.null(pred.length.integer)) pred.length.integer = floor(min(twintrawl.dn.res$pred.length)):ceiling(max(twintrawl.dn.res$pred.length))

  min.lengths = min(sapply(fall.lengths, min),sapply(spring.lengths, min))
  max.lengths = max(sapply(fall.lengths, max),sapply(spring.lengths, max))
  if(min(pred.length.integer)>min.lengths || max(pred.length.integer) < max.lengths)
  {
    warning("adjusting pred.length.integer to deal with ranges of survey data lengths")
    if(min(pred.length.integer)>min.lengths) pred.length.integer = sort(c(min.lengths:(min(pred.length.integer)-1), pred.length.integer))
    if(max(pred.length.integer) < max.lengths) pred.length.integer = sort(c(pred.length.integer, (max(pred.length.integer)+1):max.lengths))
  }

	strata = fall.strata[[stock]])$stratum.size
	strata = spring.strata[[stock]])$stratum.size

  
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