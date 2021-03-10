get.lendat.dn.TOGA.fn = function(survey,lendat=NULL,svspp=105, lens = 1:50, strata = svspp.strata, typec, operc, gearc)
{
  #source("get.survey.stratum.estimates.TOGA.fn.r")
  #print(is.null(lendat))
  #stop()
  if(is.null(lendat)) lendat <- get.survey.stratum.estimates.TOGA.fn(spp = svspp, survey = survey, strata = strata, 
    do.length = TRUE, lengths = lens, do.age = TRUE, typecode = typec, operationcode = operc, gearcode = gearc)$length.data
  if(is.null(lendat$id)) lendat$id = paste0(lendat$STRATUM,"_",lendat$CRUISE6, "_", lendat$STATION)
  lendat$GMT = paste0(lendat$GMT_MONTH,'/',lendat$GMT_DAY,'/',lendat$GMT_YEAR, " ", lendat$GMT_TIME)
  #print(head(lendat$GMT))
  lendat$GMT = as.POSIXct(lendat$GMT, format = "%m/%d/%Y %H:%M:%S", tz = "GMT")
  #print(head(lendat$GMT))
  lendat$latitude <- as.numeric(substr(lendat$BEGLAT,1,2))+as.numeric(substr(lendat$BEGLAT,3,nchar(lendat$BEGLAT)))/60
  lendat$longitude <- -as.numeric(substr(lendat$BEGLON,1,2))+as.numeric(substr(lendat$BEGLON,3,nchar(lendat$BEGLON)))/60
#  day = lendat$GMT_DAY
#  month = lendat$GMT_MONTH
#  year = lendate$GMT_YEAR
  year <- as.numeric(strftime(lendat$GMT, format="%Y", tz = "GMT"))
  month <- as.numeric(strftime(lendat$GMT, format="%m", tz = "GMT"))
  day <- as.numeric(strftime(lendat$GMT, format="%d", tz = "GMT"))
  hour <- as.numeric(strftime(lendat$GMT, format="%H", tz = "GMT"))
  minute <- as.numeric(strftime(lendat$GMT, format="%M", tz = "GMT"))
  second <- as.numeric(strftime(lendat$GMT, format="%S", tz = "GMT"))
  lat <- lendat$latitude
  long <- lendat$longitude
  #source("~/work/paired_tow_studies/R/2017/sunPosition.fn.r")
  res <- sunPosition.fn(day = day, month = month, year = year, hour = hour, minute = minute, second = second, timezone = rep(0, length(hour)), lat = lat, lon = long)
  day.night <- vector(length = length(res$PAR))
  day.night[res$PAR > 0] <- "day"
  day.night[res$PAR == 0] <- "night"
  lendat = cbind(lendat, res, day.night)
  return(lendat)
}

