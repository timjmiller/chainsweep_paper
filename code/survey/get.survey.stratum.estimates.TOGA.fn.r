get.survey.stratum.estimates.TOGA.fn <- function(spp = 105, survey = 200807, oc = sole, strata = paste0("01", 13:21, "0"), lengths, do.length = TRUE, 
	do.age = FALSE, gcf = 1, dcf = 1, vcf = 1, tow_swept_area = 0.007, len.data = NULL, 
  typecode = 1, operationcode = 3, gearcode = 2){

	#gcf.butterfish <- 1/1.35
	#gcf.butterfish.se <- 0.15
	#tow.size <- 0.01 #sq. nm
		
	stratum.sizes.q <- paste0("select stratum, stratum_area, strgrp_desc, stratum_name from svmstrata where STRATUM IN ('", 
		paste(strata, collapse = "','"), "')"," order by stratum")
	str.size <- sqlQuery(oc,stratum.sizes.q)
	str.size$NTOWS <- str.size$STRATUM_AREA/tow_swept_area
	#print(str.size)
	strata <- paste0('0', as.character(sort(unique(str.size$STRATUM))))

	#STATION
	q.sta <- paste0("select cruise6, stratum, tow, station, type_code, operation_code, gear_code, svvessel, svgear, est_year, est_month, est_day, ",
		"gmt_year, gmt_month, gmt_day, gmt_time, towdur, dopdistb, dopdistw, avgdepth, ",
		"area, bottemp, beglat, beglon from union_fscs_svsta ",
		"where cruise6 = ", survey, " and STRATUM IN ('", paste(strata, collapse = "','"), "')",
		" and TYPE_CODE <= ", typecode, 
	  " and OPERATION_CODE <= ", operationcode, 
	  " and GEAR_CODE <= ", gearcode, " order by cruise6, stratum, tow, station", sep = '')
	sta.view <- sqlQuery(oc,q.sta,stringsAsFactors= FALSE) 
	#print(dim(sta.view))
	temp <- str.size[match(sta.view$STRATUM, str.size$STRATUM),]
	sta.view <- cbind(sta.view, temp[,-1])
	
	#CATCH
	q.cat <- paste("select cruise6, stratum, tow, station, svspp, catchsex, expcatchwt, expcatchnum from union_fscs_svcat ",
		"where cruise6 = ", survey, " and STRATUM IN ('", paste(strata, collapse = "','"), "')",
		"and svspp = ", spp, " order by cruise6, stratum, tow, station", sep = '')
	cat.view <- sqlQuery(oc,q.cat,stringsAsFactors= FALSE)
	#cat.view$ID <- as.character(cat.view$ID)
	catch.data <- merge(sta.view, cat.view, by = c('CRUISE6','STRATUM','TOW','STATION'),  all.x = T, all.y=F)
	catch.data$EXPCATCHNUM[is.na(catch.data$EXPCATCHNUM)] <- 0
	catch.data$EXPCATCHWT[is.na(catch.data$EXPCATCHWT)] <- 0
	#gear conversion
	catch.data$EXPCATCHNUM[which(is.element(catch.data$SVGEAR, c(41,45)))] <- gcf * catch.data$EXPCATCHNUM[which(is.element(catch.data$SVGEAR, c(41,45)))]
	catch.data$EXPCATCHWT[which(is.element(catch.data$SVGEAR, c(41,45)))] <- gcf * catch.data$EXPCATCHWT[which(is.element(catch.data$SVGEAR, c(41,45)))]
	#door conversion
	catch.data$EXPCATCHNUM[which(catch.data$YEAR< 1985)] <- dcf * catch.data$EXPCATCHNUM[which(catch.data$YEAR< 1985)]
	catch.data$EXPCATCHWT[which(catch.data$YEAR< 1985)] <- dcf * catch.data$EXPCATCHWT[which(catch.data$YEAR< 1985)]
	#vessel conversion
	catch.data$EXPCATCHNUM[which(catch.data$SVVESSEL == 'DE')] <- vcf * catch.data$EXPCATCHNUM[which(catch.data$SVVESSEL == 'DE')]
	catch.data$EXPCATCHWT[which(catch.data$SVVESSEL == 'DE')] <- vcf * catch.data$EXPCATCHWT[which(catch.data$SVVESSEL == 'DE')]

	m <- sapply(str.size$STRATUM, function(x) sum(sta.view$STRATUM == x))
	M <- str.size$NTOWS
	samp.tot.n.w <- t(sapply(str.size$STRATUM, function(x) apply(catch.data[which(catch.data$STRATUM== x),c('EXPCATCHNUM','EXPCATCHWT')],2,sum)))
#	print(cbind(str.size$STRATUM,M,m,samp.tot.n.w))
	S.n.w.stratum <- t(sapply(str.size$STRATUM, function(x) var(catch.data[which(catch.data$STRATUM== x),c('EXPCATCHNUM','EXPCATCHWT')])))
#	print(cbind(M,m,S.n.w.stratum))
	N.W.hat.stratum <- M * samp.tot.n.w/m
	Vhat.N.W.hat.stratum <- M^2 * (1 - m/M) * S.n.w.stratum/m
	n.strata <- length(M)
	if(do.length){
		#LENGTH
		q.len <- paste("select  cruise6, stratum, tow, station, catchsex, length, expnumlen from union_fscs_svlen " ,
			"where cruise6 = ", survey, " and STRATUM IN('", paste(strata, collapse = "','"), "')",
			"and svspp = ", spp, " order by cruise6, stratum, tow, station, svspp, catchsex", sep = '')
		if(is.null(len.data))
		{
  		len.view <- sqlQuery(oc,q.len,stringsAsFactors= FALSE)
		  len.data <- merge(catch.data, len.view, by = c('CRUISE6','STRATUM','TOW','STATION','CATCHSEX'),  all.x=T, all.y = F)
		  len.data$EXPNUMLEN[is.na(len.data$EXPNUMLEN)] <- 0
		}
		if(max(len.data$LENGTH, na.rm= T) > max(lengths)) 
    {
      warning(paste0('max of lengths in data from survey ', survey, ' = ', max(len.data$LENGTH, na.rm= T), ' whereas max of lengths given is ', 
        max(lengths), ". Will change length range to accommodate.", sep = ''))
      lengths = sort(c(lengths, (max(lengths)+1):ceiling(max(len.data$LENGTH, na.rm= T))))
    }
		if(min(len.data$LENGTH, na.rm= T) < min(lengths)) 
    {
      warning(paste0('min of lengths in data from survey ', survey, ' = ', min(len.data$LENGTH, na.rm= T), ' whereas min of lengths given is ', 
      min(lengths), ". Will change length range to accommodate.", sep = ''))
      lengths = sort(c(floor(min(len.data$LENGTH, na.rm= T)):(min(lengths)-1), lengths))
    }
		
		samp.tot.nal <- sapply(lengths, function(x) sapply(str.size$STRATUM, function(y) sum(len.data$EXPNUMLEN[len.data$STRATUM == y & len.data$LENGTH == x],
			na.rm = TRUE)))
		Nal.hat.stratum <- M * samp.tot.nal/m
#		print(Nal.hat.stratum)
#		print(cbind(apply(Nal.hat.stratum,1,sum), N.W.hat.stratum))
		S.nal.stratum <- t(sapply(str.size$STRATUM, function(x){
			x.dat <- len.data[which(len.data$STRATUM == x),]
			if(dim(x.dat)[1]){
				nal.by.tow <- t(sapply(unique(x.dat$TOW), function(y){
					if(sum(x.dat$TOW == y)){
						nal.tow <- sapply(lengths, function(z) sum(x.dat$EXPNUMLEN[which(x.dat$LENGTH == z & x.dat$TOW == y)], na.rm = TRUE))
					}
					else nal.tow <- rep(0,length(lengths))
					return(nal.tow)
				}))
				cov.nal <- cov(nal.by.tow)
				return(cov.nal)
			}
			else return(matrix(NA,length(lengths),length(lengths)))
		}))
		Vhat.Nal.stratum <- M^2 * (1 - m/M) * S.nal.stratum/m
    if(any(m ==0)) warning(paste0("Some strata for survey have zero tows, will extrapolate data to these areas"))

				
		if(do.age){
			#AGE
			q.age <- paste("select  cruise6, stratum, tow, station, sex, length, age, indwt, maturity from union_fscs_svbio ",
				"where cruise6 = ", survey, " and STRATUM IN('", paste(strata, collapse = "','"), "')",
				"and svspp = ", spp, " and age is not null order by cruise6, stratum, tow, station", sep = '')
			age.view <- sqlQuery(oc,q.age,stringsAsFactors= FALSE)
      #print(q.age)
      #print(dim(age.view))
      #print(head(age.view))
			age.data <- merge(len.data, age.view, by = c('CRUISE6','STRATUM','TOW','STATION','LENGTH'),  all.x = T, all.y=F)
		}
	}

	
	
	out <- cbind.data.frame(stratum = str.size$STRATUM, M = M, m = m, mean.n.w.per.tow = N.W.hat.stratum/M, V.mean.n.w.per.tow = Vhat.N.W.hat.stratum/(M^2))
	out <- list(out = out)
	out$sta.view = sta.view
	out$stratum.size <- str.size
	out$catch.data = catch.data
  out$lengths = lengths
	if(do.length) {
		out$Nal.hat.stratum = Nal.hat.stratum
		out$V.Nal.stratum = Vhat.Nal.stratum
	  out$length.data = len.data
	}
	if(do.age) out$age.data <- age.data
  #odbcCloseAll()
	return(out)
}
