get.dfo.survey.Nal.hat.fn <- function(lengths, cal.dat, str.size, tow_swept_area = 0.011801)
{
	M <- str.size$trawl_units
	strata <- str.size$STRATA[-length(str.size$STRATA)]
	m <- sapply(strata, function(x) sum(cal.dat$STRATA == x))
	strata = strata[which(m>0)]
	m = m[which(m>0)]
	M = M[which(m>0)]
	length.cols = which(!(colnames(cal.dat) %in% c("STRATA","SLAT","SLONG","UNITAREA","SET","TOTAL","day.night")))
	lengths.present = as.numeric(substr(colnames(cal.dat),2,nchar(colnames(cal.dat)))[length.cols])
	n.strata <- length(M)

	if(max(lengths.present) > max(lengths)) warning(paste('max of lengths in length data = ', max(lengths.present), ' whereas max of lengths given is ', max(lengths), sep = ''))
	
	samp.tot.nal <- sapply(lengths, function(x) sapply(strata, function(y) 
  {
    if(x %in% lengths.present) sum(cal.dat[cal.dat$STRATA == y, length.cols[which(lengths.present == x)]])
    else 0
	}))
	Nal.hat.stratum <- M * samp.tot.nal/m
	S.nal.stratum <- t(sapply(strata, function(x)
  {
    stratum.dat <- cal.dat[which(cal.dat$STRATA == x),]
    if(dim(stratum.dat)[1])
    {
      nal.by.tow <- t(sapply(1:NROW(stratum.dat), function(y)
      {
        nal.tow <- sapply(lengths, function(z)
        { 
          if(z %in% lengths.present) 
          {
            stratum.dat[y,length.cols[which(lengths.present == z)]]
          }
          else 0
        })
      }))
      cov.nal <- cov(nal.by.tow)
      return(cov.nal)
    }
    else return(matrix(NA,length(lengths),length(lengths)))
  }))
	Vhat.Nal.stratum <- M^2 * (1 - m/M) * S.nal.stratum/m

	out <- cbind.data.frame(stratum = strata, M = M, m = m)
	out <- list(out = out)
	out$Nal.hat.stratum = Nal.hat.stratum
	out$V.Nal.stratum = Vhat.Nal.stratum
  out$length.data = cal.dat
  out$all.strata = str.size
  out$Nal.hat = apply(Nal.hat.stratum,2,sum)
  out$V.Nal.hat = matrix(apply(Vhat.Nal.stratum,2,sum), length(lengths),length(lengths))
	return(out)
}

