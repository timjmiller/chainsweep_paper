get.survey.Nal.hat.fn <- function(len.data, str.size, lengths = 1:50)
{

  M = str.size$NTOWS
  m = sapply(str.size$STRATUM, function(x) length(unique(len.data$id[len.data$STRATUM == x])))
  #print(M)
  #print(m)
	out = list(str.size = str.size)
	#print(head(len.data))
  #print(M)
  #print(lengths)
  #print(m)
  samp.tot.nal <- sapply(lengths, function(x) sapply(str.size$STRATUM, function(y) sum(len.data$EXPNUMLEN[len.data$STRATUM == y & len.data$LENGTH == x],
	  na.rm = TRUE)))
  #print(dim(samp.tot.nal))
  Nal.hat.stratum <- M * samp.tot.nal/m
  
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

	out$expand = sum(M)/sum(M[which(m>0)]) #=1 if all strata are sampled
	out$Nal.hat.stratum = Nal.hat.stratum
	out$V.Nal.stratum = Vhat.Nal.stratum
  out$length.data = len.data
  out$Nal.hat = apply(Nal.hat.stratum,2,sum, na.rm = TRUE) * out$expand
  out$V.Nal.hat = matrix(apply(Vhat.Nal.stratum,2,sum, na.rm = TRUE), length(lengths),length(lengths))  * out$expand^2
	return(out)
}
