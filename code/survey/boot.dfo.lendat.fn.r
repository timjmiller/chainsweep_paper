boot.dfo.lendat.fn = function(cal.dat, dif=0)
{
  strata = sort(unique(cal.dat$STRATA))
  boot.index = unlist(sapply(strata, function(x) {
    ih = which(cal.dat$STRATA == x)
    sample(ih, size = length(ih)-dif, replace = TRUE)
  }))
  boot.cal.dat = cal.dat[boot.index,]
  boot.cal.dat$SET = 1:length(boot.index)
  return(boot.cal.dat)
}

