boot.lendat.fn = function(lendat, dif=0)
{
  #print(head(lendat))
  #print(dim(lendat))
  if(is.null(lendat$id)) lendat$id = paste0(lendat$STRATUM,"_",lendat$CRUISE6, "_", lendat$STATION)
  strata = sort(unique(lendat$STRATUM))
  #print(strata)
  boot.ids = unlist(sapply(strata, function(x) {
    stratum.dat = lendat[which(lendat$STRATUM == x),]
    ids = unique(stratum.dat$id)
    #print(x)
    #print(ids)
    sample(ids, size = length(ids)-dif, replace = TRUE)
  }))
  #print(boot.ids)
  boot.index = unlist(sapply(boot.ids, function(x) which(lendat$id == x)))
  #print(boot.index)
  boot.nindex = sapply(boot.ids, function(x) sum(lendat$id == x))
  #print(boot.nindex)
  #print(length(boot.nindex))
  #print(length(boot.ids))
  boot.ids = rep(1:length(boot.ids), boot.nindex)
  #print(length(boot.ids))
  #print(length(boot.index))
  boot.lendat = lendat[boot.index,]
  boot.lendat$id = boot.ids
  return(boot.lendat)
}

