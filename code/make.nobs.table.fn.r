make.nobs.table.fn = function(sp.dat=sp.info)
{

  library(Hmisc)
  
  out = matrix(NA, NROW(sp.dat), 6)
  rownames(out) = sp.dat$sp.pretty.names
  colnames(out) = paste0(rep(c("Day", "Night", "Total"),2), rep(c(" Stations", " Fish"), each = 3))
  for(i in 1:NROW(out)) 
  {
    tt.dat = readRDS(paste0("data/", sp.info$sp.names[i], "_data.RDS"))
    #tt.dat = readRDS(paste0("/results/big_results/", sp.info$sp.names[i], "_fits.RDS"))
    out[i,1] = length(twintrawl.dn.res$day$model.fits$bi0$model.res$env$data$n_per_sta)
    out[i,2] = length(twintrawl.dn.res$night$model.fits$bi0$model.res$env$data$n_per_sta)
    out[i,4] = sum(twintrawl.dn.res$day$model.fits$bi0$model.res$env$data$n)
    out[i,5] = sum(twintrawl.dn.res$night$model.fits$bi0$model.res$env$data$n)
  }
  out[,3] = apply(out[,1:2], 1, sum)
  out[,6] = apply(out[,4:5], 1, sum)
  
  #x = latex(out, file = '', 
  #  rowlabel = 'Species', table.env = FALSE, rowlabel.just = "c")
  #write.csv(out, file= "nobs_table.csv")
  return(out)
}
