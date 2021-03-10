get.best = function(sp.name,pdir)
{
  out = list()
  x = readRDS(paste0(pdir, "/results/big_results/", sp.name, "_model_fits.RDS"))
  out$convergence.flags = sapply(x$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, y$model.res$opt$conv))
  #aic for the models
  out$np = sapply(x$model.fits, function(y) length(y$model.res$par))
  out$aic = sapply(x$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 
    2*(ifelse(abs(y$model.res$opt$obj)<1e8,y$model.res$opt$obj,NA) + length(y$model.res$opt$par))))
  #aic = sapply(x$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 2*(y$model.res$opt$obj + length(y$model.res$opt$par))))
  out$aic.converged = out$aic[which(out$convergence.flags==0)]
  out$best.model = names(out$aic.converged)[which(out$aic.converged == min(out$aic.converged, na.rm = TRUE))] 
  if("day" %in% names(x))
  {
    #aic.day = sapply(x$day$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 2*(y$model.res$opt$obj + length(y$model.res$opt$par))))
    out$aic.day = sapply(x$day$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 
      2*(ifelse(abs(y$model.res$opt$obj)<1e8,y$model.res$opt$obj,NA) + length(y$model.res$opt$par))))
    out$convergence.flags.day = sapply(x$day$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, y$model.res$opt$conv))
    out$aic.day.converged = out$aic.day[which(out$convergence.flags.day == 0)]
  }
  if("night" %in% names(x))
  {
    #aic.night = sapply(x$night$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 2*(y$model.res$opt$obj + length(y$model.res$opt$par))))
    out$aic.night = sapply(x$night$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 
      2*(ifelse(abs(y$model.res$opt$obj)<1e8,y$model.res$opt$obj,NA) + length(y$model.res$opt$par))))
    out$convergence.flags.night = sapply(x$night$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, y$model.res$opt$conv))
    out$aic.night.converged = out$aic.night[which(out$convergence.flags.night == 0)]
  }
  if("night" %in% names(x) & "day" %in% names(x))
  {
    #same model applied to day and night
    out$aic.day.night = sapply(x$day$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 
      2*(ifelse(abs(y$model.res$opt$obj)<1e8,y$model.res$opt$obj,NA) + length(y$model.res$opt$par)))) +
      sapply(x$night$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 
      2*(ifelse(abs(y$model.res$opt$obj)<1e8,y$model.res$opt$obj,NA) + length(y$model.res$opt$par))))
    #aic.day.night = sapply(x$day$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 2*(y$model.res$opt$obj + length(y$model.res$opt$par))))+
    #sapply(x$night$model.fits, function(y) ifelse(is.character(y$model.res$opt), NA, 2*(y$model.res$opt$obj + length(y$model.res$opt$par))))
    out$aic.day.night.converged = out$aic.day.night[which(out$convergence.flags.night == 0 & out$convergence.flags.day == 0)]
    out$best.dn.model = names(out$aic.day.night.converged)[which(out$aic.day.night.converged == min(out$aic.day.night.converged, na.rm = TRUE))] 
  }
  setwd(pdir)
  return(out)
  #return(list(best = c(best.model, best.dn.model), converge = convergence.flags, converge.d = convergence.flags.day, converge.n = convergence.flags.night,
  #  aic = aic, aic.d = aic.day, aic.n = aic.night, aic.dn = aic.day.night, np = np))
}
