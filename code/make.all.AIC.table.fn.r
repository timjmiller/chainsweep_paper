make.all.AIC.table.fn = function(sp.names, pretty.sp.names, latex.table.file = '~/work/paired_tow_studies/tex/2017/model_compare.tex')
{

  library(Hmisc)
  models = paste0(rep(c("BI","BB"), c(5,8)), "$_{", c(0:4, 0:7), "}$")
  x = get.best(sp.names[1])
  #print(x)
  out = cbind("$n_p$" = c(x$np, 2*(x$np)))
  #print(dim(out))
  #print(length(models))
  rownames(out) = rep(models, 2)
  aic = x$aic
  aic[which(!(names(aic) %in% names(x$aic.converged)))] = NA
  aic.dn = x$aic.day.night
  aic.dn[which(!(names(aic.dn) %in% names(x$aic.day.night.converged)))] = NA
  aic = c(aic, aic.dn)
  #aic[which(c(x$converge !=0, x$converge.d !=0 | x$converge.n != 0))] = NA
  aic = aic - min(aic, na.rm = TRUE)
  out = cbind(out, aic)
  #colnames(out)[2] = "Summer flounder"#sp.names[1] 
  for(i in 2:length(sp.names))
  {
    sp.name = sp.names[i]
    x = get.best(sp.names[i])
    aic = x$aic
    aic[which(!(names(aic) %in% names(x$aic.converged)))] = NA
    aic.dn = x$aic.day.night
    aic.dn[which(!(names(aic.dn) %in% names(x$aic.day.night.converged)))] = NA
    aic = c(aic, aic.dn)
    #aic[which(c(x$converge !=0, x$converge.d !=0 | x$converge.n != 0))] = NA
    aic = aic - min(aic, na.rm = TRUE)
    out = cbind(out, aic)
  }
  out = apply(out,2, round, 2)
  #colnames(out)[3:7] = c("American plaice", "Windowpane", "Winter flounder", "Yellowtail flounder", "Witch flounder")
  colnames(out)[1 + 1:length(sp.names)] = pretty.sp.names
  print(out)
  x = latex(out, file = latex.table.file, 
    rgroup = c("Combined","Day/Night-specific"),
    n.rgroup = c(13,13),
    cgroup = c("", "$\\Delta$AIC"),
    n.cgroup = c(1,length(sp.names)),
    rowlabel = 'Model', table.env = FALSE, rowlabel.just = "c")
  return(out)
}
