make.all.AIC.table.fn = function(sp.names, pretty.sp.names, latex.table.file = '', pdir)
{

  library(Hmisc)
  models = paste0(rep(c("BI","BB"), c(5,8)), "$_{", c(0:4, 0:7), "}$")
  first.names = paste0(rep(c("bi", "bb"), c(5,8)), c(0:4,0:7))
  x.all = get.best(sp.names[1], pdir = pdir)
  print(x.all)
  x = lapply(x.all, function(y) y[which(names(y) %in% first.names)])
  print(x)
  #print(x)
  out = cbind("$n_p$" = x$np)
  print(out)
  print(models)
  rownames(out) = models
  aic = x$aic
  aic[which(!(names(aic) %in% names(x$aic.converged)))] = NA
  aic = aic - min(aic, na.rm = TRUE)
  out = cbind(out, aic)
  #colnames(out)[2] = "Summer flounder"#sp.names[1] 
  for(i in 2:length(sp.names))
  {
    print(i)
    sp.name = sp.names[i]
    x.all = get.best(sp.names[i], pdir = pdir)
    x = lapply(x.all, function(y) y[which(names(y) %in% first.names)])
    print("here")
    aic = x$aic
    aic[which(!(names(aic) %in% names(x$aic.converged)))] = NA
    aic = aic - min(aic, na.rm = TRUE)
    out = cbind(out, aic)
  }
  out = apply(out,2, round, 2)
  #colnames(out)[3:7] = c("American plaice", "Windowpane", "Winter flounder", "Yellowtail flounder", "Witch flounder")
  colnames(out)[1 + 1:length(sp.names)] = pretty.sp.names
  print(out)
  #x = latex(out, file = latex.table.file, 
  #  rgroup = c("Combined","Day/Night-specific"),
  #  n.rgroup = c(13,13),
  #  cgroup = c("", "$\\Delta$AIC"),
  #  n.cgroup = c(1,length(sp.names)),
  #  rowlabel = 'Model', table.env = FALSE, rowlabel.just = "c")
  return(out)
}
