plot.results = function(sp.name, i, ymax = 6, xlim, sp.tab= sp.info, ylab = "Relative Catch Efficiency\n(Chain:Rockhopper)")
{
  parentdir = getwd()
  plot.input.fn = function(x)
  {
    y = summary(x$sdrep)
    ind = which(rownames(y) == "mean_pred_eta")
    z = y[ind,]
    z = cbind(z[,1],z[,1]-qnorm(0.975)*z[,2],z[,1]+qnorm(0.975)*z[,2])
    return(z)
  }
  #require(plotrix)
  #sp.name = sp.tab$sp.names[i]
  #sp = sp.tab$SVSPP[i]
  x = get.best(sp.name)
  best.model = x$best.model 
  print(best.model)
  #best.dn.model = x$best.dn.model #"bi3"
  print(best.model)
  setwd(sp.name)
  x = readRDS("model_fits.RDS")
  load("combined_data.RData")  
  boot = "boot.pred.eta.RDS" %in% dir()
  if(boot) boot.pred.eta = readRDS("boot.pred.eta.RDS")
  dat = combined.data$cl.less
  dat$dn = dat$day.night
  temp = aggregate(dn ~id, dat, FUN = unique)[,2]
  day.sta.ind = which(temp == "day")
  night.sta.ind = which(temp == "night")
  #no.converge = apply(boot.pred.eta, 2, function(x) sum(is.na(x)))
  temp = summary(x$model.fits[[best.model]]$model.res$sdrep)
  print(dim(temp))
  print(unique(rownames(temp)))
  plen = x$pred.length
  print(length(plen))
  #see what CV's at length are like for most complex model. CV is essentially the same as the standard error on log scale
  temp = temp[which(rownames(temp) == "mean_pred_eta"),]
  if(NROW(temp) > length(plen)) dn = TRUE else dn = FALSE #NROW(temp) should be a multiple of length(plen)
  
  y = x$model.fits[[best.model]]$model.res$rep
  z = plot.input.fn(x$model.fits[[best.model]]$model.res)
  ind = cbind(1:length(plen), length(plen) + 1:length(plen))
  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  if(missing(xlim)) xlim = range(plen)
  png(filename = paste0(sp.name,"_", best.model, "_rho_1.png"), width = 8*144, height = 8*144, res = 144, pointsize = 12, family = "Times")#,
    par(oma = c(1,1,1,1), mar = c(4,4,0,0), mfrow = c(1,1))
    plot(plen, exp(z[ind[,1],1]), type = 'n', ylim = c(0, ymax[1]), xlim = xlim, axes = FALSE, ann = FALSE)
    box(lwd = 2)
    grid(col = gray(0.7), lwd = 1, lty = 2)
    print(day.sta.ind)
    print(night.sta.ind)
    print(dim(y$station_pred_eta))
    for(i in day.sta.ind) lines(plen, exp(y$station_pred_eta[i,ind[,1]]), col = gray(0.8), lwd = 1)
    lines(plen, exp(z[ind[,1],1]), lwd = 2)
    polygon(c(plen,rev(plen)), exp(c(z[ind[,1],2],rev(z[ind[,1],3]))), col = tcol, border = "transparent")
    abline(h = 1, col = "red", lwd =2)
    #rug(jitter(dat$length[dat$dn == "day"], factor = 2), ticksize = 0.03)
    tdat = cbind(dat$length, (exp(-dat$offst) * dat$recnumlen.ch/dat$recnumlen.rh))[dat$dn == "day",]
    tdat = tdat[which(tdat[,2] > 0),]
    points(jitter(tdat[,1], factor = 2), tdat[,2], pch = 19)
    axis(1, lwd = 2, cex.axis = 1.5)
    axis(2, lwd = 2, cex.axis = 1.5)
    if(boot) lines(plen, exp(apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.025, na.rm = TRUE)), col = 'red', lty = 2)
    if(boot) lines(plen, exp(apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.975, na.rm = TRUE)), col = 'red', lty = 2)
    mtext(side = 2, "Relative Catch Efficiency\n(Chain:Rockhopper)", line = 2, cex = 1.5)
    mtext(side = 1, 'Length (cm)', line = 2.5, cex = 1.5)
  dev.off()
  if(dn)
  {
    png(filename = paste0(sp.name,"_", best.model, "_rho_2.png"), width = 8*144, height = 8*144, res = 144, pointsize = 12, family = "Times")#,
      par(oma = c(1,1,1,1), mar = c(4,4,0,0), mfrow = c(1,1))
      plot(plen, exp(z[ind[,2],1]), type = 'n', ylim = c(0, ymax[length(ymax)]), xlim = xlim, axes = FALSE, ann = FALSE)
      box(lwd = 2)
      grid(col = gray(0.7), lwd = 1, lty = 2)
      for(i in night.sta.ind) lines(plen, exp(y$station_pred_eta[i,ind[,2]]), col = gray(0.8), lwd = 1)
      lines(plen, exp(z[ind[,2],1]), lwd = 2)
      polygon(c(plen,rev(plen)), exp(c(z[ind[,2],2],rev(z[ind[,2],3]))), col = tcol, border = "transparent")
      abline(h = 1, col = "red", lwd =2)
      #rug(jitter(dat$length[dat$dn == "night"], factor = 2), ticksize = 0.03)
      tdat = cbind(dat$length, (exp(-dat$offst) * dat$recnumlen.ch/dat$recnumlen.rh))[dat$dn == "night",]
      tdat = tdat[which(tdat[,2] > 0),]
      points(jitter(tdat[,1], factor = 2), tdat[,2], pch = 19)
      axis(1, lwd = 2, cex.axis = 1.5)
      axis(2, lwd = 2, cex.axis = 1.5)
      if(boot) lines(plen, exp(apply(boot.pred.eta[,ind[,2]], 2, quantile, probs = 0.025, na.rm = TRUE)), col = 'red', lty = 2)
      if(boot) lines(plen, exp(apply(boot.pred.eta[,ind[,2]], 2, quantile, probs = 0.975, na.rm = TRUE)), col = 'red', lty = 2)
      mtext(side = 2, ylab, line = 2, cex = 1.5)
      mtext(side = 1, 'Length (cm)', line = 2.5, cex = 1.5)
    dev.off()
    png(filename = paste0(sp.name,"_", best.model, "_rho_1_2.png"), width = 16*144, height = 8*144, res = 144, pointsize = 12, family = "Times")
      par(oma = c(1,1,1,1), mar = c(4,4,3,0), mfrow = c(1,2))
      plot(plen, exp(z[ind[,1],1]), type = 'n', ylim = c(0, ymax[1]), xlim = xlim, axes = FALSE, ann = FALSE)
      box(lwd = 2)
      grid(col = gray(0.7), lwd = 1, lty = 2)
      for(i in day.sta.ind) lines(plen, exp(y$station_pred_eta[i,ind[,1]]), col = gray(0.8), lwd = 1)
      lines(plen, exp(z[ind[,1],1]), lwd = 2)
      polygon(c(plen,rev(plen)), exp(c(z[ind[,1],2],rev(z[ind[,1],3]))), col = tcol, border = "transparent")
      abline(h = 1, col = "red", lwd =2)
      #rug(jitter(x$night$cl.less$length, factor = 2), ticksize = 0.03)
      tdat = cbind(dat$length, (exp(-dat$offst) * dat$recnumlen.ch/dat$recnumlen.rh))[dat$dn == "day",]
      tdat = tdat[which(tdat[,2] > 0),]
      points(jitter(tdat[,1], factor = 2), tdat[,2], pch = 19)
      axis(1, lwd = 2, cex.axis = 1.5)
      axis(2, lwd = 2, cex.axis = 1.5)
      mtext(side = 2, ylab, line = 2, cex = 1.5)
      mtext(side = 3, "Day", line = 1, cex = 1.5)
      if(boot) lines(plen, exp(apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.025, na.rm = TRUE)), col = 'red', lty = 2)
      if(boot) lines(plen, exp(apply(boot.pred.eta[,ind[,1]], 2, quantile, probs = 0.975, na.rm = TRUE)), col = 'red', lty = 2)
      
      plot(plen, exp(z[ind[,2],1]), type = 'n', ylim = c(0, ymax[length(ymax)]), xlim = xlim, axes = FALSE, ann = FALSE)
      box(lwd = 2)
      grid(col = gray(0.7), lwd = 1, lty = 2)
      for(i in night.sta.ind) lines(plen, exp(y$station_pred_eta[i,ind[,2]]), col = gray(0.8), lwd = 1)
      lines(plen, exp(z[ind[,2],1]), lwd = 2)
      polygon(c(plen,rev(plen)), exp(c(z[ind[,2],2],rev(z[ind[,2],3]))), col = tcol, border = "transparent")
      abline(h = 1, col = "red", lwd =2)
      #rug(jitter(x$day$cl.less$length, factor = 2), ticksize = 0.03)
      tdat = cbind(dat$length, (exp(-dat$offst) * dat$recnumlen.ch/dat$recnumlen.rh))[dat$dn == "night",]
      tdat = tdat[which(tdat[,2] > 0),]
      points(jitter(tdat[,1], factor = 2), tdat[,2], pch = 19)
      axis(1, lwd = 2, cex.axis = 1.5)
      axis(2, lwd = 2, cex.axis = 1.5)
      mtext(side = 3, "Night", line = 1, cex = 1.5)
      mtext(side = 1, 'Length (cm)', line = -1, cex = 1.5, outer = TRUE)
      if(boot) lines(plen, exp(apply(boot.pred.eta[,ind[,2]], 2, quantile, probs = 0.025, na.rm = TRUE)), col = 'red', lty = 2)
      if(boot) lines(plen, exp(apply(boot.pred.eta[,ind[,2]], 2, quantile, probs = 0.975, na.rm = TRUE)), col = 'red', lty = 2)
    dev.off()
  }  
  setwd(parentdir)
  #return(cbind(no.converge, no.converge.day, no.converge.night))

}
