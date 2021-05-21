plot_biomass = function(stocks, stock.names.plot, all.boots)#i,sp.info, stock, file.loc = "paper")# = rep(6,6))
{
  temp <- col2rgb('dodgerblue')
  scol <- paste0(rgb(temp[1,],temp[2,], temp[3,], maxColorValue = 255), "50")
  scoli <- paste0(rgb(temp[1,],temp[2,], temp[3,], maxColorValue = 255), "20")
  temp <- col2rgb('orangered')
  fcol <- paste0(rgb(temp[1,],temp[2,], temp[3,], maxColorValue = 255), "50")
  fcoli <- paste0(rgb(temp[1,],temp[2,], temp[3,], maxColorValue = 255), "20")
  
  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  ymax = numeric()
  
  plot.data = lapply(stocks, function(x){
    y = read.csv(paste0("results/", x, ".biomass.csv"))
    #print(y)
    z = read.csv(paste0("results/", x, ".cv.biomass.csv"))
    #print(z)
    years = y[,1]
    #print(years)
    fall.biomass = y[,2]/1000
    fall.biomass.lo = z[,3]/1000
    fall.biomass.hi = z[,4]/1000
    spring.biomass = y[,3]/1000
    spring.biomass.lo = z[,6]/1000
    spring.biomass.hi = z[,7]/1000
    return(cbind(years, fall.biomass,fall.biomass.lo, fall.biomass.hi,spring.biomass,spring.biomass.lo, spring.biomass.hi))
  })
  ymax = c(rep(max(sapply(plot.data[c(1:12)], function(x) max(x[,-1], na.rm = TRUE))),12),
    rep(max(sapply(plot.data[c(13:15)], function(x) max(x[,-1], na.rm = TRUE))),3),
    rep(max(sapply(plot.data[c(16:17)], function(x) max(x[,-1], na.rm = TRUE))),2))
  years = plot.data[[1]][,1]
  #ymax = 200
  par(mfrow = c(6,3), oma = c(5,5,0,0), mar = c(0,0,3,0))
  for(i in 1:length(stocks)){
    spring.ind = which(!is.na(plot.data[[i]][,6]))
    fall.ind = which(!is.na(plot.data[[i]][,3]))
  #print(paste0(file.loc, "/", sp.name, "/", stock, "_biomass.png"))
  #png(filename = paste0(file.loc, "/", sp.name, "/", stock, "_biomass.png"), width = 12*144, height = 8*144, res = 144, pointsize = 12, family = "Times")#,
    plot(years, plot.data[[i]][,5], type = 'n', ylim = c(0, ymax[i]), axes = FALSE, ann = FALSE)
    box(lwd = 2)
    grid(col = gray(0.7), lwd = 1, lty = 2)
    lines(years[spring.ind], plot.data[[i]][spring.ind,5], lwd = 2, col = "dodgerblue")
    polygon(c(years[spring.ind],rev(years[spring.ind])), c(plot.data[[i]][spring.ind,6],rev(plot.data[[i]][spring.ind,7])), col = scol, border = "transparent")
    if(length(fall.ind) != NROW(plot.data[[i]]))
    {
      tyrs = 2009:2016
      ind = which(survey.years %in% tyrs)
      tyrs = c(tyrs, 2016.2)
      tdat = plot.data[[i]][ind,]
      tdat = rbind(plot.data[[i]][ind,], plot.data[[i]][ind[length(ind)],])
      lines(2009:2016, plot.data[[i]][ind,2], lwd = 2, col = "orangered")
      polygon(c(tyrs,rev(tyrs)), c(tdat[,3],rev(tdat[,4])), col = fcol, border = "transparent")

      tyrs = 2018:2019
      ind = which(survey.years %in% tyrs)
      tyrs = c(2017.8,tyrs)
      tdat = plot.data[[i]][ind,]
      tdat = rbind(plot.data[[i]][ind[1],],plot.data[[i]][ind,])
      lines(2018:2019, plot.data[[i]][ind,2], lwd = 2, col = "orangered")
      polygon(c(tyrs,rev(tyrs)), c(tdat[,3],rev(tdat[,4])), col = fcol, border = "transparent")
    }
    else{
      lines(years[fall.ind], plot.data[[i]][fall.ind,2], lwd = 2, col = "orangered")
      polygon(c(years[fall.ind],rev(years[fall.ind])), c(plot.data[[i]][fall.ind,3],rev(plot.data[[i]][fall.ind,4])), col = fcol, border = "transparent")
    }
    if(i %in% 16:18) axis(1, lwd = 2, cex.axis = 1.5)
    else axis(1, lwd = 2, cex.axis = 1.5, tick = FALSE, labels = FALSE)
    if(i %in% c(1,4,7,10,13,16)) axis(2, lwd = 2, cex.axis = 1.5)
    else axis(2, lwd = 2, cex.axis = 1.5, tick = FALSE, labels = FALSE)
    mtext(side = 3, stock.names.plot[i], line = 0, cex = 1)
#  if(i ==6) stop()
    #axis(2, lwd = 2, cex.axis = 1.5)
    #axis(1, lwd = 2, cex.axis = 1.5)
    #axis(2, lwd = 2, cex.axis = 1.5, labels = FALSE)
  }

#  plot(years, fall.biomass, type = 'n', ylim = c(0, ymax[i]), axes = FALSE, ann = FALSE)
#  box(lwd = 2)
#  grid(col = gray(0.7), lwd = 1, lty = 2)
  #mtext(side = 3, "Fall", line = 1, cex = 1.5)
  #stock.name.text = use.stock.names[i]
  #mtext(side = 4, stock.name.text, line = 3, cex = 1.5)
  mtext(side = 2, parse(text = paste0("Biomass",  "~(10^3", "~ mt)")), line = 2, cex = 1.5, outer = TRUE)
  mtext(side = 1, 'Year', line = 3, cex = 1.5, outer = TRUE)
  #dev.off()
  
}
