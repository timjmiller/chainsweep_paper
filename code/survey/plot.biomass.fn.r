plot.biomass.fn = function(i,sp.info, stock, plot.stock.name, file.loc = "~/work/paired_tow_studies/tex/2017")# = rep(6,6))
{
  sp.pretty.name = sp.info$sp.pretty.names[i]
  #use.sp.names = c("Summer flounder", "American plaice", "Windowpane", "Winter flounder", "Yellowtail flounder", "Witch flounder")
  #sp.names.ind = rep(1:6, c(1,1,2,3,3,1))
  sp.name = sp.info$sp.names[i]
  #sp.names = c("fluke", "plaice", "windowpane", "winter_flounder", "yellowtail_flounder", "witch_flounder")
  spp = sp.info$SVSPP[i]
  #spps = c(103, 102, 108, 106, 105, 107)
  #stock.names = c("fluke", "plaice", "gbgom_windowpane", "snemab_windowpane", "gb_winter_flounder", "gom_winter_flounder", "sne_winter_flounder", 
  #  "gb_yellowtail_flounder", "snema_yellowtail_flounder", "ccgom_yellowtail_flounder", "witch_flounder")
  use.stock.name = plot.stock.name
  #use.stock.names = c("Summer \n flounder", "American \n plaice", "GB-GOM \n windowpane", "SNE-MAB \n windowpane", "GB winter \n flounder", "GOM winter \n flounder", "SNE winter \n flounder", 
  #  "GB yellowtail \n flounder", "SNE-MA yellowtail \n flounder", "CC-GOM yellowtail \n flounder", "Witch \n flounder")

  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  ymax = numeric()
  y = read.csv(paste0(sp.name,"/2019_", stock, ".biomass.1.csv"))
  print(y)
  z = read.csv(paste0(sp.name,"/2019_", stock, ".cv.biomass.1.csv"))
  print(z)
  years = y[,1]
  #print(years)
  fall.biomass = y[,2]/1000
  fall.biomass.lo = z[,3]/1000
  fall.biomass.hi = z[,4]/1000
  spring.biomass = y[,3]/1000
  spring.biomass.lo = z[,6]/1000
  spring.biomass.hi = z[,7]/1000
  spring.ind = which(!is.na(spring.biomass.lo))
  fall.ind = which(!is.na(fall.biomass.lo))
  ymax[i] = max(spring.biomass.hi,fall.biomass.hi, na.rm = TRUE)
  print(paste0(file.loc, "/", sp.name, "/", stock, "_biomass_1.png"))
  png(filename = paste0(file.loc, "/", sp.name, "/", stock, "_biomass_1.png"), width = 12*144, height = 8*144, res = 144, pointsize = 12, family = "Times")#,
  par(mfrow = c(1,2), oma = c(4,5,3,4), mar = c(1,1,1,1))
  plot(years, spring.biomass, type = 'n', ylim = c(0, ymax[i]), axes = FALSE, ann = FALSE)
  box(lwd = 2)
  grid(col = gray(0.7), lwd = 1, lty = 2)
  lines(years[spring.ind], spring.biomass[spring.ind], lwd = 2)
  polygon(c(years[spring.ind],rev(years[spring.ind])), c(spring.biomass.lo[spring.ind],rev(spring.biomass.hi[spring.ind])), col = tcol, border = "transparent")
  axis(1, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5)
  mtext(side = 3, "Spring", line = 1, cex = 1.5)

  plot(years, fall.biomass, type = 'n', ylim = c(0, ymax[i]), axes = FALSE, ann = FALSE)
  box(lwd = 2)
  grid(col = gray(0.7), lwd = 1, lty = 2)
  lines(years[fall.ind], fall.biomass[fall.ind], lwd = 2)
  polygon(c(years[fall.ind],rev(years[fall.ind])), c(fall.biomass.lo[fall.ind],rev(fall.biomass.hi[fall.ind])), col = tcol, border = "transparent")
  axis(1, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5, labels = FALSE)
  mtext(side = 3, "Fall", line = 1, cex = 1.5)
  #stock.name.text = use.stock.names[i]
  #mtext(side = 4, stock.name.text, line = 3, cex = 1.5)
  mtext(side = 2, parse(text = paste0("Biomass",  "~(10^3", "~ mt)")), line = 2, cex = 1.5, outer = TRUE)
  mtext(side = 1, 'Year', line = 2, cex = 1.5, outer = TRUE)
  dev.off()
  
  png(filename = paste0(file.loc, "/", sp.name, "/", stock, "_biomass_2.png"), width = 12*144, height = 8*144, res = 144, pointsize = 12, family = "Times")#,
  par(mfrow = c(1,2), oma = c(4,5,3,4), mar = c(1,1,1,1))
  ymax = numeric()
  y = read.csv(paste0(sp.name,"/2019_", stock, ".biomass.2.csv"))
  z = read.csv(paste0(sp.name,"/2019_", stock, ".cv.biomass.2.csv"))
  years = y[,1]
  #print(years)
  fall.biomass = y[,2]/1000
  fall.biomass.lo = z[,3]/1000
  fall.biomass.hi = z[,4]/1000
  spring.biomass = y[,3]/1000
  spring.biomass.lo = z[,6]/1000
  spring.biomass.hi = z[,7]/1000
  spring.ind = which(!is.na(spring.biomass.lo))
  fall.ind = which(!is.na(fall.biomass.lo))
  ymax[i] = max(spring.biomass.hi,fall.biomass.hi, na.rm = TRUE)
    
  plot(years, spring.biomass, type = 'n', ylim = c(0, ymax[i]), axes = FALSE, ann = FALSE)
  box(lwd = 2)
  grid(col = gray(0.7), lwd = 1, lty = 2)
  lines(years[spring.ind], spring.biomass[spring.ind], lwd = 2)
  polygon(c(years[spring.ind],rev(years[spring.ind])), c(spring.biomass.lo[spring.ind],rev(spring.biomass.hi[spring.ind])), col = tcol, border = "transparent")
  axis(1, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5)
  mtext(side = 3, "Spring", line = 1, cex = 1.5)

  plot(years, fall.biomass, type = 'n', ylim = c(0, ymax[i]), axes = FALSE, ann = FALSE)
  box(lwd = 2)
  grid(col = gray(0.7), lwd = 1, lty = 2)
  lines(years[fall.ind], fall.biomass[fall.ind], lwd = 2)
  polygon(c(years[fall.ind],rev(years[fall.ind])), c(fall.biomass.lo[fall.ind],rev(fall.biomass.hi[fall.ind])), col = tcol, border = "transparent")
  axis(1, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5, labels = FALSE)
  mtext(side = 3, "Fall", line = 1, cex = 1.5)

  #stock.name.text = use.stock.names[i]
  #mtext(side = 4, stock.name.text, line = 3, cex = 1.5)
  #stop()
  mtext(side = 2, parse(text = paste0("Biomass",  "~(10^3", "~ mt)")), line = 2, cex = 1.5, outer = TRUE)
  mtext(side = 1, 'Year', line = 2, cex = 1.5, outer = TRUE)
  dev.off()

}
