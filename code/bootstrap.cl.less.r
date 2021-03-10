bootstrap.cl.less = function(cl.less)
{
  #bootstrap is conditioned on the number of day and night tows.
  cl.less.day = cl.less[which(cl.less$day.night == "day"),]
  cl.less.night = cl.less[which(cl.less$day.night == "night"),]
  
  ids = unique(cl.less.day$id)
  ids.boot = sample(ids, replace = TRUE)
  n.per.sta.boot = sapply(ids.boot, function(x) sum(cl.less.day$id == x))
  id.boot = rep(1:length(ids.boot), n.per.sta.boot)
  index.boot = unlist(sapply(ids.boot, function(x) which(cl.less.day$id == x)))
  cl.less.day = cl.less.day[index.boot,]
  cl.less.day$id = id.boot
  
  ids = unique(cl.less.night$id)
  ids.boot = sample(ids, replace = TRUE)
  n.per.sta.boot = sapply(ids.boot, function(x) sum(cl.less.night$id == x))
  id.boot = rep(1:length(ids.boot), n.per.sta.boot) + length(unique(cl.less.day$id))
  index.boot = unlist(sapply(ids.boot, function(x) which(cl.less.night$id == x)))
  cl.less.night = cl.less.night[index.boot,]
  cl.less.night$id = id.boot
  
  cl.less = rbind(cl.less.day, cl.less.night)
  return(cl.less)

}
