#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
#
#Author: Guido Cervone (cervone@psu.edu)
#        Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#        Department of Geography and Institute for CyberScience
#        The Pennsylvania State University
#

library(raster)  
library(rgdal)
library(lubridate)
library(dplyr)
require('grid')
require('gridBase') 
library(mapview)

# generate pie for clocks
pie.clock <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
                       init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
                       col = NULL, border = NULL, lty = NULL, main = NULL, ...){
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1) * P$x, c(1, 1) * P$y)
      text(1.15 * P$x, 1.13 * P$y, labels[i], xpd = TRUE, ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}
# functions to draw the clock
drawClock <- function(hour, minute, second)
{
  pie.clock(rep(1,12),1:12,
      radius=1.0,
      init.angle=75,
      clockwise = T,
      border=NA,
      col="lightgrey",cex=2)
  t <- seq(0, 2*pi, length=13)[-13]
  x <- cos(t)
  y <- sin(t)
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  # ticks
  grid.segments(x, y, x*.9, y*.9, default="native",gp=gpar(lex=.5))
  # Hour hand
  hourAngle <- pi/2 - (hour + minute/60)/12*2*pi
  grid.segments(0, 0,
                .6*cos(hourAngle), .6*sin(hourAngle),
                default="native", gp=gpar(lex=3))
  # Minute hand
  minuteAngle <- pi/2 - (minute)/60*2*pi
  grid.segments(0, 0,
                .8*cos(minuteAngle), .8*sin(minuteAngle),
                default="native", gp=gpar(lex=2))   
  # Second hand
  secondAngle <- pi/2 - (second)/60*2*pi
  grid.segments(0, 0,
                .8*cos(secondAngle), .8*sin(secondAngle),
                default="native", gp=gpar(lex=1, col="red"))
  popViewport(3)
}
# function to smooth the time
# when the second changes from 59 to 0???
smoothclock.time <- function(clocktime){
  # first column, hour; Second column, minute; Third column, second
  unique.index <- which(!duplicated(clocktime))
  for(i in 1:(length(unique.index)-1)){
    duplicatelength <- unique.index[i+1]-unique.index[i]
    if(duplicatelength>1){
      
      starting <- clocktime[unique.index[i]]
      ending <- clocktime[unique.index[i+1]]
      
      for(j in 1:(duplicatelength-1)){ # skip the first element in the each duplicate group since it's the starting value
        clocktime[unique.index[i]+j]  <- clocktime[unique.index[i]+j-1] + (ending-starting)/duplicatelength
      }
    }
  }
  # find whether the last unique element has duplicates or not
  if(unique.index[length(unique.index)] < length(clocktime)){
    duplicatelength <- length(clocktime)-unique.index[length(unique.index)]
    for(k in 1:length(duplicatelength)){
      clocktime[unique.index[length(unique.index)]+k]=clocktime[unique.index[length(unique.index)]+k-1]+0.5
    }
  }
  return(clocktime)
}

topdir <- "~/disk10TB/DARPA/Mosaic_VIS/cragMosaics"
# "Vis_2019-04-18_17-40-34_03122.ntf"
filenames <- dir(topdir,recursive = T,pattern="ntf$",full.names = T)
# extract the time from all files
re <- regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}-[0-9]{2}-[0-9]{2}", filenames)
times <- regmatches(filenames, re)
times <- as.POSIXct(times,format="%Y-%m-%d_%H-%M-%S", tz="UTC")
times <- with_tz(times,"EST5EDT")
# reformat the second to the 3 digits millisecond
op <- options(digits.secs=3)
times <- smoothclock.time(times)

# generate data table for the time index
clock.time <- cbind(
  hours <- as.numeric(format(times,format="%H")),
  minutes <- as.numeric(format(times,format="%M")),
  seconds <- as.numeric(format(times,format="%OS")) )
colnames(clock.time) <- c("hours", "minutes", "seconds")



#################################################################################
#################################################################################
# # let's solve the ls[[7]] problem in the footprints.kml.RData
# read locations and their spatial polygon 
locations <- readOGR("~/disk10TB/DARPA/BH_analysis/Shapefiles/NittanyRadiancePoints.shp")
times.loc<-as.POSIXct(locations$time)
load("footprints.kml.RData")
xlim=c(-77.88, -77.85)
ylim=c(40.784, 40.810)

# remove the ring circle observed
i = 1107 #7, 1107
plot(ls[[i]])
sppdf <- ls[[i]]
coords <- click(n = 4)
# close the points to polygon
coords <- rbind(coords, coords[1, ])
spp.coords <- coords2Polygons(coords, ID = 'bbox')
crs(spp.coords) <- crs(sppdf)
sppts <- SpatialPoints(sppdf@polygons[[1]]@Polygons[[1]]@coords, proj4string = crs(sppdf))
plot(sppts)
plot(spp.coords, add = T)
sppts.crop <- sppts[which(over(sppts, spp.coords) == 1), ]
# plot cropped polygons
plot(sppdf)
plot(sppts, add = T, cex = 0.5, col = 'grey')
plot(spp.coords, add = T)
plot(sppts.crop, add = T, col = 'red')
spp.crop <- coords2Polygons(sppts.crop@coords, ID = 'cropped')
ls[[i]]@polygons[[1]]@Polygons[[1]]@coords <- spp.crop@polygons[[1]]@Polygons[[1]]@coords
plot(ls[[i]])


###############################################################################################
# full scene images and video
for (i in 1:length(clock.time)) {
  fname <- paste("~/disk10TB/DARPA/Mosaic_VIS/full_scene/Plot",sprintf("%06d",i),".png",sep="")
  print(fname)
  if (!file.exists(fname)) {
    temp <- brick(filenames[i])
    # find the closest time from Blue heron data 
    idx <- which.min(abs(difftime(times[i], times.loc, units="secs")))
    png(fname,width=1920,height=1080)
    layout(matrix(c(3,3,3,1,
                    3,3,3,2,
                    3,3,3,2,
                    3,3,3,3),ncol=4,byrow=T))

    par(mar=c(2,10,2,2))
    drawClock(hour = clock.time[i,1], minute = clock.time[i,2], second=clock.time[i,3])

    par(mar=c(0,7,0,7))
    plot(locations, col="grey",pch=19,cex=.5,xlim=xlim,ylim=ylim)
    plot(ls[[idx]],add=T,col="lightgrey")
    plot(locations[idx,],add=T,col="red",pch=19,cex=1.5)
    
    plotRGB(temp,stretch="lin")
    legend('topleft',legend="Nittany Radiance 2019 \n\n April 18, 2019",bty="n",cex=4)
    
    dev.off()
  }
}




### Zoom in scence
i=1
plotRGB(brick(filenames[i]), stretch="lin")
ext<- drawExtent(show=TRUE, col="red")
plotRGB(brick(filenames[i]), stretch="lin", ext=ext)

for (i in 10794:length(clock.time)) {
  fname <- paste("~/disk10TB/DARPA/Mosaic_VIS/zoom_scene/PlotSmall",sprintf("%06d",i),".png",sep="")
  print(fname)
  if (!file.exists(fname)) {
    temp <- brick(filenames[i])
    # find the closest time from Blue heron data 
    idx <- which.min(abs(difftime(times[i], times.loc, units="secs")))
    png(fname,width=1920,height=1080)
    layout(matrix(c(3,3,3,1,
                    3,3,3,2,
                    3,3,3,2,
                    3,3,3,3),ncol=4,byrow=T))
    
    par(mar=c(2,10,2,2))
    drawClock(hour = clock.time[i,1], minute = clock.time[i,2], second=clock.time[i,3])
    
    par(mar=c(0,7,0,7))
    plot(locations, col="grey",pch=19,cex=.5,xlim=xlim,ylim=ylim)
    plot(ls[[idx]],add=T,col="lightgrey")
    plot(locations[idx,],add=T,col="red",pch=19,cex=1.5)
    
    plotRGB(temp,stretch="lin",ext=ext)
    legend('topleft',legend="Nittany Radiance 2019 \n\n April 18, 2019",bty="n",cex=4)
    
    dev.off()
  }
}








