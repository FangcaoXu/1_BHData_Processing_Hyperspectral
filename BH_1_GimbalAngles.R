#
#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
#
# Author: Guido Cervone (cervone@psu.edu) and Fangcao Xu (xfangcao@psu.edu)
#         Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#         Department of Geography and Institute for CyberScience
#         The Pennsylvania State University
#

# EST5EDT: given the time, decide whether its
# This time zone changed to daylight saving time at 2:00 AM on Mar 8, 2020. The GMT offset is UTC/GMT -4 hours (EDT)
# It will change back to standard time at 2:00 AM on Sunday, Nov 1, 2020. The GMT offset is UTC/GMT -5 hours (EST)

# read hdr files in lwir_1/cal; lwir_2/cal folders
library(lubridate)
library(sp)
library(rgdal)
library(raster)
source("Z_plot_functions.R")


my.crs <- CRS("+proj=longlat")
dir.name <- "/amethyst/s0/fbx5002/NittanyRadiance2019/20190418/bh"
# dir.name <- "/Volumes/Data8TB/NittanyRadiance2019/20190418/bh/"
# read all hdr files
hdr.all <- dir(path=dir.name,pattern=".hdr$",recursive = T,include.dirs =  T,full.names=T)
hdr.lwir <- hdr.all[grep("lwir",hdr.all) ]
hdr.lwir <- hdr.lwir[ grep("cal",hdr.lwir) ]
length(hdr.lwir)
# remove .hdr in such the format that likes 20190418_182409_basic_badpixelmap_0000.hdr
filenames <- hdr.lwir[grep("basic",hdr.lwir,invert=T) ]

mat <- matrix(0, ncol=19, nrow=length(filenames))
for ( f in 1:length(filenames) ) {
  print(f/length(filenames)*100)
  data <- scan(filenames[f],what="Character",sep="\n")
  
  # time info
  line.id <- grep("unixtime",data)
  temp <- data[line.id] 
  re <-  regexpr("=?[0-9]+",temp)
  time <- as.numeric(regmatches(temp, re))/1e+9

  # extract geometric info
  line.id <- grep("Angle", data)
  temp <- data[line.id]  
  re <- regexpr("-?[0-9]{2}\\.[0-9]{2}",temp)
  angle <- as.numeric(regmatches(temp, re))
  # altitude
  temp <- data[line.id-1]  
  re <- regexpr("-?[0-9]{1,4}\\.?[0-9]*",temp)
  altitude <- as.numeric( regmatches(temp, re))
  # longitude and latitude at oimage four cornors
  temp <- data[line.id+1]  
  re <- regexpr("[0-9]+,\\s+[0-9]+,\\s+-?[0-9]{1,2}\\.[0-9]+,\\s+-?[0-9]{1,3}\\.[0-9]+",temp)
  v1 <- as.numeric( strsplit(regmatches(temp,re),",")[[1]])
  ##
  temp <- data[line.id+2]  
  re <- regexpr("[0-9]+,\\s+[0-9]+,\\s+-?[0-9]{1,2}\\.[0-9]+,\\s+-?[0-9]{1,3}\\.[0-9]+",temp)
  v2 <- as.numeric( strsplit(regmatches(temp,re),",")[[1]])
  ##
  temp <- data[line.id+3]  
  re <- regexpr("[0-9]+,\\s+[0-9]+,\\s+-?[0-9]{1,2}\\.[0-9]+,\\s+-?[0-9]{1,3}\\.[0-9]+",temp)
  v3 <- as.numeric( strsplit(regmatches(temp,re),",")[[1]])
  ##
  temp <- data[line.id+4]  
  re <- regexpr("[0-9]+,\\s+[0-9]+,\\s+-?[0-9]{1,2}\\.[0-9]+,\\s+-?[0-9]{1,3}\\.[0-9]+",temp)
  v4 <- as.numeric(strsplit(regmatches(temp,re),",")[[1]])
  ##
  mat[f,] <- c(time, angle,altitude,v1,v2,v3,v4)
}

# as.POSIXct(1555609137436638900/1e+9, origin = '1970-01-01', tz="EST5EDT")
### combine all information together
lon.id <- c(7,11,15,19)
# longitude and latitude are the mean values of the four cornors of the image
df <- data.frame(filename=basename(filenames), time = as.POSIXct(mat[,1], origin = '1970-01-01', tz="EST5EDT"), 
                 angle=mat[,2], altitude=mat[,3], side=1, lon=apply(mat[,lon.id],1,mean), 
                 lat=apply(mat[,(lon.id-1)],1,mean))
# get the index of the images took by the second camera
side2 <- grep("lwir_2",filenames)
df$side[side2] <- 2  # Determine which angles are relative to side2
df$index <- 1:nrow(df)

# Spatial Points
pts <- SpatialPointsDataFrame(cbind(df$lon, df$lat), data=df,proj4string =my.crs)
writeOGR(pts, "~/disk10TB/DARPA/BH_analysis/Shapefiles", "NittanyRadiancePoints", driver="ESRI Shapefile", overwrite_layer = T) 
# Spatial Polygons
ls <- list()
for ( i in 1:nrow(df) ) {
  coord <- matrix(mat[i, c(6,5,10,9,14,13,18,17)+1],ncol=2,byrow=T)
  p  <-  Polygon( coord)
  ps <-  Polygons(list(p), i) # i is the id
  ls[[i]] <- ps
}
ps.sp <- SpatialPolygons(ls, proj4string =my.crs)
ps.spdf <- SpatialPolygonsDataFrame(ps.sp, data=df)
writeOGR(ps.spdf, "~/disk10TB/DARPA/BH_analysis/Shapefiles", "NittanyRadiancePolygons", driver="ESRI Shapefile",overwrite_layer = T) 

# plot the polygons created from the hdr files
valid <- df$side==1
plot(df$time[valid], df$angle[valid])
col <- rainbow(100)
plot(df$lon[valid], df$lat[valid], col=val2col(df$angle[valid], min=floor(min(df$angle[valid])), max=ceiling(max(df$angle[valid])), col=col),
     pch=19, xlim=c(-77.89, -77.84), ylim=c(40.78, 40.815))
points(df$lon[!valid], df$lat[!valid], col=val2col(df$angle[!valid], min=floor(min(df$angle[valid])), 
                                                   max=ceiling(max(df$angle[valid])), col=col),pch=18,cex=1)


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
################################################## Here is the Dividing Line #####################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

# Now read the polygons associated with BH collections
# "nav" folder: georeferencing products generated for each calibrated HSI scene
# nav: 1436 files
# hdr: lwir, cal, remove basic_bad_pixel, 1436 files.
kml.all <- dir(path=dir.name,pattern="\\.kml$",recursive = T,include.dirs =  T,full.names=T)
kml.nav <- kml.all[grep("nav",kml.all)]
re <- regexpr("[0-9]{8}_[0-9]{6}", kml.nav)
times.nav <- regmatches(kml.nav, re)
times.nav <- as.POSIXct(times.nav, format="%Y%m%d_%H%M%S", tz="UTC")
times.nav <- with_tz(times.nav,"EST5EDT")
df.nav <- data.frame(time=times.nav, filename=basename(kml.nav))
# represent for the same file in the same sequence
all(sub('\\..*$', '', basename(kml.nav))==sub('\\..*$', '', df$filename))

# read each kml file
ls.nav <- list()
for ( f in 1:length(kml.nav) ) {
  print(f/length(kml.nav)*100)
  ls.nav[[f]] <- readOGR(kml.nav[f])
}

# ls[[]]: An object of class "Polygons"
# ps.spdf: SpatialPolygonsDataFrame
# ps.spdf@polygons: List of Polygons
# ps.spdf@polygons[[]]: Polygons
# ls.nav[[]]: SpatialPolygonsDataFrame
# ls.nav[[1]]@polygons: one element list
# ls.nav[[1]]@polygons[[1]]: Polygons
# ls.nav[[1]]@polygons[[1]]@Polygons[[1]] : Polygon
# a list of Polygons
temp.nav <- lapply(ls.nav, function(x){x@polygons[[1]]@Polygons[[1]]})
polygons.nav <-  lapply(1:nrow(df), function(x) Polygons(list(temp.nav[[x]]), x))
sp.nav <- SpatialPolygons(polygons.nav, proj4string =my.crs)
spdf.nav <- SpatialPolygonsDataFrame(sp.nav, data=df) # after checking the kml and hdr are representing for the same file

plot(spdf.nav)
plot(ls.nav[[1]], xlim=c(-77.92167, -77.81421), ylim=c(40.66608, 40.80969))
lapply(ls.nav, function(x) plot(x, add=T))
############### Conclusion, some KML files have the data problem, which is shown in the plots above
