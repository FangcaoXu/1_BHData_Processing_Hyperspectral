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



# Define a longitude and latitude
# Extract spectra for all files in BH relative to this location
# -- Read the navigation file to extract the elevation angle, altitude, extent of location
# -- If it is included, then read the BH geo file to extract long and lat and identify matrix location
# -- Read the spectra associated with the location in the matrix
library(tools)
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(rgl)
library(geosphere)
library(plotKML)
source("Z_plot_functions.R")
source("Z_general_functions.R")


# createresult directory
res.dir <- create.results.dir(base.dir = res.dirname, prefix=res.prefix)
result <- data.frame()
# Read data
crs.ll <- CRS("+proj=longlat")
dir.name <- "/amethyst/s0/fbx5002/NittanyRadiance2019/20190418/bh"
footprints.fname <- "~/disk10TB/DARPA/BH_analysis/Shapefiles/NittanyRadiancePolygons.shp"
res.dirname      <- "~/disk10TB/DARPA/BH_analysis"
res.prefix       <- "results_"


footprints.spdf <- readOGR(footprints.fname) # 1436*8
rownames(footprints.spdf@data) <- 1:nrow(footprints.spdf@data)
crs(footprints.spdf) <- crs.ll
# target
# 40.804733, -77.865312
target    <- cbind(-77.865312, 40.804733)
target.sp       <- SpatialPoints(target, proj4string = crs.ll)


# Extract all the scenes that include this target
# I assume only a single target at a time
o   <- over(target.sp, footprints.spdf,returnList = T)[[1]] 
sids <- as.numeric(rownames(o))

# expand the extent of the sp object 
ext <- extent(gBuffer(target.sp, width=.005))
#jpeg("./Footprints.jpg", width=750,height=500, res = 100)
plot(ext,col="NA", xlab="Longitude", ylab="Latitude")
plot(footprints.spdf,add=T,col="grey")
plot(footprints.spdf[sids,],col="blue",add=T)
plot(target.sp,add=T,col="red",lwd=5)
#dev.off()
footprints.spdf@data[sids[1:5],]


######################################################################################################################
######################################################################################################################
######################################################################################################################

######################################################################################################################
######################################################################################################################
######################################################################################################################
# HSI images
sids <- footprints.spdf@data$sid
npixels <- vector()
for ( i in 1:length(sids)) {
  print( i )
  # Identify the filename for the geo location information asssociate with this scene
  scene <- footprints.spdf@data[sids[i],]
  hdr.fname <- scene$filename
  geo.fname <- gsub("hdr$","geo",hdr.fname)
  img.fname <- gsub("hdr$","img",hdr.fname)
  jpg.fname <- gsub("hdr$","jpg",hdr.fname)
  
  # location and dem file for 94686 cells
  geo.full.fname <- paste(dir.name,"/BH_",substring(geo.fname,1,15),"/nav/",geo.fname,sep="")
  # hyperspectral image of 256 bands for 94686 cells 
  img.full.fname <- paste(dir.name,"/BH_",substring(geo.fname,1,15),"/hsi_lwir_",scene$side,"/cal/",img.fname,sep="")
  jpg.full.fname <- paste(res.dir,jpg.fname,sep="")
  # Read the geolocation ifnormation
  # dimensions: 367, 258, 94686, 3  (nrow, ncol, ncell, nlayers), here the nrow, ncol, ncell may vary for different images
  geo.ras <- brick(geo.full.fname)
  npixels <- c(npixels, ncell(geo.ras))
  # nlayers: represent for lon, lat, dem
  lon <- as.vector(geo.ras[[1]]) # as.matrix(geo.ras[[1]]) #(367 258)
  lat <- as.vector(geo.ras[[2]]) # lat[1] == 40.80743; FALSE
  dem <- as.vector(geo.ras[[3]])

  # Convert to Spatial Points
  geo.spdf <- SpatialPointsDataFrame(cbind(lon,lat), data=data.frame(dem=dem), proj4string = crs.ll)
  #plot(geo.spdf, col=val2col(unlist(geo.spdf@data), col=rainbow(70)))

  # Compute the distance between the target to all cell to find which cell contains the target
  dmat <- distm(geo.spdf, target.sp, fun=distGeo)    # dim: (94686,1)
  cellid   <- which.min(dmat) # this cell id is same fro all (367,258) data matrix
  coords <- geo.spdf@coords[cellid,]  # xyFromCell(geo.ras, cellid)

  # nlayers: represent for radiance observed at 256 wavelength bands
  img.ras <- brick(img.full.fname) # 367, 258, 94686, 256 (256 bands)
  target.spectra <- as.vector(extract(img.ras,cellid)) # extract value based on cell id, 256
  nrow(img.ras); ncol(img.ras); ncell(img.ras); nlayers(img.ras); names(img.ras)

  # import as image data for one specific wavelength
  jpeg(jpg.full.fname, width=nrow(img.ras)*2,height=ncol(img.ras)*2)
  # read all values for a specific wavelength
  # all(as.vector(img.ras[[100]])==values(img.ras[[100]]))
  BH.spdf <- SpatialPointsDataFrame(cbind(lon,lat), data=data.frame(BH=as.vector(img.ras[[100]])), proj4string = crs.ll)
  plot(ext,col="white", xlab="Longitude", ylab="Latitude")
  plot(BH.spdf, col=val2col(BH.spdf$BH, col=rev(brewer.pal(11,"Spectral"))), pch=20,add=T,cex=0.5)
  plot(target.sp, col="red",pch=1, cex=1,add=T)
  dev.off()
  writeOGR(BH.spdf, res.dir, file_path_sans_ext(jpg.fname), driver="ESRI Shapefile", overwrite_layer = T)
  # angle is elevation angle from the sensor looking down to the target
  df <- data.frame(time=scene$time, filename=scene$filename, angle=scene$angle, altitude=scene$altitude,
                   side=scene$side, cell.lat=lat[cellid], cell.lon=lon[cellid], cell.dem=dem[cellid], range=(scene$altitude-dem[cellid])/sin(abs(scene$angle)*pi/180),
                   spectra=t(target.spectra))

  result <- rbind(result, df)
}

# ffmpeg -framerate 5 -pattern_type glob -i '*.jpg' BHExtract.mp4
# # wavelengths are slightly different for each image
# # Adjust the column names to reflect the spectra
# names(result)[grep("spectra",names(result))] = names(img.ras)
# result.path <- paste(res.dir,"result.csv", sep="")
# write.csv(result, result.path)
# result <- read.csv(result.path)[,-1]

# # If you want to read exisiting result, do this
result<- read.csv('~/disk10TB/DARPA/BH_analysis/results_20200511_082121/result.csv')[,-1]
ids <- grep("microns",names(result)) # here the id is the column id of original result matrix
rad.result <- t(result[,ids]*10^-6)  #'rad.result' is the transposed matrix of the 'result' above
rad.geoinfo <- t(result[,-ids])      # first  column is first row of original result matrix
dim(rad.result); dim(rad.geoinfo)
wavelength <- as.numeric(gsub("X([0-9]+\\.[0-9]+).*$","\\1",rownames(rad.result)))
rownames(rad.result) <- wavelength

range(result$range) # (3354, 5297)
# let's first select radiance except for outliers
outlier.cols <- unique(unname(which(rad.result>=0.00095, arr.ind=TRUE)[,2]))
rad.result <- rad.result[,-outlier.cols]
rad.geoinfo <- rad.geoinfo[, -outlier.cols]
angles <- abs(as.numeric(rad.geoinfo["angle",]))
colnames(rad.result) <- angles
angles.sorted <- sort(angles, index.return = TRUE) # list, $x; $ix
rad.result <- rad.result[, angles.sorted$ix]

# First try, find nearest observed angles to simulated andgles and copy largest angles (55.9) to all missing angles after 55 until 90.
cols.selected <- unique(findnearest(angles.sorted$x, seq(30,60, 5)))
angles.selected <- angles.sorted$x[cols.selected]
write.csv(rad.result[, cols.selected], "BHExtracted.csv")


######################################################################################################################################
# plot all radiance lines including or excluding outliers
jpeg("./RadiaceExtractedAtSamplePoints.jpg", width = 800, height=600, res=100)
matplot(rad.result,type="l", xlab=expression(paste("Wavelength (", mu, "m)")), xaxt = 'n', 
        ylab=expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), col = rev(rainbow(ncol(rad.result))))
axis(1, at=seq(1,256, length.out=8), labels=round(seq(wavelength[1], wavelength[256], length.out=8), digits=2))
dev.off()

jpeg("./RadiaceExtractedAtSamplePoints1.jpg", width = 800, height=600, res=100)
matplot(rad.result, type="l", xlab=expression(paste("Wavelength (", mu, "m)")), xaxt = 'n', 
        ylab=expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")), ylim=c(0.0005,0.00095),col = rev(rainbow(ncol(rad.result))))
axis(1, at=seq(1,256, length.out=8), labels=round(seq(wavelength[1], wavelength[256], length.out=8), digits=2))
dev.off()



##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
################################################## Here is the Dividing Line #####################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
# Start to research the outliers
#### Select out those pink lines and their associated image ids
outliers.id <- unique(which(rad.result >=0.00097, arr.ind=TRUE)[,2])
outliers <- rad.result[,outliers.id]
outliers.col <- rev(rainbow(ncol(rad.result)))[outliers.id]
matplot(outliers, type="l", xlab=expression(paste("Wavelength (", mu, "m)")), xaxt = 'n', 
        ylab=expression(paste("Radiance (W cm-2 sr-1 ", mu, "m-1)")),col = outliers.col)
axis(1, at=seq(1,256, length.out=8), labels=round(seq(wavelength[1], wavelength[256], length.out=8), digits=2))
# get the outlier file names
outlier.files <- result$filename[outliers.id]
# get all files taking at the same time and also cover the target point with the outlier files
# grep: search for matches within each element of a character vector and return a vector of indices
allfiles.insamefolder <- result$filename[grep(substring(outlier.files,1,15)[1], result$filename)]
# check the visible original images
img.fname <- gsub("hdr$","img",outlier.files)
img.fname <- gsub("hsi_s[1,2]","vis",img.fname)
# find associated images in visible bands
img.full.fname <- paste(dir.name,"/BH_",substring(img.fname,1,15),"/hri_vis/cal/",img.fname,sep="")
brick(img.full.fname[1])
plotRGB(brick(img.full.fname[1]), stretch='hist')

# find the pixels similar to a specific radiance
# ids: columns
# water   <- as.numeric(result[1,ids])
# img.mat <- matrix(as.numeric(values(img.ras)),ncol=dim(img.ras)[3])
# f <- function(x, y) {
#      return( mean( (x-y)^2 ) )
#    }
# diff <- apply(img.mat, 1, f,as.numeric(water))
# diff <- sapply(img.mat[,150], f, as.numeric(water[150]))
# diff[diff>2]=NA
# plot(geo.spdf, col=val2col(diff, col=rainbow(100),na.rm = T),pch=19)





