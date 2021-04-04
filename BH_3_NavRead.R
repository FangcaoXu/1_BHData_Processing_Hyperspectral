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

### PowerPC: big endian
### x86: little endian
library(raster)
library(mcga)
library(rgdal)
library(sp)


######################################################################################################################
# Read given packets from Hexadecimal to double
NavReadDouble <- function(s, bin, interval=7) {
  packets.bin <- (bin[(s):(s+interval)])
  # Convert to Decimal
  packets.dec <- as.numeric(paste("0x",packets.bin,sep=""))
  # Convert to double
  double <- BytesToDouble(rev(packets.dec)) 
  return(double)
}
# info of the target, bits starting at 0
# range 184-191; heading 192-199; elevation 200-207
# Lat 208-215; Lon 216-223; Altitude 224 - 231
NavGround <- function(s, bin){
  # trackid <- sapply(s+176, NavReadDouble, bin)
  range <- sapply(s+184, NavReadDouble, bin)
  heading <- sapply(s+192, NavReadDouble, bin)
  ele <- sapply(s+200, NavReadDouble, bin)
  lats <- sapply(s+208, NavReadDouble, bin)
  longs <- sapply(s+216, NavReadDouble, bin)
  alt <- sapply(s+224, NavReadDouble, bin)
  ground.df <- data.frame(range, heading, ele, longs, lats, alt)
  return(ground.df)
}
# sigma is the error of measured variables
# sigma Lat 232 - 239; sigma Long 240 - 247
# sigma Down 248 - 255; sigma Range 256 - 263
NavSigma <- function(s, bin){
  lats.sigma <- sapply(s+232, NavReadDouble, bin)
  longs.sigma <- sapply(s+240, NavReadDouble, bin)
  down.sigma <- sapply(s+248, NavReadDouble, bin)
  range.sigma <- sapply(s+256, NavReadDouble, bin)
  sigma.df <- data.frame(longs.sigma, lats.sigma, down.sigma, range.sigma)
  return(sigma.df)
}
# GPS Gimbal: info of the plane/gimbal on the ground
# gps time: 284 - 291 gps week: 292-297
# Gmb Lat: 300 - 307; Gmb Long: 308 - 315; Gmb Altitude: 316 - 323
# Gmb Elevation 324 - 331; Gmb Roll 332 -339; Gmb Heading 340 - 347
NavGmb <- function(s, bin){
  gpstime  <- sapply(s+284, NavReadDouble, bin) 
  gpsweek <- sapply(s+292, NavReadDouble, bin) 
  time <- as.POSIXct((gpsweek + 1024)*3600*24*7 + gpstime, origin="1980-01-06")
  lats.gmb <- sapply(s+300, NavReadDouble, bin)
  longs.gmb <- sapply(s+308, NavReadDouble, bin)
  alt.gmb <- sapply(s+316, NavReadDouble, bin)
  ele.gmb <- sapply(s+324, NavReadDouble, bin)
  roll.gmb <- sapply(s+332, NavReadDouble, bin)
  heading.gmb <- sapply(s+340, NavReadDouble, bin)
  gmb.df <- data.frame(time=time, longs.gmb, lats.gmb, alt.gmb, ele.gmb, roll.gmb, heading.gmb)
  return(gmb.df)
}
# combine three functions above and return designed value 
readNav <- function(turret.file, error.return = FALSE){
  bin <-  readBin(turret.file, n=10000000, what="raw", endian='big')
  # Packets start wit 0x00 0x3c 0x03 0x38 0x6a 0x6a
  szz <- which(bin==0x00)
  start <- szz[which(bin[szz+1]== 0x3C & bin[szz+2]== 0x03 & bin[szz+3]== 0x38 & bin[szz+4]== 0x6a & bin[szz+5]==0x6a)]
  # Define the interval_bit
  interval_bits <- start[-1] - start[1:(length(start)-1)]
  if (any(interval_bits!= 824)) { # bits for each packet will be 824
    stop("The packet sizes are different")
  }
  # read plane and target info given the turrent file
  ground.info <- NavGround(start, bin)
  gmb.info <- NavGmb(start, bin)
  if(error.return){
    errors.info <- NavSigma(start, bin)
    return(cbind(gmb.info, ground.info, errors.info))
  }
  else{
    return(cbind(gmb.info, ground.info))
  }
}
# read all turret files in the same folder
turretReader_per_folder <- function(folderid, id_turret){
  turretinfo <- NULL
  tids <-  which(id_turret$folderid == folderid)
  for (tid in tids){
    info.nav <- readNav(id_turret$turret.all[tid])
    turret.info <- cbind(tid, folderid, foldername=id_turret$foldername[tid], info.nav)
    turretinfo <- rbind(turretinfo, turret.info)
  }
  return(turretinfo)
}
# decide the timelag between turret files and scene files
timelagcal <- function(turretinfo, sceneinfo, id=c(1,2)){
  # layout can be distributed. Then using the second scene is better
  # layout is crowd, using first scene is better
  targetsp <- SpatialPoints(turretinfo[, c("longs", "lats")], proj4string = crs.ll)
  for(i in id){
    over_twoside <- which(!is.na(over(targetsp, footprints.spdf[sceneinfo$sid[c(i, i+nrow(sceneinfo)/2)],])[1]))
    timediff1 <-  sceneinfo$time[i] - range(turretinfo[over_twoside,]$time)[1]     # start overlap time
    timediff2 <-  sceneinfo$time[i+1] - range(turretinfo[over_twoside,]$time)[2]   # end overlap time
    if (i == 1){                               # when scene 1 and 2 are spatially continous with each other
      if(abs(timediff2) > 1.2*abs(timediff1)){ # points overlap with this scene during other scenes collection
        timelag <- timediff1                   # use start
      }
      else{
        timelag <- timediff2                   # use end as default, 
      }
    }
    else{                                      # i = 2, when scene 1 and 2 are not spatially continous 
      timelag <- timediff1                     # use start as default for scene 2 to calculate time lag
    }
    # this relationship is roughly true for scenes in all the folders
    if(timelag < -18 & timelag > -20){
      return(timelag)
    }
  }
  # assign it to be -18, with spatial overlay and continous time points check, it's accurate enough 
  return(-18)
}
# decide the extent for plot
decideExtent <- function(turretinfo, scene.extent){
  plane.extent <- extent(min(turretinfo$longs.gmb), max(turretinfo$longs.gmb), 
                         min(turretinfo$lats.gmb),  max(turretinfo$lats.gmb))
  # scene_extent out of plane extent or not
  extent.logi <- (scene.extent[1] < plane.extent[1]) | (scene.extent[2] > plane.extent[2]) | 
    (scene.extent[3] < plane.extent[3]) | (scene.extent[4] > plane.extent[4])
  if(extent.logi){
    ext <- extent(min(scene.extent[1], plane.extent[1]), max(scene.extent[2], plane.extent[2]), 
                  min(scene.extent[3], plane.extent[3]),  max(scene.extent[4], plane.extent[4]))
  }else{ext <- plane.extent}
  return(ext)
}
# write to shapefile
gc.writeOGR <- function(data, fname, driver="ESRI Shapefile", ...) {
  writeOGR(data, dirname(fname), file_path_sans_ext(basename(fname)), driver=driver , ...)
}
# return nav points over one scene (idx) in one folder
return_point_perscene_perfolder <- function(points_scenes, idx, sceneinfo, turretinfo) {
  if (length(points_scenes[[idx]]) == 0){
    # write one line of NA with sid and folder id
    tmp <- cbind(sid=sceneinfo$sid[idx], sid_2=sceneinfo$sid[idx+nrow(sceneinfo)/2], turretinfo[points_scenes[[idx]],][1,])
    tmp$folderid <- sceneinfo$folderid[1]
    tmp$foldername <- sceneinfo$foldername[1]
    return(tmp)
  }else{
    # if time difference of two continous turret points exceeds 10 times of average, trim them out
    tmp <- turretinfo[points_scenes[[idx]],]
    aver.time <- (tmp$time[nrow(tmp)] - tmp$time[1])/(nrow(tmp)-1)
    for(i in 2:nrow(tmp)){
      if(tmp$time[i]-tmp$time[i-1] > 7*aver.time){
        t_break = i
        break
      }else{
        t_break = 1
      }
    }
    if(t_break > nrow(tmp)/2) tmp <- tmp[1:(t_break-1),] else tmp <- tmp[t_break:nrow(tmp),]
    return(cbind(sid=sceneinfo[idx, ]$sid, sid_2=sceneinfo[idx+nrow(sceneinfo)/2, ]$sid, tmp))
  }
}
# write all extracted scene_turret_folder_info to csv
return_scene_turret_folder_info  <- function(folderid, id_turret, id_scene){
  # target info
  turretinfo <- turretReader_per_folder(folderid, id_turret)
  targetsp <- SpatialPoints(turretinfo[, c("longs", "lats")], proj4string = crs.ll)
  # scene info
  sids <- which(id_scene$folderid == folderid) # scene id
  sceneinfo <- id_scene[sids, ]
  # time lag
  timelag <- timelagcal(turretinfo, sceneinfo)
  # return points by spatial overlay and time lag
  points_scenes = list()
  for(i in 1:(nrow(sceneinfo)/2)){
    if(i != nrow(sceneinfo)/2){
      turrettemp <- which((turretinfo$time + timelag) >= sceneinfo$time[i] & (turretinfo$time+timelag) < sceneinfo$time[i+1])
    }else{
      turrettemp <- which((turretinfo$time + timelag) >= sceneinfo$time[i])
    }
    # find continous points overlap with the scene
    over_twoside <- which(!is.na(over(targetsp, footprints.spdf[sceneinfo$sid[c(i, i+nrow(sceneinfo)/2)],])[1]))
    points_scenes[[i]] <- intersect(over_twoside,turrettemp)
  }
  # calculate for each scene in the current folder 
  scene_turretinfo <- lapply(1:length(points_scenes), function(i) return_point_perscene_perfolder(points_scenes, i, sceneinfo, turretinfo))
  scene_turretinfo <- do.call(rbind, scene_turretinfo)
  return(list(timelag, scene_turretinfo))
}

######################################################################################################################
######################################################################################################################
######################################################################################################################
########################################## Codes processing the data start here ######################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
# define the projection
crs.ll <- CRS("+proj=longlat")
# Read all the data
dir.name <- "/amethyst/s0/fbx5002/NittanyRadiance2019/20190418/bh"

# read all nav turret files which represent for plane geographic information and target geographic information
turret.all <- dir(path=dir.name, pattern="turret.*img$",recursive = T, include.dirs =  T, full.names=T)
footprints.fname <- "~/disk10TB/DARPA/BH_analysis/Shapefiles/NittanyRadiancePolygons.shp"
footprints.spdf <- readOGR(footprints.fname) # 1436
crs(footprints.spdf) <- crs.ll
rownames(footprints.spdf@data) <- 1:nrow(footprints.spdf@data)
scene.all <- footprints.spdf@data
scene.all$time <- as.POSIXct(scene.all$time)
length(turret.all); length(scene.all$filename) # 85 turret img in nav folder

# extract the folder pattern for turret
turret.pattern <- gsub(".*([0-9]{8}_[0-9]{6}).*", "\\1", turret.all)
# some folders doesn't have any hdr/geo/img in their lwir folders
folder.pattern <- regmatches(list.files(dir.name), regexpr("[0-9]{8}_[0-9]{6}", list.files(dir.name))) #65
scene.pattern <- regmatches(scene.all$filename, regexpr("[0-9]{8}_[0-9]{6}", scene.all$filename))      #1436
length(unique(scene.pattern))    # 55 unique
# relationship between hyperspectral with turret file based on pattern, id is the folder id
id_turret <- as.data.frame(cbind(tid = 1:length(turret.all), folderid = match(turret.pattern, folder.pattern), 
                                 foldername=folder.pattern[match(turret.pattern, folder.pattern)], turret.all), stringsAsFactors = FALSE)
id_scene <-  as.data.frame(cbind(sid=scene.all[,1], folderid = match(scene.pattern, folder.pattern), 
                                 foldername=folder.pattern[match(scene.pattern, folder.pattern)],
                                 scene.all[,2:ncol(scene.all)]), stringsAsFactors = FALSE)
# some folders only have turret files but no hyperspectral files
valid_folderids <- unique(id_scene$folderid)

######################################################################################################################
######################################################################################################################
######################################################################################################################
# # Observation 1: hyperspectral images of side 1 and 2 are taken at the same time
# sceneinfo$time[1:(nrow(sceneinfo)/2)] == sceneinfo$time[(nrow(sceneinfo)/2+1): nrow(sceneinfo)]
# # Observation 2: turret time slight changes even some look the same
# turretinfo$time[3] == turretinfo$time[2]; turretinfo$time[4] - turretinfo$time[3]
# # Observation 3: the turretinfo lag behind the sceneinfo
######################################################################################################################
################################################# Plot All Scenes and Nav ############################################
######################################################################################################################
# sceneinfo and turretinfo are info extracted for the files in a specific same folder
for (folderid in valid_folderids){
  turretinfo <- turretReader_per_folder(folderid, id_turret)
  sids <- which(id_scene$folderid == folderid)
  sceneinfo <- id_scene[sids, ]
  # decide the extent for plot
  ext <- decideExtent(turretinfo, extent(footprints.spdf[sids,]))
  # plot the scenes and turret together
  plot(ext, xlab="Longitude", ylab="Latitude", col="white")
  plot(footprints.spdf[sids,], col="grey", add =T)
  points(turretinfo[, c("longs", "lats")], col="red", cex = 0.3) # targets on the ground
  points(turretinfo[, c("longs.gmb", "lats.gmb")], col='blue', cex = 0.3)  # plane
}
############################# 
jpeg(file="navHyperspectral.jpeg", width = 3000, height = 5500, res=300)
par(mfrow=c(11,5), mar=c(1,1.5,1,1), oma = c(0.5, 1, 0.5, 0.5))
for(folderid in valid_folderids){ # 55: folders having valid nav and scenes
  turretinfo <- turretReader_per_folder(folderid, id_turret)
  sids <- which(id_scene$folderid == folderid)
  sceneinfo <- id_scene[sids, ]
  # decide the extent for plot
  ext <- decideExtent(turretinfo, extent(footprints.spdf[sids,]))
  x1 = ceiling(ext[1]*100)/100; x2 = floor(ext[2]*100)/100
  y1 = ceiling(ext[3]*100)/100; y2 = floor(ext[4]*100)/100
  # -77.82:6 chars; -77.825:7 chars
  if(nchar((x1 + x2)/2) <= 6) xlab <- c(x1, (x1 + x2)/2, x2) else xlab <- c(x1, x2)
  ylab <- c(y1, y2)
  # plot the scenes and turret together
  plot(ext, xlab="Longitude", ylab="Latitude", col="white", xaxt = "n", yaxt = "n")
  axis(1, at= xlab, labels=xlab, padj = -0.7)  
  axis(2, at= ylab, labels=ylab, padj = 0.7)  
  plot(footprints.spdf[sids,], col="grey", add =T)
  points(turretinfo[, c("longs", "lats")], col="red", cex = 0.3) # targets on the ground
  points(turretinfo[, c("longs.gmb", "lats.gmb")], col='blue', cex = 0.3)  # plane
}
dev.off()

###################################################################################################################
################################################# Dynamic Plot Each Folder ########################################
###################################################################################################################
# 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13  okay
folderid <-valid_folderids[55]
turretinfo <- turretReader_per_folder(folderid, id_turret)
sids <- which(id_scene$folderid == folderid)
sceneinfo <- id_scene[sids, ]
ext <- decideExtent(turretinfo, extent(footprints.spdf[sids,]))
targetsp <- SpatialPoints(turretinfo[, c("longs", "lats")], proj4string = crs.ll)
timelag <- timelagcal(turretinfo, sceneinfo)

timelag <- -18
#############################
# plot the scenes and turret together
plot(ext, xlab="Longitude", ylab="Latitude", col="white")
plot(footprints.spdf[sids,], col="grey", add =T)
points(turretinfo[, c("longs", "lats")], col="red", cex = 0.3) # targets on the ground
points(turretinfo[, c("longs.gmb", "lats.gmb")], col='blue', cex = 0.3)  # plane
#############################
# dynamic plotting 
plot(ext, xlab="Longitude", ylab="Latitude", col="white")
for(i in 1:(nrow(sceneinfo)/2)){
  if(i != nrow(sceneinfo)/2){
    turrettemp <- which((turretinfo$time + timelag) >= sceneinfo$time[i] & (turretinfo$time+timelag) < sceneinfo$time[i+1])
  }else{
    turrettemp <- which((turretinfo$time + timelag) >= sceneinfo$time[i])
  }
  if(i == 1 & turrettemp[1] > 1){
    # first plot starting points which are not over the scene
    points(turretinfo[1:(turrettemp[1]-1), c("longs", "lats")], col="red", cex = 0.3)
    points(turretinfo[1:(turrettemp[1]-1), c("longs.gmb", "lats.gmb")], col='blue', cex = 0.3)
    Sys.sleep(1)
  }
  # then start to plot each scene and its associated points
  plot(footprints.spdf[sceneinfo$sid[c(i,i+nrow(sceneinfo)/2)],], col="grey", add =T)    #side 1/2
  points(turretinfo[turrettemp, c("longs", "lats")], col="red", cex = 0.3)
  points(turretinfo[turrettemp, c("longs.gmb", "lats.gmb")], col='blue', cex = 0.3)
  Sys.sleep(1)
}
###############################
# correspond turret points to their correct scenes
points_scenes = list()
for(i in 1:(nrow(sceneinfo)/2)){
  if(i != nrow(sceneinfo)/2){
    turrettemp <- which((turretinfo$time + timelag) >= sceneinfo$time[i] & (turretinfo$time+timelag) < sceneinfo$time[i+1])
  }else{
    turrettemp <- which((turretinfo$time + timelag) >= sceneinfo$time[i])
  }
  # find continous points overlap with the scene
  over_twoside <- which(!is.na(over(targetsp, footprints.spdf[sceneinfo$sid[c(i, i+nrow(sceneinfo)/2)],])[1]))
  points_scenes[[i]] <- intersect(over_twoside,turrettemp)
}
# calculate for each scene in the current folder 
length(points_scenes) == nrow(sceneinfo)/2
scene_turretinfo <- lapply(1:length(points_scenes), function(i) 
  return_point_perscene_perfolder(points_scenes, i, sceneinfo, turretinfo))
scene_turretinfo <- do.call(rbind, scene_turretinfo)

################################
# check the continouty of points generated based on 4_BH_StatisticNavwithScene for folderid 8
tmp <- turretinfo[points_scenes[[4]],]
aver.time <- (tmp$time[nrow(tmp)] - tmp$time[1])/(nrow(tmp)-1)
for(i in 2:nrow(tmp)){
  if(tmp$time[i]-tmp$time[i-1] > 5*aver.time){
    t_break = i
    break
  }else{
    t_break = 0
  }
}
tmp$time[t_break] - tmp$time[t_break-1]
# examine overlay and timelag relationship
# only depend on overlay, points identified for a specific scene is not correct
# only depends on timelag is not accurate for all the scene
i = 1 # here i are the scene index
turrettemp <- which((turretinfo$time + timelag) >= sceneinfo$time[i] & (turretinfo$time+timelag) < sceneinfo$time[i+1])
over_twoside <-  which(!is.na(over(targetsp, footprints.spdf[sceneinfo$sid[c(i, i+nrow(sceneinfo)/2)],])[1]))
par(mfrow = c(1,2))
plot(ext, xlab="Longitude", ylab="Latitude", col="white")
plot(footprints.spdf[sceneinfo$sid[c(i,i+nrow(sceneinfo)/2)],], col="grey", add=T)
points(turretinfo[turrettemp, c("longs", "lats")], col="red", cex = 0.3)
points(turretinfo[turrettemp, c("longs.gmb", "lats.gmb")], col='blue', cex = 0.3)
plot(ext, xlab="Longitude", ylab="Latitude", col="white")
plot(footprints.spdf[sceneinfo$sid[c(i,i+nrow(sceneinfo)/2)],], col="grey", add=T)
points(turretinfo[over_twoside, c("longs", "lats")], col="red", cex = 0.3)
points(turretinfo[over_twoside, c("longs.gmb", "lats.gmb")], col='blue', cex = 0.3)
dev.off()
# depends on both time lag and overlay
plot(footprints.spdf[sceneinfo$sid[c(1, 1+nrow(sceneinfo)/2)],], col="grey")
points(turretinfo[points_scenes[[1]], c("longs", "lats")], col="red", cex = 0.3)
plot(footprints.spdf[sceneinfo$sid[c(2, 2+nrow(sceneinfo)/2)],], col="grey")
points(turretinfo[points_scenes[[2]], c("longs", "lats")], col="red", cex = 0.3)
# check index between two adjacent scenes
plot(footprints.spdf[sceneinfo$sid[c(i:(i+2), i:(i+2)+nrow(sceneinfo)/2)],], col="grey")
points(turretinfo[points_scenes[[1]], c("longs", "lats")], col="red", cex = 0.3)
start <- points_scenes[[1]][length(points_scenes[[1]])]+1
end <-  points_scenes[[2]][1]-1
points(turretinfo[start:end, c("longs", "lats")], col="red", cex = 0.3)
points(turretinfo[over_twoside, c("longs", "lats")], col="blue", cex = 0.3)

###################################################################################################################
##################################### Write Extracted Nav, Scene, Folder info to CSV ##############################
###################################################################################################################
timelag_info <- c()
scene_turret_folder_info <- NULL
for(i in valid_folderids){
  print(i)
  tmp <- return_scene_turret_folder_info(i, id_turret, id_scene)
  timelag_info <- c(timelag_info, tmp[[1]])
  scene_turret_folder_info <-rbind(scene_turret_folder_info, tmp[[2]])
}

write(timelag_info, '~/disk10TB/DARPA/BH_analysis/scene_turret_folder_timelag.txt', ncolumns =1)
save(id_scene, file='~/disk10TB/DARPA/BH_analysis/id_scene.RData')
save(id_turret, file='~/disk10TB/DARPA/BH_analysis/id_turret.RData')
save(scene_turret_folder_info, file='~/disk10TB/DARPA/BH_analysis/scene_turret_folder_info.RData')





######################################################################################################################
######################################################################################################################
######################################################################################################################
# Run codes below only for checing accuracy of measured info given by nav turret files
# take one as example to check the range calculated based on ground points and range given by gimbal info
info.nav <- readNav(id_turret$turret.all[10])
class(info.nav[, 11:12])
# calculate range based on horizontal distance and altitudes
cal.ranges <- sqrt((info.nav[, 4] - info.nav[,13])^2 +  pointDistance(info.nav[, 2:3], info.nav[, 11:12], lonlat= TRUE)^2)
range(cal.ranges)
range(info.nav[, 8])
#########################################################################################################################
# create spatial points 
coords.all <- rbind(info.nav[, 2:3], unname(info.nav[, 11:12]))
sp <- SpatialPoints(coords.all, proj4string = crs.ll)
spdf <- SpatialPointsDataFrame(sp, data=data.frame(val=coords.all))
# gc.writeOGR(spdf,"targets.shp")
# create soatial lines
lines=list()
for ( i in 1:nrow(coords)) {
  r <- spRbind(spdf[i,],spdf[i+nrow(coords),])
  rl <- as(r,"SpatialLines")
  lines[[i]] <- rl  
}
merged.lines <- do.call(rbind, lines)
lines.spdf <- SpatialLinesDataFrame(merged.lines,data=data.frame(id=1:length(merged.lines)),match.ID = F)
#gc.writeOGR(lines.sldf,"lines.shp")