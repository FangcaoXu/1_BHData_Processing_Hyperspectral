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

library(data.table)
library(geosphere)
library(raster)
library(rgdal)
library(sp)
source("Z_BH_functions.R")

######################################################################################################################
########################################## Codes processing the data start here ######################################
######################################################################################################################
crs.ll <- CRS("+proj=longlat")
dir.name <- "/amethyst/s0/fbx5002/NittanyRadiance2019/20190418/bh"
# Read scene turret info
load('~/disk10TB/DARPA/BH_analysis/scene_turret_folder_info.RData')
load('~/disk10TB/DARPA/BH_analysis/id_scene.RData')
load('~/disk10TB/DARPA/BH_analysis/id_turret.RData')
# Here the scene_turret_folder_info sid only contain scenes at side 1, however, turret points correspond to the scenes at two sides.
length(unique(id_scene$sid)); length(unique(scene_turret_folder_info$sid))
# Read footprints
footprints.fname <- "~/disk10TB/DARPA/BH_analysis/Shapefiles/NittanyRadiancePolygons.shp"
footprints.spdf <- readOGR(footprints.fname) # 1436
rownames(footprints.spdf@data) <- 1:nrow(footprints.spdf@data)
crs(footprints.spdf) <- crs.ll
# group by scene id
sceneturretinfo <- scene_turret_folder_info[!is.na(scene_turret_folder_info$tid),]
by_folder <- split(sceneturretinfo, sceneturretinfo$folderid)
folder_ids <- unique(sceneturretinfo$folderid)
################################
############### write points and center points in each folder to shapefiles
# res.dir <- "~/disk10TB/DARPA/BH_analysis/Shapefiles/Centers/"
# for (i in 1:length(folder_ids)){
#   scenes <- by_folder[[i]]
#   folderid <- scenes$folderid
#   scenes <- split(scenes, scenes$sid)
#   scenepoints2shapefile(scenes,folderid, res.dir)
# }

###################################################################################################################
######################################### Record the video focusing on campus #####################################
###################################################################################################################
################################
# by_folder_campus: list, Each data frame is all scenes info in one folder
# by_scene_campus:  list, Each data frame is one scene info
# scenes_by_folder_campus: list of list, each data frame is one scene info in one folder
tmp <- sceneturretinfo[which(sceneturretinfo$lats.gmb>40.75),]
# tableID_scene <- unique(tmp[,c(1,2,4,5)])
# colnames(tableID_scene) <- 1:nrow(tableID_scene)
# save(tableID_scene, file= "~/disk10TB/DARPA/BH_analysis/SpectralGeoDataTable/tableID_scene.RData")
by_folder_campus <- split(tmp, tmp$folderid)
folder_ids_campus <- unique(tmp$folderid)
by_scene_campus <- split(tmp, tmp$sid)
ext.campus <-  decideExtent(tmp, extent(footprints.spdf[c(unique(tmp$sid), unique(tmp$sid_2)),]))
# extract scenes by folder
extinfo <- vector()
scenes_by_folder_campus <- list()
for (i in 1:length(by_folder_campus)){
  scenes <- by_folder_campus[[i]]
  sids <- c(unique(scenes$sid),unique(scenes$sid_2))
  extinfo <- c(extinfo, decideExtent(scenes, extent(footprints.spdf[sids,])))
  scenes_by_folder_campus[[i]] <- split(scenes, scenes$sid)
}
sum(unlist(lapply(scenes_by_folder_campus, length))); length(by_scene_campus)
length(by_scene_campus)
################################
# plot one folder scenes
# 5, 6, 7, 8, 10, 13, 14
i = 5
plot_scenes(scenes_by_folder_campus[[i]], extinfo[[i]], foot.col="bisque", target.col="black", centerline =T)

################################
# Plot all scenes by folder
x11()
par(mgp = c(2, 0.7, 0))
for (i in 1:length(by_folder_campus)){
  plot_scenes(scenes_by_folder_campus[[i]], extinfo[[i]], foot.col="bisque", target.col="black", sleeptime = 0.4, centerline=T)
}
dev.off()
################################
# Plot all scenes
x11()
par(mgp = c(2, 0.7, 0))
plot_scenes(by_scene_campus, ext.campus, cex.point= 0.03, cex.center=0.4, foot.col="bisque", target.col="black", sleeptime= 0.4)
dev.off()


###################################################################################################################
########################### Colorcoded pixels by range and angle to the plane average #############################
###################################################################################################################
################################
# extract info for one scene
# read "lon", "lat", "dem",  "center.dist", "center.alt" ,"center.range", "center.angle", spectral info into scenes
scene <- by_scene_campus[[3]]
# tmp is a list of two dataframe, since each scene has two sides sid and sid_2
# it takes 21 secs to read one scene-pair using data.frame
# it takes 19 secs to read one scene-pair using data.table
start_time <- Sys.time()
# here the idx is from 1:256 spectrum
tmp <- dem_spec_scene(scene, dir.name, 110)
tmp <- lapply(tmp, dem_spec_scene2center, unlist(scene[round(nrow(scene)/2), c("longs.gmb", "lats.gmb", "alt.gmb")]))
end_time <- Sys.time()
end_time-start_time
View(tmp[[1]])
# here idx is 8+256:  8 geometries
plot_dem_spec_twosides(tmp, 240, col=gray.colors(50))
plot_dem_spec_twosides(tmp, 8, col=rainbow(10, alpha = 0.05), add = T) 
plot_dem_spec_twosides(tmp, 9)

# generate a jpg for one folder at one wavelength 10.010258 um
i = 5; tmps <- list()
values <- vector()
for(j in 1:length(scenes_by_folder_campus[[i]])){
  scene <- scenes_by_folder_campus[[i]][[j]]
  tmp <- dem_spec_scene(scene, dir.name, 110)
  # since the dem_spec_scene2center hasn't been applied, 4 geometries + 1 spectral
  values <- c(values, unname(unlist(tmp[[1]][,5])), unname(unlist(tmp[[2]][,5])))
  tmps[[j]] <- tmp
} 
boxplot(values)
save(tmps, file='~/disk10TB/DARPA/BH_analysis/folder5Spectral.RData')
load('~/disk10TB/DARPA/BH_analysis/folder5Spectral.RData')

minmax<- c(min(values), max(values))  # 546, 2196
i <- 5; minmax <- c(550, 1500)
ext <- extent(footprints.spdf[c(unique(by_folder_campus[[i]]$sid),unique(by_folder_campus[[i]]$sid_2)), ])
ext[2] <- ext[2]+0.003

plottime<-cbind(c(7,13,18,23), c(3,1,1,3), c(1,1.5,1.5,1), c(45,135,225,315))
jpeg("hyperscenes.jpg", width=4000, height=2000, res=300)
par(mfrow=c(1,2), mar=c(3,3,1.7,0), mgp=c(2,1,0), oma=c(0, 0, 0, 0.5))
plot_scenes(scenes_by_folder_campus[[i]], extinfo[[i]], foot.col="bisque", 
            target.col= alpha("black", 0.05), plottime=plottime, centerline =T)
plot_dem_spec_twosides(tmps[[1]], 5,  minmax=minmax, ext=ext, legend=T)
for (tmp in tmps[2:length(tmps)]){
  plot_dem_spec_twosides(tmp, 5, minmax=minmax, add =T)
}
dev.off()




# save wavelength
# all scenes at side 1 or side 2 have the same wavelengths, but different sides have different wavelengths 
wv1 <- c("sid_1", as.numeric(gsub("X([0-9]+\\.[0-9]+).microns", '\\1', colnames(tmp[[1]])[9:264])))
wv2 <- c("sid_2", as.numeric(gsub("X([0-9]+\\.[0-9]+).microns", '\\1', colnames(tmp[[2]])[9:264])))
wavelength_side <- rbind(wv1, wv2)
colnames(wavelength_side) <- c('sid', paste('wv', 1:256, sep=''))
save(wavelength_side, file = '~/disk10TB/DARPA/BH_analysis/SpectralGeoDataTable/wavelength_side.RData')


################################
# read all scenes geo and spectral info into a list, where each element is also a list with 2 datafranes at 2 sides for each scene
for(scene in by_scene_campus){
  gc()   # A call of gc() causes a garbage collection to take place
  tmp <- dem_spec_scene(scene, dir.name)
  chunk <- rbindlist(lapply(tmp, dem_spec_scene2center, unlist(scene[round(nrow(scene)/2), c("longs.gmb", "lats.gmb", "alt.gmb")])), 
                     use.names=F)
  save(chunk, file = paste('~/disk10TB/DARPA/BH_analysis/SpectralGeoDataTable/', paste(scene[1,1], scene[1,2], sep='_'), '.RData', sep=''))
}




