library(data.table)
library(geosphere)
library(raster)
library(rgdal)
library(rgeos)

source("Z_BH_functions.R")

##### Radiance Unit: flick [W cm^-2 sr^-1 um^-1]
##### BH Radiance Unit: microflicks, 10^-6*[W cm^-2 sr^-1 um^-1], value*10^-6 to flick
##### blackbody_radiance_calculation Unit: [W m^-2 sr^-1 um^-1],  value*10^-4 to flick

######################################################################################################################
######################################################################################################################
######################################################################################################################
########################################## Codes processing the data start here ######################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
crs.ll <- CRS("+proj=longlat")
# load info for datatable in SpectralGeoDataTable folder
load('~/disk10TB/DARPA/BH_analysis/SpectralGeoDataTable/wavelength_side.RData')
load( "~/disk10TB/DARPA/BH_analysis/SpectralGeoDataTable/tableID_scene.RData")

# Read footprints
footprints.fname <- "~/disk10TB/DARPA/BH_analysis/Shapefiles/NittanyRadiancePolygons.shp"
footprints.spdf <- readOGR(footprints.fname) # 1436
rownames(footprints.spdf@data) <- 1:nrow(footprints.spdf@data)
crs(footprints.spdf) <- crs.ll

# find one target on the ground
target <- rbind(c(-77.865033, 40.804466), c(-77.865143, 40.804403), 
                c(-77.864803, 40.804965), c(-77.872292, 40.806454), 
                c(-77.871576, 40.806950), c(-77.870616,  40.807391),
                # old main
                c(-77.861911, 40.795954), c(-77.860965, 40.795642),
                # site 3 and golfing
                c(-77.871460, 40.795407),  c(-77.871750, 40.795862), 
                c( -77.873295, 40.794896), c(-77.873381, 40.793304), 
                c(-77.869845, 40.794539), c(-77.870189, 40.794302)
                )
target.sp  <- SpatialPoints(target, proj4string = crs.ll)
o   <- over(target.sp, footprints.spdf, returnList = T)
sids <- unique(unlist(lapply(o, function(x) as.numeric(rownames(x)))))
sides <- ifelse(sids %in% tableID_scene$sid, 1, 2)

# get the sids to read associated datatables, each datatable contains info of two-side images
over_scenes <- tableID_scene[which(tableID_scene[, 1] %in% sids |  tableID_scene[, 2] %in% sids), c(1,2)]
specfiles <- paste('~/disk10TB/DARPA/BH_analysis/SpectralGeoDataTable/', over_scenes[,1], "_", over_scenes[,2], ".RData", sep = "")

########## Plot one scene and pixels found
# Observation: longitude/latitute have 8-13 decimal digits. Google only returns 6 decimal digits, given any target
# Compare to find nearest pixel, search pixels within the range of given target
# One degree in lat is not equal to one degree in lon in meters, however, we want the search box to be square
# Oberservation: gBuffer is also using lat/lon rather than meters. 
# http://www.csgnetwork.com/degreelenllavcalc.html,  ratio <- 84.389/111.05
# Length in meters of 1° of latitude ~ 111.32 km
# Length in meters of 1° of longitude ~ 40075*cos(latitude)/360 km; 40075*cospi(40.8/180)/360: 84.26827
plot_extractedPixels(target.sp, specfiles, 2, "extractedPixels.jpg", threshold=5e-6,
                     textoffset=c(-0.2,0.5))

########## Extract all pixels from blue heron scenes over the target
start_time <- Sys.time()
result <- pixelsReader(target.sp, specfiles, sides)
end_time <- Sys.time()
end_time-start_time
result.side1 <- result[which(result$side==1),]
result.side2 <- result[which(result$side==2),]
nrow(result.side1); nrow(result.side2)
# check the outliers in extracted pixels, IQR(result1$wv1); sd(result1$wv1)
boxplot(result.side1$wv1)
# remove outliers, 8: center.range; 10: wv1
result.side1 <- remove.outliers(result.side1, usecols= c(8,10)) 
result.side2 <- remove.outliers(result.side2, usecols=c(8,10))
nrow(result.side1); nrow(result.side2)
# write.csv(result.side1, "~/disk10TB/DARPA/BH_analysis/TrainingData/Grass_side1.csv", row.names = F)
# write.csv(result.side2, "~/disk10TB/DARPA/BH_analysis/TrainingData/Grass_side2.csv", row.names = F)












######################################################################################################################
################################## Calculate Downwelling/Upwelling/Transmission ######################################
######################################################################################################################
########## Generate grass spectrum via MODTRAN ######
reflec  <- read.csv("~/disk10TB/DARPA/BH_analysis/MaterialsReflectivity/Grass_reflec.csv")

# calculate blackbody thermal radiance L_T, assuming temperature is 310 Kelvin, 36.85 for grass
# change to microflicks, same as blue heron unit
L_T1 <- blackbody_radiance_calculation(wavelength_side[1, -1], 300)*10^2
L_T2 <- blackbody_radiance_calculation(wavelength_side[2, -1], 300)*10^2
ylim <- c(0,round(max(result.side1[, 10:265], L_T1)+50))
matplot(t(result.side1[, 10:265]), type='l', col="blue", xlab =expression(paste("Wavelength (", mu, "m)")), 
        ylim= ylim, ylab = expression(paste("Microflick 10^-6*(W cm-2 sr-1 ", mu, "m-1)")))
lines(L_T1, type='l', col="red")
legend("bottomright", c("Blackbody Thermal", "Total Radiance"), col=c("red", "blue"), bty='n')

######################################################################################################################
########################################## Derive Algorithm ######################################
######################################################################################################################
########## Final Version of the algorithm
# tau <- exp(-p1/sin(angle)*range)
# upwelling <- p2*range*(1-tau)*tau #p2(λ): blackbody thermal radiance unit factor
# totalRa <- (down*reflec + (1-reflec)*L_T)*tau + upwelling
# totalRa <- (down*reflec + (1-reflec)*L_T + p2*range*(1-exp(-p1/sin(angle)*range)))*exp(-p1/sin(angle)*range)
##################

##################
# side 1 for lambda 1, radiance starts from column 10
deg2rad <- pi/180
angle = result.side1$center.angle*deg2rad; range=result.side1$center.range
# change to microflicks
down.vec <- c(8e-5, 9e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3)*10^6
tau.vec <- seq(0, 1, 0.02)
p1.vec <- -log(tau.vec)*sin(mean(angle))/mean(range)
upwelling.vec <- c(8e-5, 9e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3)*10^6
# p2: length(tau)*length(upwelling) matrix. # has duplicates for each column due to (1-tau)*tau expression
p2.mat <- data.frame(sapply(upwelling.vec, function(x) x/(mean(range)*(1-tau.vec)*tau.vec)))  
nrow(p2.mat); ncol(p2.mat)
##################
# result.side1[, 9+i] <- (down*reflec.side1[i] + (1-reflec.side1[i])*L_T1[i] + p2*range*(1-exp(-p1/sin(angle)*range)))*exp(-p1/sin(angle)*range)
# for ith p1, try all unique p2 at ith row. They are corresponding to the same transmission
total.rad <- result.side1[, 10]
range(total.rad)
reflec <- reflec.side1[1]; LT <- L_T1[1]
possible.mat <- NULL
for (i in 1:length(down.vec)){
  down <- down.vec[i]
  for (j in 1:length(p1.vec)){
    p1 <- p1.vec[j]
    for(p2 in p2.mat[j,]){  
      poss <- (down*reflec + (1-reflec)*LT + p2*range*(1-exp(-p1/sin(angle)*range)))*exp(-p1/sin(angle)*range)
      possible.mat <- cbind(possible.mat, poss)
    }
  }
}

##################
# for a specidifc lambda, downwelling is same 
# self-emitted is same for same pixels
# similar range and angle, the upwelling should be similar
# 8: center.range; 9: center.angle
# t1 = exp(-p1*range1/sin(angle1)).  t2 = exp(-p1*range2/sin(angle2))
# t1/t2 = exp(-p1*range1/sin(angle1) - (-p1*range2/sin(angle2))) = exp(p1[range2/sin(angle2) -range1/sin(angle1)])
# log(t1/t2) = p1[range2/sin(angle2)-range1/sin(angle1)]
caltransPara <- function(obs, row1, row2, wv.col){
  exponent <- log(obs[row1, ..wv.col]/obs[row2, ..wv.col])
  if(exponent == 0){
    print("observations are same, cannot calculate p")
    return()
  }
  r.over.a <- obs[row2, 8]/sin(obs[row2, 9]) - obs[row1, 8]/sin(obs[row1, 9])
  p1 <- unname(unlist(exponent/r.over.a))
  return(p1)
}
caltransPara(result.side1,1,2,10)
caltransPara(result.side1,3,4,10)
caltransPara(result.side1,4,5,10)
caltransPara(result.side1,5,6,10)
caltransPara(result.side1,6,7,10)
caltransPara(result.side1,7,8,10)

# distribution of all p1,
# first, group by similar observations
ncluster <- 12
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
l = range(result.side1[,8])[1]; r = range(result.side1[,8])[2]
centers <- seq(range(result.side1[,8])[1], range(result.side1[,8])[2], (r-l)/(ncluster+1))[-c(1,ncluster+2)]
fit <- kmeans(result.side1[,8], centers = centers)
plot(result.side1[,c(8,9)], col=colfunc(ncluster)[fit$cluster])
# for each cluster with similar observations
p1.values <- c()
for(i in 1:ncluster){
  idx <- which(fit$cluster == i)
  # comb <- t(comb)
  for (i in 2:length(idx)){
    p1.values <- c(p1.values, caltransPara(result.side1,i-1,i,10))
  }
}
boxplot(p1.values)
p1.values.filter <- p1.values[which(p1.values < 0.15 & p1.values > -0.15)]
boxplot(p1.values.filter)



######################################################################################################################
######################################################################################################################
######################################################################################################################
########## Refer to DARPA Notes 3.7 Atmoshpheric Propagation
# total radiance = (down*reflec + (1-reflec)*L_T)*tau + upwelling 
# downwelling: same for all pixels
# L_T : blackbody self-emitted radiance, depending on T
# tau: transmission between target and sensor, depending on p1, range, angle
# LOS upwelling : upwelling solar and thermal radiance, depending on p2, range, angle
##########
# Lambert’s/Bouguer’s law for a homogeneous medium: τ=e^(−β_a*z); optical depth δa = β_a*z
# β(λ,θv)=p1(λ)*csc(angle)=p1(λ)/sin(angle): fractional amount of flux absorbed per unit length
# p1(λ): the absorption factor along the vertical atmosphere
# The total optical depth δ(λ) =β_a*z decreases as the elevation angle increases.
##########
# Assumption: tau increases when elevation angle increase; angle is from 30 to 90 so sin(angle) is from 0.5 to 1
# tau: exp(1)^(-p1/sin(angle)*z): tau increases while p1/sin(angle)*z decreases when elevation angle increases
# Assumption: path thermal radiance decreases when elevation angle increase. This is cause by two factors
# tau decreases but thermal emission increases when angle decreases, however, the thermal emission effect is dominated 
# Luελ(θ,φ) = sum(L_Ti(λ)[τ_i(λ) − τ_i+1(λ)]) = sum(ith layers' (thermal radiance *  emissivity * transmission to the sensor)
tmp <-blackbody_radiance_calculation(wavelength_side[1,2:ncol(wavelength_side)], 320)
solar_bbradiance <- blackbody_radiance_calculation(wavelength_side[1, 2:ncol(wavelength_side)], temp = 5778)
solar_bbradiance <- as.data.frame(cbind(wv, radiance=solar_bbradiance*10^-4)) # [W cm^-2 sr^-1 um^-1]
