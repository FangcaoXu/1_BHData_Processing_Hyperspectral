# simulate different range and elevation angle
# MODTRAN takes the zenith as input to its json file.

theta <- seq(30, 60, 1)
range <- seq(3000, 6500, 100)
deg2rad <- pi/180

length(theta); length(range)
sensor <- NULL
for(a in theta){
  for(r in range){
    ele <- r*sin(a*deg2rad)
    sensor <- rbind(sensor, c(90-a, ele, r)) # save zenith rather than elevation angle
  }
}

colnames(sensor) <- c('theta', 'elevation', 'range')
  
write.csv(sensor, '/amethyst/s0/fbx5002/PythonWorkingDir/sensorLocation/sensorAngleEle.csv')

