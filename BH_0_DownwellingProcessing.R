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

library(gtools)
source("Z_MODTRAN_functions.R")
source("Z_global_variables.R")

# Downwelling .dwr files are processed via the matlab codes
down.folder <- '~/disk10TB/DARPA/BH_analysis/NittanyRadiance_Downwelling/D_P'
sam.files <- list.files(down.folder, pattern='.sam', full.names = T)
cbb.files <- list.files(down.folder, pattern='.cbb', full.names = T)
wbb.files <- list.files(down.folder, pattern='.wbb', full.names = T)
dwr.files <- list.files(down.folder, pattern='.dwr', full.names = T)
length(sam.files); length(cbb.files); length(wbb.files); length(dwr.files)

# list files of txt data generated from dwr files
downtxt.files <- list.files('~/disk10TB/DARPA/BH_analysis/NittanyRadiance_Downwelling/DWR_txt', full.names = T)

# Date: 4; spectrum: 49; wavelength:50
read.downtxt <- function(file){
  data <- scan(file,what="Character",sep="\n")
  # time
  re <-  regexpr(":?[0-9]{2}:[0-9]{2}:[0-9]{2}",data[4])
  time <- regmatches(data[4], re)
  # spectrum
  spectrum <- sub(".*: ", "", data[49])
  spectrum <- as.numeric(strsplit(spectrum, ' ')[[1]][-1])
  # wavelength
  wavelength <- sub(".*: ", "", data[50])
  wavelength <- as.numeric(strsplit(wavelength, ' ')[[1]])
  info <- as.data.frame(cbind(wavelength, spectrum))
  return(list(time, info))
}

file <- downtxt.files[4]
temp <- read.downtxt(file)
time <- temp[[1]]
info <- temp[[2]]

plot(info, type='l', xlab=expression(paste("Wavelength (", mu, "m)")), 
     ylab=expression(paste("Downwelling Radiance (W cm-2 sr-1 ", mu, "m-1)")))

# here guess the wavelength is wavenumber
info <- info[which(info$wavelength>=400),]
plot(info, type='l', xlab=expression(paste("Wavelength (", mu, "m)")), 
     ylab=expression(paste("Downwelling Radiance (W cm-2 sr-1 ", mu, "m-1)")))
info$wavelength = 10000/info$wavelength
plot(info, type='l', xlab=expression(paste("Wavelength (", mu, "m)")), 
     ylab=expression(paste("Downwelling Radiance (W cm-2 sr-1 ", mu, "m-1)")))


temp.scan <-  mixedsort(list.files('~/disk10TB/DARPA/MODTRANSimulated/Stage1_7.5_12um/Polyethylene',pattern="*scan.csv",full.names=TRUE))
simulated <- MODTRAN.readfiles.scan(temp.scan)
table.scan <-  MODTRAN.combineGeo(simulated, eles)
table.scan<- MODTRAN.calculate.downwelling(table.scan)
table.scan <- table.scan[which(table.scan$theta==90),]
lines(table.scan$wavelength, table.scan$downwelling)
lines(table.scan$wavelength, 5*pi*table.scan$downwelling, col='red')
