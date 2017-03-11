#Clear workspace
rm(list=ls())

###################TABLE OF CONTENTS

###[1] Creating grid that will be used for connectivity
###[2] Identifying release locations and their place in grid
###[3] Identifying settlement locations and linking to release locations


#load in required packages
require(data.table)
require(tidyverse)
require(readr)
require(rgdal)
require(rgeos)
require(maptools)
#setwd("C:/Christopher_MSc/Github/ParticleTracking")
getwd()
#testing adding code from home computer

############################################################
############################################################
### [1] Creating connectivity grid

#Loading a grid that you will divide BC ocean into
grid <- readOGR("./cuke_present/StudyExtent", "grid")
#Dissolve into one polygon since so you can change grid dimensions
My_BC_projection <- CRS("+proj=aea +lat_1=47 +lat_2=54 +lat_0=40 +lon_0=-130 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")
ConPoly <- spTransform(grid, My_BC_projection) #use your custom BC projection for this
ConPoly <- gUnaryUnion(ConPoly)
plot(ConPoly)
rm(grid)

#Adding dataframe to convert back to SpatialPolygonsDataFrame
ConPoly_ID <- row.names(ConPoly)
ConPoly_ID <- as.data.frame(ConPoly_ID)
# And add the data back in
ConPoly <- SpatialPolygonsDataFrame(ConPoly, ConPoly_ID)
rm(ConPoly_ID)

#Make your Connectivity grid into a raster, then convert back into a shapefile so you can see how many larvae move from release site to final location
library(raster)
r <- raster(extent(ConPoly))
projection(r) <- proj4string(ConPoly)
res(r) <- 3000

ptm <- proc.time()
ConGrid <- rasterize(ConPoly, r) #5 seconds
proc.time() - ptm
plot(ConGrid)
rm(r)

#Back into polygon
ConPoly <- rasterToPolygons(ConGrid, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) #couple of seconds, can be much longer depending on file size
head(ConPoly@data)
#Adding in unique IDs for each polygon
ConPoly@data$Poly_ID <- as.numeric(row.names(ConPoly))

writeOGR(ConPoly, dsn = "./output", "ConPoly",
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)



########################################################################
########################################################################
### [2] Loading up larval release points

#Acquiring files
filenames <- list.files(path = "./cuke_present/ReleaseLocations", pattern="rl_.", full.names=TRUE,recursive=T)

# load all files into a list, read_csv is much faster than read.csv
rllist <- lapply(filenames, read_csv,
                 col_names = c("long0","lat0","Z0","delay","site0"),
                 col_types = cols("d","d","i","i","i")
                 )

# set the names of the items in the list, so that you know which file it came from
rllist <- setNames(rllist,filenames)

# rbind the list
rl <- rbindlist(rllist, idcol="filename")

rl$bin <- gsub(".*rl_|.txt.*", "",rl$filename)
head(rl)

#Creating csv file ith all starting locations
write.csv(rl, file="./output/release_locations.csv", row.names = F)


### [2b] Creating map of release locations
release_points <- subset(rl, select = c(long0, lat0, Z0))
release_points <- as.data.frame(release_points)

xy <- release_points[,c(1,2)]
NAD_projection <- proj4string(grid)
My_BC_projection <- CRS("+proj=aea +lat_1=47 +lat_2=54 +lat_0=40 +lon_0=-130 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")

release_larvae <- SpatialPointsDataFrame(coords = xy, data = release_points, proj4string = CRS(NAD_projection))
release_larvae <- spTransform(release_larvae, My_BC_projection) #use your custom BC projection for this
plot(release_larvae)

#Associate points with where they were released from
require(spatialEco)
release_larvae_map <- point.in.poly(release_larvae, ConPoly)
head(release_larvae_map@data)
plot(release_larvae_map)

#write out points
writeOGR(release_larvae_map, dsn = "./output", layer = "release_points", driver = "ESRI Shapefile", overwrite = TRUE)



########################################################################
########################################################################
#[3] Identifying settlement locations and linking to release locations

# List the particle tracking files for that particular year and pld
year <- 1999
pld <- 120

#Acquiring files
filenames <- list.files(path=paste0("./cuke_present/ConData/G",year), pattern=glob2rx(paste0("*para    1",formatC(pld+1, width = 3, format = "d", flag = "0"),"*")), full.names=TRUE,recursive=T)

# load all files into a list, read_csv is much faster than read.csv
datalist <- lapply(filenames, read_csv,
                   col_names = c("long","lat","Z","Out","site"),
                   col_types = cols("d","d","d","i","i")
                   )

# set the names of the items in the list, so that you know which file it came from
datalist <- setNames(datalist,filenames)

# rbind the list
dataset <- rbindlist(datalist, idcol="filename")
dataset$site <- NA
rm(datalist)

#Reshaping dataset to take filename info and turning it into columns
dataset <- dataset %>%
  mutate(temp=substr(filename,24,nchar(filename))) %>%
  separate(temp,c("temp_type_year","rday","bin","time"),"/",convert=TRUE) %>% 
  separate(temp_type_year,c("type","year"),sep=1,convert=TRUE) %>% 
  mutate(time=as.integer(substr(time,12,15)))

#Linking release locations to settlement locations based on bin
for(i in unique(dataset$bin)){
  x <- rl$bin==i
  y <- dataset$bin==i
  dataset$long0[y] <- rl$long0[x]
  dataset$lat0[y] <- rl$lat0[x]
  dataset$Z0[y] <- rl$Z0[x]
  dataset$delay[y] <- rl$delay[x]
  dataset$site0[y] <- rl$site0[x]
}
head(dataset)
View(dataset)

#Giving each larvae a unique ID
Con_df <- dataset
Con_df <- subset(Con_df, select = -c(filename, Out, site, time))
Con_df$larvae_ID <- row.names(Con_df) 

Release_df <- subset(Con_df, select = c(long0, lat0, Z0, larvae_ID))
Settle_df <- subset(Con_df, select = c(long, lat, Z, larvae_ID))





























#####################extra TO BE REINSERTED




###Clipping release map locations to intertidal, nearshore, and offshore habitat
head(release_larvae@data)
#Intertidal
Intertidal_release <- release_larvae[release_larvae@data$Z0 > -10,]
Intertidal_release <- Intertidal_release[Intertidal_release@data$Z0 < 5,]

writeOGR(Intertidal_release, dsn = "C:/Christopher_MSc/Remi_data/ParticleTracking/present/intertidal", layer = "Intertidal_release", driver = "ESRI Shapefile")

#Nearshore
Nearshore_release <- release_larvae[release_larvae@data$Z0 > -60,]
Nearshore_release <- Nearshore_release[Nearshore_release@data$Z0 < -10,]

writeOGR(Nearshore_release, dsn = "C:/Christopher_MSc/Remi_data/ParticleTracking/present/nearshore", layer = "Nearshore_release", driver = "ESRI Shapefile")

#Offshore
Offshore_release <- release_larvae[release_larvae@data$Z0 > -250,]
Offshore_release <- Offshore_release[Offshore_release@data$Z0 < -60,]

writeOGR(Offshore_release, dsn = "C:/Christopher_MSc/Remi_data/ParticleTracking/present/offshore", layer = "Offshore_release", driver = "ESRI Shapefile")

