#Clear workspace
rm(list=ls())

###################TABLE OF CONTENTS

###[1] Loading up larval release points
###[2] Identifying settlement locations and linking to release locations
###[3] Setting up study extent you will be using to clip your larval release points to your BC study extent
###[4] Creating connectivity matrices from Con_df dataframe
###

#load in required packages
require(data.table)
require(tidyverse)
require(readr)
require(rgdal)
require(rgeos)
require(maptools)
#setwd("C:/Christopher_MSc/Github/ParticleTracking")
getwd()


########################################################################
########################################################################
### [1] Loading up larval release points

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

rl$bin <- as.numeric(gsub(".*rl_|.txt.*", "",rl$filename))
head(rl)
rm(rllist)

#Creating csv file ith all starting locations
write.csv(rl, file="./output/release_settlement/Release_lat_long.csv", row.names = F)


### [1b] Creating map of release locations
#release_points <- subset(rl, select = c(long0, lat0, Z0))
#release_points <- as.data.frame(release_points)

#xy <- release_points[,c(1,2)]
#grid <- readOGR("./cuke_present/StudyExtent/Starting_grid", "grid") #Using Remi's grid to get NAD projection
#NAD_projection <- proj4string(grid)
#My_BC_projection <- CRS("+proj=aea +lat_1=47 +lat_2=54 +lat_0=40 +lon_0=-130 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")


#release_larvae <- SpatialPointsDataFrame(coords = xy, data = release_points, proj4string = CRS(NAD_projection))
#release_larvae <- spTransform(release_larvae, My_BC_projection) #use your custom BC projection for this
#rm(grid, xy)
#plot(release_larvae)

#write out points
#writeOGR(release_larvae, dsn = "./output/shapefiles", layer = "Released_locations", driver = "ESRI Shapefile", overwrite = TRUE)


########################################################################
########################################################################
#[2] Identifying settlement locations and linking to release locations

# List the particle tracking files for that particular year and pld
year <- 1999
pld <- 120

#Acquiring files
filenames <- list.files(path=paste0("./cuke_present/ConData/G",year), pattern=glob2rx(paste0("*para    1",formatC(pld+1, width = 3, format = "d", flag = "0"),"*")), full.names=TRUE,recursive=T)
# filenames <- list.files(path=paste0("F:/MPA_particles/output/G",year), pattern=glob2rx(paste0("*para    1",formatC(pld+1, width = 3, format = "d", flag = "0"),"*")), full.names=TRUE,recursive=T)

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


###This process takes a long time ~ 5 - 10 minutes

#Reshaping dataset to take filename info and turning it into columns
dataset <- dataset %>%
  mutate(temp=substr(filename,24,nchar(filename))) %>%
  # mutate(temp=substr(filename,25,nchar(filename))) %>% # you probably want this back to 24?
  separate(temp,c("temp_type_year","rday","bin","time"),"/",convert=TRUE) %>% 
  separate(temp_type_year,c("type","year"),sep=1,convert=TRUE) %>% 
  mutate(time=as.integer(substr(time,9,13))-1001)


#Testing combining dataset with release locations
dataset2 <- dataset
rl2 <- rl

dataset2$long0[dataset2$bin==1] <- rl2$long0[rl2$bin==1]
dataset2$lat0[dataset2$bin==1] <- rl2$lat0[rl2$bin==1]
dataset2$Z0[dataset2$bin==1] <- rl2$Z0[rl2$bin==1]
dataset2$delay[dataset2$bin==1] <- rl2$delay[rl2$bin==1]
dataset2$site0[dataset2$bin==1] <- rl2$site0[rl2$bin==1]

write.csv(dataset2, "./output/connectivity_tables/testing_dataset.csv", row.names = F)

#Linking release locations to settlement locations based on bin
for(i in unique(dataset$bin)){
  x <- rl$bin==i
  y <- dataset$bin==i
  dataset$long0[y] <- rl$long0[x]
  dataset$lat0[y] <- rl$lat0[x]
  dataset$Z0[y] <- rl$Z0[x]
  dataset$delay[y] <- rl$delay[x]
  dataset$site0[y] <- rl$site0[x]
  print(paste(i,sum(x),sum(y),sum(is.na(dataset$long0)))) # this is just to show its working
}
head(dataset)
rm(filenames,x,y,i)

#Write out dataset - takes a long time
write.csv(dataset, "./output/connectivity_tables/dataset.csv", row.names = F) #6 minutes



#END OF CODE RELEVANT TO REMI
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################






#Add larvae IDs to dataset
Con_df <- dataset
Con_df <- subset(Con_df, select = c(long0, lat0, Z0, long, lat, Z, year, rday))
Con_df$larvae_ID <- row.names(Con_df)

#Write out Con_df - takes a long time
write.csv(Con_df, "./output/connectivity_tables/Con_df.csv", row.names = F)
#rm(dataset)

########################################################################
########################################################################
#[3] Setting up study extent you will be using to clip your larval release points to your BC study extent

#Clipping to your study extent
Ecozone_mask <- readOGR("./cuke_present/StudyExtent/BC_EcozonesMask", "Ecozone_mask")

#Loading Remi's grid where larvae were released
grid <- readOGR("./cuke_present/StudyExtent/Starting_grid", "grid")
#Dissolve into one polygon since so you can change grid dimensions
grid <- spTransform(grid, My_BC_projection) #I think this works fine, if not use: Ecozone_mask@proj4string
grid <- gUnaryUnion(grid)

#Intersecting
ConPoly <- gIntersection(grid, Ecozone_mask, byid = FALSE, drop_lower_td = TRUE) #This works, but you'll have to choose a shapefile that includes islands and doesn't cut-off at rivers 

#Adding dataframe so you can create a shapefile of new study extent
ConPoly_ID <- 1 
ConPoly <- SpatialPolygonsDataFrame(ConPoly, as.data.frame(ConPoly_ID))
ConPoly@data
plot(ConPoly)
rm(grid, Ecozone_mask, ConPoly_ID)


#Turn your study extent into raster to analyse movement between cells
require(raster)
r <- raster(extent(ConPoly))
projection(r) <- proj4string(ConPoly)
res(r) <- 3000

ConGrid <- rasterize(ConPoly, r) #5 seconds
rm(r)

#Convert back to polygon for useable shapefile
ConPoly <- rasterToPolygons(ConGrid, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) #couple of seconds, can be much longer depending on file size
rm(ConGrid)
#Adding in unique IDs for each polygon
ConPoly@data$Poly_ID <- as.numeric(row.names(ConPoly))
#And removing useless "layer" column from attribute table
ConPoly <- ConPoly[,!(names(ConPoly) %in% "layer")]
head(ConPoly@data)
plot(ConPoly)

########################################################################
########################################################################
#[4] Creating connectivity matrices from Con_df dataframe

#Showing where each larvae begings and ends
Release_df <- subset(Con_df, select = c(long0, lat0, Z0, larvae_ID))
Settle_df <- subset(Con_df, select = c(long, lat, Z, larvae_ID, year, rday))

#Associate points with where they were released from
xy <- subset(Release_df, select = c(long0, lat0))

write.csv(xy, "C:/Christopher_MSc/temp/xy.csv", row.names = F)
Released_larvae <- SpatialPointsDataFrame(coords = xy, data = Release_df, proj4string = CRS(NAD_projection))
Released_larvae <- spTransform(Released_larvae, My_BC_projection) #use your custom BC projection for this

#write out points
writeOGR(release_larvae, dsn = "./output/shapefiles", layer = "release_points", driver = "ESRI Shapefile", overwrite = TRUE)

require(spatialEco)
Release_df <- point.in.poly(release_larvae, ConPoly)
#head(release_larvae_map@data)



##########
#########
########
#######
######
#####
####
###
##
#END

















#####################extra TO BE REINSERTED




############################################################
############################################################
### [1] Creating connectivity grid



writeOGR(ConPoly, dsn = "./output", "ConPoly",
         driver = "ESRI Shapefile", verbose = TRUE, overwrite = TRUE, morphToESRI = TRUE)




################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################






My_BC_projection <- CRS("+proj=aea +lat_1=47 +lat_2=54 +lat_0=40 +lon_0=-130 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")







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

