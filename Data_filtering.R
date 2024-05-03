##############################################################################
# Climate and land use change leading to niche expansion and shifts in birds #
# Pablo M. Avidad, Miguel Clavero, Duarte S. Viana                           #
##############################################################################


# Get BBS data and environmental data
# Some code from Harris et al. 2018 (https://doi.org/10.7717/peerj.4278)

# Load functions
source("Data_prep_functions.R")

# Load libraries
library(tidyverse)


#############################################################################################

# Prepare BBS data
path <- "~/Dropbox (Personal)/Paper_sCom_birds/Data/2020Release_Nor"
setwd(path)
knitr::opts_knit$set(root.dir = path)

# Unzip the data
unzipBBS(path)

# Read in the data
bbs <- readBBS(path = path)

# Summarize the data
bbs <- summarizeBBS(bbs)
bbs$aou <- as.integer(bbs$aou)
bbs$count <- as.integer(bbs$count)

# Change AOU codes to species names
bbs <- aouToSp(bbs, path = path)

# filter species
#' Removes waterbirds, shorebirds, owls, kingfishers, knightjars,
#' dippers. These species are poorly sampled due to their aquatic or
#' noctural nature. Also removes taxa that were either partially unidentified
#' (e.g. "sp.") or were considered hybrids (e.g. "A x B") or were listed as more
#' than one species (e.g. "A / B")
valid_taxa <-  bbs %>%
  filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010)
bbs <- bbs[bbs$aou %in% valid_taxa$aou,]

# Get route data
locs <- read.csv(paste(path, "routes.csv", sep="/"), colClasses = "character") %>%
  mutate(id = paste0(StateNum, Route)) %>%
  dplyr::select(id, Latitude, Longitude)
bbs <- left_join(bbs, locs, by = "id")

# Filter data
weather <- read.csv(paste(path, "weather.csv", sep="/"), colClasses = "character") 
weather <-  weather %>% 
  mutate(id = paste0(StateNum, Route), year = Year) %>%
  dplyr::select(id, year, ObsN, RunType, RPID)
weather <- weather[weather$RunType == 1 & weather$RPID ==101,]
bbs <- inner_join(bbs, weather, by = c("id", "year"))

bbs2 <- combine_subspecies(bbs)
bbs <- as.data.frame(bbs2)

bbs <- bbs[!is.na(bbs$aou),]
bbs$year <- as.integer(bbs$year)
bbs$Latitude <- as.numeric(bbs$Latitude)
bbs$Longitude <- as.numeric(bbs$Longitude)

setwd("~/Dropbox (Personal)/Paper_sCom_birds/Data")
save(bbs, file="BBS_data.RData")
load("BBS_data.RData")


#----------------------------------------------------

# Data filtering

# Only routes below lat 55
bbs <- bbs[bbs$Latitude<55,]

# Common routes in 3-year periods for the entire time span (1980-2019)
ss <- list()
t <- seq(1980,2016,3)
for(i in 1:length(t)) ss[[i]] <- unique(bbs$id[bbs$year %in% t[i]:(t[i]+2)])
common <- Reduce(intersect, ss)
bbs2 <- bbs[bbs$id %in% common,]


#----------------------------------------------------

# Visualisation and exploration
library(sp)
library(raster)

# All sites surveyed across all years
sites <- bbs2[!duplicated(bbs2$id),c("Longitude","Latitude","id")]
row.names(sites) <- sites$id
sites.sp<-SpatialPointsDataFrame(sites[,1:2],data=data.frame(sites$id),
                                 proj4string=CRS("+proj=longlat +datum=WGS84"))

# Map sites
library("mapdata")
map('worldHires')
points(sites.sp, col='red')
(e<-bbox(extent(sites.sp)*1.2))

quartz(height=4,width=6)
par(mar=c(0,0,0,0))
map('worldHires', xlim = e[1, ], ylim = e[2, ])
map.axes()	
points(sites.sp, col='darkgreen', pch=16, cex=0.6)



#############################################################################################

# Add climate data

library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(daymetr)


bbs2$id <- as.character(as.integer(bbs2$id))

# load routes shapefile
routes <- readOGR("nabbs02_mis_alb", "nabbs02_mis_alb")
proj4string(routes) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
routes <- routes[routes$RTENO %in% unique(bbs2$id),]

dup.routes <- routes$RTENO[duplicated(routes$RTENO)]
routes.out <- c()
for(i in unique(dup.routes)){
  routesi <- routes@data[routes$RTENO==i,]
  routes.out <- c(routes.out, routesi[-which.max(routesi$ARC_LENGTH),"SEQNO"])
}
routes <- routes[!(routes$SEQNO %in% routes.out),]

# Export as shapefile
setwd("~/Dropbox (Personal)/Paper_sCom_birds/Data/BBS_routes")
writeOGR(routes, ".", layer='BBS_routes', driver="ESRI Shapefile") 

setwd("~/Dropbox (Personal)/Paper_sCom_birds/Data")
#save(routes,file="routes_shapefile.RData")
load("routes_shapefile.RData")
load("clim.RData")


# Extract climate for each route
buff <- list()
for(i in 1:nrow(routes)){
  routei <- routes[i,]
  centri <- SpatialPoints(coordinates(as(extent(routei), "SpatialPolygons")))
  proj4string(centri) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  buffi <- gBuffer(centri, width=20000, id=routes$RTENO[i])
  buff[[i]] <- spTransform(buffi, CRS("+proj=longlat +datum=WGS84"))
}

buff.dt <- routes@data
row.names(buff.dt) <- routes$RTENO

# Make a single satialPolygonsDataFrame object
list.df <- lapply(seq_along(buff), function(i){buff[i][[1]]})
buff.sp <- do.call("rbind", c(args = list.df, makeUniqueIDs = TRUE))
buff.sp <- SpatialPolygonsDataFrame(buff.sp, buff.dt[,c("RTENO","SEQNO","SRTENAME")], match.ID = TRUE)
# Export as shapefile
setwd("~/Dropbox (Personal)/Paper_sCom_birds/Data/route_buffers")
writeOGR(buff.sp, ".", layer='route_buffers', driver="ESRI Shapefile") 



nc <- length(1980:2019)*12*nrow(routes)*5*3
#clim <- as.matrix(data.frame(year=integer(nc),month=integer(nc),value=numeric(nc),clim_var=character(nc),site_id=character(nc)))
clim <- matrix(NA,nrow=0,ncol=5)
# RUN ON REMOTE PLATFORM OR CLUSTER
for(i in 773:nrow(routes)){ 
  buffi <- buff[[i]]
  bbi <- bbox(buffi)
  location <- c(bbi[2,2],bbi[1,1],bbi[2,1],bbi[1,2])
  download_daymet_ncss(location = location, start = 1980, end = 2019, frequency = "monthly", 
                       param = c("tmin","tmax","prcp"), path = tempdir(), silent = TRUE)
  all.r <- list.files(tempdir(), pattern = ".nc")
  for(j in 1:length(all.r)){
    r <- raster::stack(file.path(tempdir(),all.r[j]))
    projection(r) <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
    r <- raster::projectRaster(r, crs = "+proj=longlat +datum=WGS84")
    r2 <- crop(r, buffi)
    r3 <- mask(r2, buffi)
    mean.r3 <- data.frame(value=cellStats(r3, "mean"), site_id=routes@data$RTENO[i])
    mean.r3$year <- as.integer(substr(row.names(mean.r3),2,5))
    mean.r3$month <- as.integer(substr(row.names(mean.r3),7,8))
    mean.r3$clim_var <- substr(all.r[j],1,4)
    mean.r3$site_id <- routes@data$RTENO[i]
    mean.r3 <- as.matrix(mean.r3[,c('year','month','value','clim_var','site_id')])
    clim <- rbind(clim,mean.r3)
  }
}

clim <- as.data.frame(clim)
#  put zeros left of  in site_id (due conversion of character to integer when using as.matrix...)
sb <- unique(bbs$site_id)
sc <- unique(clim$site_id)
not_complete <- sc[which(!(sc %in% sb))]
clim$site_id[clim$site_id %in% not_complete] <- 
  paste("0",clim$site_id[clim$site_id %in% not_complete],sep="")
sc[which(!(unique(clim$site_id) %in% sb))]

clim$year <- as.integer(clim$year)
clim$month <- as.integer(clim$month)
clim$value <- as.numeric(clim$value)
clim$clim_var[clim$clim_var=="prcp"] <- "ppt"

save(buff,clim,file="clim.RData")

length(unique(clim[,5]))
table(clim[,5])


# Calculate BIOCLIM variables

clim2 <- spread(clim, key=clim_var, value=value)
clim2$tmean <- apply(clim2[,4:6],1,mean)
clim3 <- gather(clim2, key=clim_var, value=value, ppt, tmax, tmin, tmean)
clim3 <- clim3[,c(1,2,5,4,3)]

bioclim <- process_bioclim_data(clim3)
bioclim <- as.data.frame(bioclim)

# Match BBS and climate data
names(bbs)[1] <- "site_id"
bbs.clim <- inner_join(bbs, bioclim, by = c('site_id', 'year')) %>% data.frame()
write.table(bbs.clim, file="BBS_BIOCLIM.txt", sep="\t", row.names=FALSE)




#############################################################################################

# Filter species and aggregate counts in 3-year periods

# read land cover data
lc <- read.table("BBS_LC_1980_2019.txt", header=T, sep="\t")


# read BBS data
setwd("~/Dropbox (Personal)/Paper_sCom_birds/Data")
bbs.clim <- read.table("BBS_BIOCLIM.txt", header=T, sep="\t")

# Merge land cover data
bbs <- left_join(bbs.clim, lc, by = c('site_id', 'year')) %>% data.frame()

setwd("~/Dropbox (Personal)/Paper_sCom_birds/Data/BBS_data_final")
write.table(bbs, file="BBS_ALL_DATA.txt", sep="\t", row.names=FALSE)

bbs <- bbs[bbs$year!=2019,]
t <- seq(1980,2019,3)
bbs$time <- cut(bbs$year, t, labels=seq(1981,2017,3), include.lowest = T, right=F)
# Summarise data by averaging in 3-year periods and across observers
bbs.agg <- bbs %>% 
  group_by(site_id, time, aou) %>%
  summarise(count=ceiling(mean(count)),across(names(bbs)[c(4:9,11,12)],unique),across(names(bbs)[16:48],mean)) %>% 
  data.frame()

# Checking...
tapply(bbs.agg$site_id,bbs.agg$time,function(x) length(unique(x)))


# choose species with a given occupancy
total.routes <- length(unique(bbs.agg$site_id))
sps <- unique(bbs.agg$aou)
sp.occ <- data.frame(aou=sps, occ=NA)
for(i in 1:length(sps)){
  spi <- sps[i]
  bbsi <- bbs.agg[bbs.agg$aou==spi,]
  sp.occ[i,"occ"] <- min(tapply(bbsi$site_id, bbsi$time, function(x) length(x)))
}

# number of species occurring in at least 20 routes per time period across years
sp.occ2 <- na.exclude(sp.occ)
nrow(sp.occ2[sp.occ2$occ>=20,])
sp.in <- sp.occ2$aou[sp.occ2$occ>=20]
# Filter species
bbs.agg <- bbs.agg[bbs.agg$aou %in% sp.in,]

setwd("~/Dropbox (Personal)/Paper_sCom_birds/Data/BBS_data_final")
write.table(bbs.agg, file="BBS_AGGREGATTED.txt", sep="\t", row.names=FALSE)






