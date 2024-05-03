##############################################################################
# Climate and land use change leading to niche expansion and shifts in birds #
# Pablo M. Avidad, Miguel Clavero, Duarte S. Viana                           #
##############################################################################

# Data preparation

# Load libraries
library(dplyr)
library(tidyr)
library(sp)
library(sf)
library(raster)
library(mapdata)
library(ggplot2)


# Load data
# Bird data and some environmental variables (climate and land use)
bbs <- read.table("BBS_AGGREGATTED.txt", sep="\t", header=T)


#############################################################################
#############################################################################
# DON'T RUN AGAIN
# Add elevation
library(elevatr)
routes <- readOGR("route_buffers", "route_buffers")
proj4string(routes)
routes <- st_as_sf(routes)

(cos(40 * pi/180) * 2 * pi * 6378137) / (256 * 2^10) # resolution for z=10 (~110 m)
elev <- data.frame(site_id=as.integer(routes$RTENO), elev.mean=NA, elev.max=NA)
for(i in 1:nrow(routes)){
  elevi <- get_elev_raster(routes[i,], z = 10)
  elevi.m <- crop(elevi, extent(routes[i,]))
  elevi.crop <- mask(elevi.m, routes[i,])
  elevi.mean <- cellStats(elevi.crop, 'mean')
  elevi.max <- cellStats(elevi.crop, 'max')
  elev[i,c(2,3)] <- c(elevi.mean,elevi.max)
}

save(elev,file="Elevation_routes.RData")

bbs <- left_join(bbs, elev, by="site_id")
write.table(bbs, file="BBS_AGGREGATTED.txt", sep="\t", row.names=FALSE)

#############################################################################
#############################################################################

# Get species
bbs.sp <- data.frame(aou=unique(bbs$aou))
bbs.sp <- left_join(bbs.sp, bbs[!duplicated(bbs$aou),c("aou","species")], by="aou")
bbs.sp$species <- gsub(' ', '_', bbs.sp$species)

#############################################################################
#############################################################################

# Species traits

# Load trait data
traits <- read.table("Data_AVONET.txt",header=T,sep="\t")
traits$species <- gsub(' ', '_', traits$Species1)

# Merge niche and trait data
sp.traits <- dplyr::left_join(bbs.sp,traits,by="species")

# Correct taxonomy
bbs.sp[bbs.sp$species=="Dryobates_villosus","species"] <- "Leuconotopicus_villosus"
bbs.sp[bbs.sp$species=="Dryocopus_pileatus","species"] <- "Hylatomus_pileatus"
bbs.sp[bbs.sp$species=="Icterus_bullockii","species"] <- "Icterus_bullockiorum"
bbs.sp[bbs.sp$species=="Coccothraustes_vespertinus","species"] <- "Hesperiphona_vespertina"
bbs.sp[bbs.sp$species=="Centronyx_henslowii","species"] <- "Passerculus_henslowii"

# Merge niche and trait data again 
sp.traits <- dplyr::left_join(bbs.sp,traits,by="species")

# Delete subspecies
sp.traits<- sp.traits[!is.na(sp.traits$Species1),]

# Choose traits
sp.traits <- sp.traits[,c("aou","species","Species1","Family1","Order1",
                          "Hand.Wing.Index","Mass","Habitat.Density",
                          "Habitat","Migration","Min.Latitude","Max.Latitude")]

#------------------------------

# Calculate range size and proportion

load("Range_polygons.Rdata")

###############
# DO NOT RUN AGAIN
# Add two species that were missing
Cpla <- st_read("Cistothorus_platensis/data_0.shp")
Cpla <- as_Spatial(Cpla)
proj4string(Cpla) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
names(Cpla@data) <- names(ranges@data)
Cpal <- st_read("Cistothorus_palustris/data_0.shp")
Cpal <- as_Spatial(Cpal)
proj4string(Cpal) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
names(Cpal@data) <- names(ranges@data)
ranges <- rbind(ranges,Cpla,Cpal)
#save(ranges,file="Range_polygons.Rdata")
###############

# Check if all species are included
all.sp <- gsub('_', ' ', unique(sp.traits$species))
sp.rg <- as.character(unique(ranges@data$BINOMIAL))
all(all.sp %in% sp.rg)
# Subset and prepare
ranges <- ranges[ranges$BINOMIAL %in% all.sp,]
ranges <- st_as_sf(ranges)
ranges <- st_transform(ranges, crs=3310)

# Make BBS routes total polygon (MCP)
routes <- st_read("route_buffers/route_buffers.shp")
routes <- st_transform(routes, crs=3310)
bbs.points <- st_union(routes)
bbs.pol <- st_convex_hull(bbs.points)

ggplot() + 
  geom_sf(data = bbs.pol) +
  geom_sf(data = routes)

range.prop <- data.frame(species=all.sp, range.size=NA, bbs.size=NA, range.prop=NA, range.cov=NA, range.n=NA)
for(i in all.sp){
  range.sp <- ranges[ranges$BINOMIAL==i,]
  if(length(range.sp)>1 && 2 %in% range.sp$SEASONAL) range.sp <- range.sp[range.sp$SEASONAL==2,]
  if(length(range.sp)>1 && !(2 %in% range.sp$SEASONAL)) range.sp <- range.sp[range.sp$SEASONAL==1,]
  range.sp1 <- range.sp[range.sp$ORIGIN %in% 1:3 & range.sp$PRESENCE==1,]
  try({
    sp.bbs <- st_intersection(bbs.pol, range.sp1)
    range.sp1.area <- sum(st_area(range.sp1))
    sp.bbs.area <- sum(st_area(sp.bbs))
    sp.cov <- st_contains(range.sp1,routes)
    range.prop[range.prop$species==i,2] <- range.sp1.area/1e06
    range.prop[range.prop$species==i,3] <- sp.bbs.area/1e06
    range.prop[range.prop$species==i,4] <- sp.bbs.area/range.sp1.area
    range.prop[range.prop$species==i,5] <- length(sp.cov[[1]])/(sp.bbs.area/1e06)
    range.prop[range.prop$species==i,6] <- length(sp.cov[[1]])
  })
}
# Metadata: https://nc.iucnredlist.org/redlist/content/attachment_files/Legend_combinations_Dec2018.pdf
# ggplot() + 
#   geom_sf(data = range.sp1) +
#   geom_sf(data = bbs.pol,fill=NA) +
#   geom_sf(data = sp.bbs,fill="green")


# Merge niche and trait data again 
range.prop$species <- gsub(' ', '_', range.prop$species)
sp.traits <- left_join(sp.traits, range.prop, by="species")

# Remove species whose sampled range <50%
sp.traits <- sp.traits[sp.traits$range.prop>=0.5,]
sp.traits <- sp.traits[sp.traits$range.n>=100,]

# Number of routes needed to keep coverage constant
sp.traits[sp.traits$range.cov==min(sp.traits$range.cov),]
sp.traits$n.adj <- (105*sp.traits$bbs.size)/3660047
sp.traits$n.adj/sp.traits$bbs.size



#------------------------------


#Specialization Index
load("specialization_index.Rdata")
ssi$species <- gsub(' ', '_', ssi$sci)
sp.traits <- dplyr::left_join(sp.traits,ssi[,c("species","SSI")],by="species")
length(sp.traits$SSI[!is.na(sp.traits$SSI)]) # SSI missing for a few species

# Territoriality
traits2 <- read.table("Traits_Sheard_HWI.txt",header=T,sep="\t")
traits2$species <- gsub(' ', '_', traits2$IUCN.name)
traits2$species[traits2$species=="Troglodytes_troglodytes"] <- "Troglodytes_hiemalis"
sp.traits <- dplyr::left_join(sp.traits,traits2[,c("species","Territoriality")],by="species")

#ACAD 
acad <- read.csv("ACAD Global 2021.02.05-filtered.csv")
#Correct taxonomy
acad[acad$Scientific.Name=="Dryobates villosus","Scientific.Name"] <- "Leuconotopicus villosus"
acad[acad$Scientific.Name=="Dryocopus pileatus","Scientific.Name"] <- "Hylatomus pileatus"
acad[acad$Scientific.Name=="Icterus bullockii","Scientific.Name"] <- "Icterus bullockiorum"
acad[acad$Scientific.Name=="Coccothraustes vespertinus","Scientific.Name"] <- "Hesperiphona vespertina"
acad[acad$Scientific.Name=="Centronyx henslowii","Scientific.Name"] <- "Passerculus henslowii"
sp.traits <- dplyr::left_join(sp.traits, acad[, c("Scientific.Name","PS.g","PT.c", "TB.c","TN.c")],
                                  by = c("Species1" = "Scientific.Name"))



#############################################################################
#############################################################################

# Transform data table to wide format
bbs <- bbs[bbs$aou %in% sp.traits$aou,]
sp.data <- bbs[,c(1:4)]
sp.data <- sp.data %>%
  complete(site_id, time, aou, fill=list(count = 0)) %>%
  spread(aou, count) %>% 
  data.frame()
site_year <- paste(bbs$site_id, bbs$time, sep="")
env <- bbs[!duplicated(site_year),c(1,2,13:ncol(bbs))]
geo <- bbs[!duplicated(site_year),c(1,2,11,12)]

bbs.env <- left_join(sp.data, env, by=c("site_id","time"))
bbs.geo <- left_join(sp.data, geo, by=c("site_id","time"))

# Randomly choose a subset of routes to rarefy range coverage 
nsp <- nrow(sp.traits)
di1 <- bbs.env[bbs.env$time==1981,]
lsrand <- list()
for(s in 1:nsp){
  n.adj <- round(sp.traits$n.adj[sp.traits$aou==as.integer(sub('.', '', names(bbs.env)[s+2]))])
  mr <- matrix(nrow=100,ncol=n.adj)
  for(i in 1:100) mr[i,] <- sample(di1$site_id, n.adj)
  lsrand[[s]] <- mr
}

# Save all pertinent data
save(bbs.env,bbs.geo,sp.traits,lsrand,file="Data_all.RData")


