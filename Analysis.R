##############################################################################
# Climate and land use change leading to niche expansion and shifts in birds #
# Pablo M. Avidad, Miguel Clavero, Duarte S. Viana                           #
##############################################################################

# Analysis


# Load libraries
library(plyr)
library(dplyr)
library(tidyr)
library(sp)
library(sf)
library(raster)
library(mapdata)
library(ggplot2)
library(doParallel)


# Load data
load("Data_all.RData")

#############################################################################
#############################################################################

# Analysis
library(ecospat)
library(ade4)

# Choose final set of environmental variables for niche characterisation
# Remove very collinear variables
nsp <- nrow(sp.traits)
nper <- 13
cormat <- cor(bbs.env[,c((nsp+3):ncol(bbs.env))])
bbs.env <- subset(bbs.env, select=-c(bio3,bio5,bio6,bio8,bio9,bio13,bio14,bio18,bio19,
                                     secdn,secmb,c3nfx,elev.max))

# global environmental PCA 
period <- unique(bbs.env$time)
env.global <- bbs.env[,(nsp+3):ncol(bbs.env)]
pca.env <- dudi.pca(env.global,scannf=F,nf=2)
summary(pca.env)
# Plot Variables Contribution with ecospat.plot.contrib()
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
# PCA scores for the whole study area across periods
scores.globclim <- pca.env$li


#############################################################################
#############################################################################

# Niche change components
# fixing either the geographical distribution or the environment 
# at the reference period (1981)

# Initial period
di1 <- bbs.env[bbs.env$time==1981,]

cl <- makeCluster(8)
registerDoParallel(cl)
clusterExport(cl, list("nsp","period","nper","bbs.env","di1","lsrand","env.global","pca.env","scores.globclim"))
clusterCall(cl,function(){library(ecospat);library(ade4)})

# Environmental change per species in original distribution
# considering a fixed geographical distribution
lec <- foreach(i=1:100) %dopar% {
  ec <- data.frame(sp=rep(names(bbs.env[,3:(3+nsp-1)]),each=nper-1), 
                   time=rep(period[-1],nsp), ec=NA)
  for(s in 3:(3+nsp-1)) {
    try({
      spi <- bbs.env[,s]
      env.spi <- env.global[spi>0,]
      scores.spi <- suprow(pca.env,env.spi)$li
      #rs <- di1$site_id[di1[,s]>0]
      rs <- lsrand[[s-2]][i,]
      spi1 <- di1[di1$site_id %in% rs,s]
      env1 <- di1[di1$site_id %in% rs,(nsp+3):ncol(di1)]
      env1 <- env1[rep(1:nrow(env1),spi1),] # unfold abundance to multiple presences in the same site
      scores.spi1 <- suprow(pca.env,env1)$li
      grid.env1 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                         glob1=scores.spi,
                                         sp=scores.spi1, R=100, th.sp=0)
      for(p in period[-1]) {
        dip <- bbs.env[bbs.env$time==p,]
        spip <- spi1
        envp <- dip[dip$site_id %in% rs,(nsp+3):ncol(dip)]
        envp <- envp[rep(1:nrow(envp),spip),]
        scores.spip <- suprow(pca.env,envp)$li
        grid.envp <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                           glob1=scores.spi,
                                           sp=scores.spip, R=100, th.sp=0)
        # Compute Schoener's D
        D <- ecospat.niche.overlap (grid.env1, grid.envp, cor = TRUE)$D
        ec[ec$sp==names(bbs.env)[s] & ec$time==p, "ec"] <- D
      }
    })
  }
  ec
}

#save(lec,file="lec.RData")
leci <- lec[[1]]$ec
for(i in 1:length(lec)) leci <- cbind(leci,lec[[i]]$ec)
ec.mean <- apply(leci,1,mean,na.rm=T)
ec <- data.frame(sp=rep(names(bbs.env[,3:(3+nsp-1)]),each=nper-1), 
                 time=rep(period[-1],nsp), ec=ec.mean)

# Estimate temporal trends of Denv
ec.trends <- data.frame(a.ec=rep(NA,nsp),b.ec=rep(NA,nsp), p.ec=rep(NA,nsp))
for(s in 1:nsp){
  try({
    lmi <- lm(ec ~ time, data=ec[ec$sp==unique(ec$sp)[s],])
    ec.trends[s,1:2] <- coef(lmi)
    ec.trends[s,3] <- summary(lmi)$coefficients[2,4]
  })
}

# Add species names
ec.trends$aou <- as.integer(sub('.', '', unique(ec$sp)))
sp.traits <- left_join(sp.traits, ec.trends, by="aou")



# Niche change due to geographical change alone (fix environment)
lgc <- foreach(i=1:100) %dopar% {
  gc <- data.frame(sp=rep(names(bbs.env[,3:(3+nsp-1)]),each=nper-1), 
                     time=rep(period[-1],nsp), gc=NA)
  for(s in 3:(3+nsp-1)) {
    try({
      spi <- bbs.env[,s]
      env.spi <- env.global[spi>0,]
      scores.spi <- suprow(pca.env,env.spi)$li
      #rs <- di1$site_id[di1[,s]>0]
      rs <- lsrand[[s-2]][i,]
      spi1 <- di1[di1$site_id %in% rs,s]
      env0 <- di1[di1$site_id %in% rs,(nsp+3):ncol(di1)]
      env1 <- env0[rep(1:nrow(env0),spi1),] # unfold abundance to multiple presences in the same site
      scores.spi1 <- suprow(pca.env,env1)$li
      grid.env1 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                         glob1=scores.spi,
                                         sp=scores.spi1, R=100, th.sp=0)
      for(p in period[-1]) {
        dip <- bbs.env[bbs.env$time==p,]
        spip <- dip[dip$site_id %in% rs,s]
        envp <- env0[rep(1:nrow(env0),spip),]
        scores.spip <- suprow(pca.env,envp)$li
        grid.envp <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                           glob1=scores.spi,
                                           sp=scores.spip, R=100, th.sp=0)
        # Compute Schoener's D
        D <- ecospat.niche.overlap (grid.env1, grid.envp, cor = TRUE)$D
        gc[gc$sp==names(bbs.env)[s] & gc$time==p, "gc"] <- D
      }
    })
  }
  gc
}

stopCluster(cl)


lgci <- lgc[[1]]$gc
for(i in 1:length(lgc)) lgci <- cbind(lgci,lgc[[i]]$gc)
gc.mean <- apply(lgci,1,mean,na.rm=T)
gc <- data.frame(sp=rep(names(bbs.env[,3:(3+nsp-1)]),each=nper-1), 
                 time=rep(period[-1],nsp), gc=gc.mean)

# Estimate temporal trends of Denv
gc.trends <- data.frame(a.gc=rep(NA,nsp),b.gc=rep(NA,nsp),p.gc=rep(NA,nsp))
for(s in 1:nsp){
  try({
    lmi <- lm(gc ~ time, data=gc[gc$sp==unique(gc$sp)[s],])
    gc.trends[s,1:2] <- coef(lmi)
    gc.trends[s,3] <- summary(lmi)$coefficients[2,4]
  })
}

# Add species names
gc.trends$aou <- as.integer(sub('.', '', unique(gc$sp)))
sp.traits <- left_join(sp.traits, gc.trends, by="aou")

save(lec,lgc,sp.traits,file="Niche_components.RData")

#############################################################################
#############################################################################

# 1. Calculate niche change


#------------------------------

# Shoener's D metric
# D is calculated by comparing the first period with all subsequent periods

cl <- makeCluster(8)
registerDoParallel(cl)
clusterExport(cl, list("nsp","period","nper","bbs.env","di1","lsrand","env.global","pca.env","scores.globclim"))
clusterCall(cl,function(){library(ecospat);library(ade4)})


lDres <- foreach(i=1:100) %dopar% {
  Dres <- data.frame(sp=rep(names(bbs.env[,3:(3+nsp-1)]),each=nper-1), 
                     time=rep(period[-1],nsp), D=NA)
  for(s in 3:(3+nsp-1)) {
    try({
      spi <- bbs.env[,s]
      env.spi <- env.global[spi>0,]
      scores.spi <- suprow(pca.env,env.spi)$li
      #rs <- di1$site_id[di1[,s]>0]
      rs <- lsrand[[s-2]][i,]
      spi1 <- di1[di1$site_id %in% rs,s]
      env1 <- di1[di1$site_id %in% rs,(nsp+3):ncol(di1)]
      env1 <- env1[rep(1:nrow(env1),spi1),] # unfold abundance to multiple presences in the same site
      scores.spi1 <- suprow(pca.env,env1)$li
      grid.env1 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                         glob1=scores.spi,
                                         sp=scores.spi1, R=100, th.sp=0)
      for(p in period[-1]) {
        dip <- bbs.env[bbs.env$time==p,]
        spip <- dip[dip$site_id %in% rs,s]
        envp <- dip[dip$site_id %in% rs,(nsp+3):ncol(dip)]
        envp <- envp[rep(1:nrow(envp),spip),]
        scores.spip <- suprow(pca.env,envp)$li
        grid.envp <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                           glob1=scores.spi,
                                           sp=scores.spip, R=100, th.sp=0)
        # Compute Schoener's D
        D <- ecospat.niche.overlap (grid.env1, grid.envp, cor = TRUE)$D
        Dres[Dres$sp==names(bbs.env)[s] & Dres$time==p, "D"] <- D
      }
    })
  }
  Dres
}

stopCluster(cl)


lDresi <- lDres[[1]]$D
for(i in 1:length(lDres)) lDresi <- cbind(lDresi,lDres[[i]]$D)
D.mean <- apply(lDresi,1,mean,na.rm=T)
Dres <- data.frame(sp=rep(names(bbs.env[,3:(3+nsp-1)]),each=nper-1), 
                 time=rep(period[-1],nsp), D=D.mean)

# Estimate temporal trends of D
Dtrends <- data.frame(a.D=rep(NA,nsp),b.D=rep(NA,nsp),se.D=rep(NA,nsp), p.D=rep(NA,nsp))
for(s in 1:nsp){
  try({
    lmi <- lm(D ~ time, data=Dres[Dres$sp==unique(Dres$sp)[s],])
    Dtrends[s,1:2] <- coef(lmi)
    Dtrends[s,3] <- summary(lmi)$coefficients[2,2]
    Dtrends[s,4] <- summary(lmi)$coefficients[2,4]
  })
}

# Add species names
Dtrends$aou <- as.integer(sub('.', '', unique(Dres$sp)))

hist(Dtrends$b.D)
summary(Dtrends$b.D)
hist(Dtrends$p.D)
summary(Dtrends$p.D)
hist(Dtrends$b.D[Dtrends$p.D<0.05])
length(Dtrends$b.D[Dtrends$p.D<0.05])


#------------------------------

# Niche breadth and position
# For each period:
#' 1. Build data.frame with PCs and species abundances as columns (see help)
#' 2. Apply function ecospat.nichePOSNB

# Calculate metrics
sp.env <- cbind(bbs.env[,1:2], scores.globclim, bbs.env[,3:(3+nsp-1)])
# subsample by assigning 0 to routes left out
lsp.env <- list()
for(i in 1:100){
  sp.env.i <- sp.env
  for(s in 5:(5+nsp-1)) sp.env.i[!(sp.env.i$site_id %in% lsrand[[s-4]][i,]),s] <- 0
  lsp.env[[i]] <- sp.env.i
}


llniche <- list()
llnb <- list()
for(p in 1:nper){
  lniche <- list()
  lnb <- list()
  for(i in 1:100){
    dti <- lsp.env[[i]]
    dtp <- dti[dti$time==period[p],]
    nichei <- ecospat.nichePOSNB(dtp, 3:4, 5:ncol(dtp))
    lniche[[i]] <- nichei[,1:2]
    lnb[[i]] <- ecospat.nicheNBmean(nichei)
  }
  llniche[[p]] <- lniche
  llnb[[p]] <- lnb
}

#save(llniche,llnb,file="lniche.RData")

lnb <- list()
for(p in 1:nper){
  mp <- do.call("cbind",llnb[[p]])
  lnb[[p]] <- apply(mp,1,mean,na.rm=T)
} 

lniche <- list()
for(p in 1:nper){
  lniche[[p]] <- as.data.frame(aaply(laply(llniche[[p]], as.matrix), 
                                     c(2, 3), mean, na.rm=T))
} 



# Temporal trends of niche breadth and position
# Calculate summary statistics per species:
#' trend of niche breath (linear regression)
#' persistence (straightness) of niche position (angle differences)
#' distance between first and last position

# Niche breadth
nbd <- as.data.frame(do.call("cbind",lnb))
colnames(nbd) <- period
nb.trends <- data.frame(a.nb=numeric(0),b.nb=numeric(0),se.nb=numeric(0),p.nb=numeric(0))
for(i in 1:nrow(nbd)){
  try({
    lmi <- lm(as.numeric(nbd[i,]) ~ period)
    nb.trends[i,1:2] <- coef(lmi)
    nb.trends[i,3] <- summary(lmi)$coefficients[2,2]
    nb.trends[i,4] <- summary(lmi)$coefficients[2,4]
  })
}
# Add species
nb.trends$aou <- as.integer(sub('.', '', row.names(nbd)))


# Niche position

# Trends of distance to first position
dist.trends <- data.frame(a.dist=numeric(0),b.dist=numeric(0),se.dist=numeric(0),p.dist=numeric(0))
for(s in 1:nsp){
  ms <- do.call(rbind, (lapply(lniche, function(x) x[s,])))
  ds <- as.matrix(dist(ms))[-1,1]
  lmi <- lm(ds ~ period[-1])
  dist.trends[s,1:2] <- coef(lmi)
  dist.trends[s,3] <- summary(lmi)$coefficients[2,2]
  dist.trends[s,4] <- summary(lmi)$coefficients[2,4]
}
# Add species
dist.trends$aou <- as.integer(sub('.', '', row.names(lniche[[1]])))


# Straightness and net distance
library(trajr)
npos.trends <- data.frame(straightness=numeric(0), dist=numeric(0))
centroid.pc <- data.frame(pc1=numeric(0),pc2=numeric(0))
pos1.time <- data.frame(period=period)
pos2.time <- data.frame(period=period)
for(i in 1:nsp){
  x.pos <- unlist(lapply(lniche, function(x) x[i,1]))
  y.pos <- unlist(lapply(lniche, function(x) x[i,2]))
  centroid.pc[i,] <- c(mean(x.pos), mean(y.pos))
  coords <- data.frame(x = x.pos, y = y.pos, times = period)
  pos1.time[,i+1] <- x.pos
  pos2.time[,i+1] <- y.pos
  trj <- TrajFromCoords(coords)
  # plot(trj)
  npos.trends[i,1] <- TrajStraightness(trj)
  npos.trends[i,2] <- TrajDistance(trj)
}

# Per species niche position change
pos1.linear <- data.frame(a.pc1=numeric(nsp),b.pc1=numeric(nsp),p.pc1=numeric(nsp))
pos2.linear <- data.frame(a.pc2=numeric(nsp),b.pc2=numeric(nsp),p.pc2=numeric(nsp))
for(i in 2:ncol(pos1.time)){
  lm1 <- lm(pos1.time[,i] ~ period)
  lm2 <- lm(pos2.time[,i] ~ period)
  pos1.linear[i-1,] <- c(coef(lm1),summary(lm1)$coefficients[2,4])
  pos2.linear[i-1,] <- c(coef(lm2),summary(lm2)$coefficients[2,4])
}

pos1.linear$aou <- as.integer(sub('.', '', row.names(lniche[[1]])))
pos2.linear$aou <- as.integer(sub('.', '', row.names(lniche[[1]])))


# Join niche breadth and position metrics
centroid.pc$aou <- as.integer(sub('.', '', row.names(lniche[[1]])))
npos.trends$aou <- as.integer(sub('.', '', row.names(lniche[[1]])))
names(pos1.time)[-1] <- as.integer(sub('.', '', row.names(lniche[[1]])))
names(pos2.time)[-1] <- as.integer(sub('.', '', row.names(lniche[[1]])))
niche <- left_join(Dtrends, nb.trends, by="aou")
niche <- left_join(niche, dist.trends, by="aou")
niche <- left_join(niche, npos.trends, by="aou")
niche <- left_join(niche, pos1.linear, by="aou")
niche <- left_join(niche, pos2.linear, by="aou")
niche <- left_join(niche, centroid.pc, by="aou")


#############################################################################

# Calculate geographical change
# Methodologically analogous to niche change

# Estimate geographical overlap (Dgeo) by correlating abundance in each route between two years?
# size ("breadth") with occupancy
# position with weighted centroid

# As in https://plantarum.ca/2021/12/02/schoenersd/


cl <- makeCluster(8)
registerDoParallel(cl)
clusterExport(cl, list("nsp","period","nper","bbs.geo","lsrand"))
clusterCall(cl,function(){library(ecospat)})


lDgeo <- foreach(i=1:100) %dopar% {
  Dgeo <- data.frame(sp=rep(names(bbs.geo[c(3:(3+nsp-1))]),each=nper-1), 
                     time=rep(period[-1],nsp), Dgeo=NA)
  for(s in 3:(3+nsp-1)){
    rs <- lsrand[[s-2]][i,]
    spi1 <- bbs.geo[bbs.geo$time==1981 & bbs.geo$site_id %in% rs,s]
    spi1 <- spi1/sum(spi1)
    for(p in period[-1]){
      spip <- bbs.geo[bbs.geo$time==p & bbs.geo$site_id %in% rs,s]
      spip <- spip/sum(spip)
      # Compute Schoener's D
      D <- 1 - sum(abs(spi1 - spip))/2
      Dgeo[Dgeo$sp==names(bbs.geo)[s] & Dgeo$time==p, "Dgeo"] <- D
    }
  }
  Dgeo
}
stopCluster(cl)

lDgeoi <- lDgeo[[1]]$Dgeo
for(i in 1:length(lDgeo)) lDgeoi <- cbind(lDgeoi,lDgeo[[i]]$Dgeo)
Dgeo.mean <- apply(lDgeoi,1,mean,na.rm=T)
Dgeo <- data.frame(sp=rep(names(bbs.geo[,3:(3+nsp-1)]),each=nper-1), 
                   time=rep(period[-1],nsp), Dgeo=Dgeo.mean)


# Estimate temporal trends of D
Dgeo.trends <- data.frame(a.D=rep(NA,nsp),b.D=rep(NA,nsp),p.D=rep(NA,nsp))
for(s in 1:nsp){
  try({
    lmi <- lm(Dgeo ~ time, data=Dgeo[Dgeo$sp==unique(Dgeo$sp)[s],])
    Dgeo.trends[s,1:2] <- coef(lmi)
    Dgeo.trends[s,3] <- summary(lmi)$coefficients[2,4]
  })
}

# Add species names
Dgeo.trends$aou <- as.integer(sub('.', '', unique(Dgeo$sp)))
#Dgeo.trends <- left_join(Dgeo.trends, sp.traits[,c("aou","species")], by="aou")
#Dgeo.trends$species <- gsub(' ', '_', Dgeo.trends$species)

hist(Dgeo.trends$b.D)
summary(Dgeo.trends$b.D)
hist(Dgeo.trends$p.D)
summary(Dgeo.trends$p.D)
hist(Dgeo.trends$b.D[Dgeo.trends$p.D<0.05])
length(Dgeo.trends$b.D[Dgeo.trends$p.D<0.05])


plot(Dgeo.trends$b.D,Dtrends$b.D)
points(Dgeo.trends$b.D[Dtrends$p.D<0.05],Dtrends$b.D[Dtrends$p.D<0.05],pch=16,col=3)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(a=0,b=1)

cor.test(Dgeo.trends$b.D,Dtrends$b.D)
cor.test(Dgeo.trends$b.D[Dtrends$p.D<0.05],Dtrends$b.D[Dtrends$p.D<0.05],pch=16,col=3)


#------------------------------

# Occupancy and position


sp.geo <- bbs.geo[,c(1:2,(ncol(bbs.geo)-1):ncol(bbs.geo),3:(3+nsp-1))]
# subsample by assigning 0 to routes left out
lsp.geo <- list()
for(i in 1:100){
  sp.geo.i <- sp.geo
  for(s in 5:(5+nsp-1)) sp.geo.i[!(sp.geo.i$site_id %in% lsrand[[s-4]][i,]),s] <- 0
  lsp.geo[[i]] <- sp.geo.i
}

# Calculate metrics
llrange <- list()
llocc <- list()
for(p in 1:nper){
  lrange <- list()
  locc <- list()
  for(i in 1:100){
    dti <- lsp.geo[[i]]
    dtp <- dti[dti$time==period[p],]
    rangei <- ecospat.nichePOSNB(dtp, 3:4, 5:ncol(dtp))
    lrange[[i]] <- rangei[,1:2]
    occ <- c()
    for(s in 5:ncol(dtp)){
      dtp.sp <- dtp[,s]
      occ[s-4] <- length(dtp.sp[dtp.sp>0])/length(dtp.sp)
    }
    occ <- as.data.frame(occ)
    row.names(occ) <- names(dti[,5:ncol(dti)])
    locc[[i]] <- occ
  }
  llrange[[p]] <- lrange
  llocc[[p]] <- locc
}


locc <- list()
for(p in 1:nper){
  mp <- do.call("cbind",llocc[[p]])
  locc[[p]] <- apply(mp,1,mean,na.rm=T)
} 

lrange <- list()
for(p in 1:nper){
  lrange[[p]] <- as.data.frame(aaply(laply(llrange[[p]], as.matrix), 
                                     c(2, 3), mean, na.rm=T))
} 



# Temporal trends of range occupancy and position

# Occupancy
occd <- do.call("cbind",locc)
colnames(occd) <- period
occ.trends <- data.frame(a.occ=numeric(0),b.occ=numeric(0),p.occ=numeric(0))
for(i in 1:nrow(occd)){
  lmi <- lm(as.numeric(occd[i,]) ~ period)
  occ.trends[i,1:2] <- coef(lmi)
  occ.trends[i,3] <- summary(lmi)$coefficients[2,4]
}

# Add species
occ.trends$aou <- as.integer(sub('.', '', row.names(occd)))

hist(occ.trends$b.occ,xlab="Occupancy",cex.lab=1.2, font.lab=2 )
summary(occ.trends$b.occ)
hist(occ.trends$p.occ)
summary(occ.trends$p.occ)
hist(occ.trends$b.occ[occ.trends$p.occ<0.05])
summary(occ.trends$b.occ[occ.trends$p.occ<0.05])


# Geographical position

# Trends of distance to first position
gdist.trends <- data.frame(a.gdist=numeric(0),b.gdist=numeric(0),se.gdist=numeric(0),p.gdist=numeric(0))
for(s in 1:nsp){
  ms <- do.call(rbind, (lapply(lrange, function(x) x[s,])))
  coords <- st_as_sf(ms, coords = c("Longitude_pos", "Latitude_pos"), crs = 4326)
  coords <- st_transform(coords, crs=3310)
  ds <- as.matrix(st_distance(coords))[-1,1]/1000
  lmi <- lm(as.numeric(ds) ~ period[-1])
  gdist.trends[s,1:2] <- coef(lmi)
  gdist.trends[s,3] <- summary(lmi)$coefficients[2,2]
  gdist.trends[s,4] <- summary(lmi)$coefficients[2,4]
}
# Add species
gdist.trends$aou <- as.integer(sub('.', '', row.names(lrange[[1]])))


# Straightness and net distance
library(trajr)
gpos.trends <- data.frame(straightness=numeric(0), dist=numeric(0))
centroid.geo <- data.frame(lat=numeric(0),lon=numeric(0))
gpos1.time <- data.frame(period=period)
gpos2.time <- data.frame(period=period)
for(i in 1:nsp){
  x.pos <- unlist(lapply(lrange, function(x) x[i,2]))
  y.pos <- unlist(lapply(lrange, function(x) x[i,1]))
  centroid.geo[i,] <- c(mean(x.pos), mean(y.pos))
  coords0 <- data.frame(x = x.pos, y = y.pos)
  coords0 <- st_as_sf(coords0, coords = c("x", "y"), crs = 4326)
  coords0 <- st_transform(coords0, crs=3310)
  coords <- do.call(rbind, st_geometry(coords0)) %>% 
    as.data.frame() %>% setNames(c("x","y"))
  coords <- data.frame(x = coords$x, y = coords$y, times = period)
  gpos1.time[,i+1] <- x.pos
  gpos2.time[,i+1] <- y.pos
  trj <- TrajFromCoords(coords)
  # plot(trj)
  gpos.trends[i,1] <- TrajStraightness(trj)
  gpos.trends[i,2] <- TrajDistance(trj)
}

# Join niche breadth and position metrics
centroid.geo$aou <- as.integer(sub('.', '', row.names(lrange[[1]])))
gpos.trends$aou <- as.integer(sub('.', '', row.names(lrange[[1]])))
names(gpos1.time)[-1] <- as.integer(sub('.', '', row.names(lrange[[1]])))
names(gpos2.time)[-1] <- as.integer(sub('.', '', row.names(lrange[[1]])))
geospace <- left_join(Dgeo.trends, occ.trends, by="aou")
geospace <- left_join(geospace, gdist.trends, by="aou")
geospace <- left_join(geospace, gpos.trends, by="aou")
names(geospace) <- c("a.Dgeo","b.Dgeo","p.Dgeo","aou","a.occ","b.occ","p.occ",
                     "a.gdist","b.gdist","se.gdist","p.gdist","straightness.geo","dist.geo")




save(Dcli, Dlu, Dres, Dgeo, cluc, niche, geospace, centroid.pc, centroid.geo,
     pos1.time, pos2.time, gpos1.time, gpos2.time, file="Temporal_trends.RData")




#############################################################################
#############################################################################

# Niche drivers (species traits and other attributes)
setwd("/Users/viana/Dropbox/TFM_birds_BBS/Data")
load("Temporal_trends.RData")

niche.drivers <- dplyr::left_join(sp.traits,niche,by="aou")
niche.drivers <- dplyr::left_join(niche.drivers,geospace,by="aou")

# Variable transformation
niche.drivers$SSI <- as.numeric(niche.drivers$SSI)
niche.drivers$dist <- sqrt(niche.drivers$dist)
niche.drivers$Hand.Wing.Index <- log(niche.drivers$Hand.Wing.Index)
niche.drivers$Mass <- log(niche.drivers$Mass)
niche.drivers$range.size <- log(niche.drivers$range.size)
#niche.drivers$Migration <- as.factor(niche.drivers$Migration)
niche.drivers$Habitat <- as.factor(niche.drivers$Habitat)
niche.drivers$Trophic.Level <- as.factor(niche.drivers$Trophic.Level)
niche.drivers$fHabitat.Density <- as.factor(niche.drivers$Habitat.Density)
niche.drivers$range.lat <- log(niche.drivers$Max.Latitude-niche.drivers$Min.Latitude)
niche.drivers$Min.Latitude <- abs(niche.drivers$Min.Latitude)
niche.drivers$Max.Latitude <- abs(niche.drivers$Max.Latitude)
niche.drivers$Territoriality[niche.drivers$Territoriality=="none"] <- 0
niche.drivers$Territoriality[niche.drivers$Territoriality=="weak"] <- 1
niche.drivers$Territoriality[niche.drivers$Territoriality=="strong"] <- 2
niche.drivers$Territoriality <- as.integer(niche.drivers$Territoriality)



# Some results
nrow(niche.drivers[niche.drivers$p.D<0.05 & niche.drivers$b.D<0,])
nrow(niche.drivers[niche.drivers$p.D<0.05 & niche.drivers$b.D>0,])

nrow(niche.drivers[niche.drivers$p.nb<0.05 & niche.drivers$b.nb<0,])
nrow(niche.drivers[niche.drivers$p.nb<0.05 & niche.drivers$b.nb>0,])

nrow(niche.drivers[niche.drivers$p.dist<0.05 & niche.drivers$b.dist<0,])
nrow(niche.drivers[niche.drivers$p.dist<0.05 & niche.drivers$b.dist>0,])


# Variation in D explained by environmental change in initial distribution
mec <- lm(b.D ~ b.ec, data=niche.drivers)
summary(mec) # half of the variation explained by env change in original distribution
mgc <- lm(b.D ~ b.gc, data=niche.drivers)
summary(mgc)
megc <- lm(b.D ~ b.ec + b.gc, data=niche.drivers)
summary(megc)

# Environmental change is not related to geographical change
cor.test(niche.drivers$b.Dgeo,niche.drivers$b.ec)



table1 <- niche.drivers[,c("Species1","b.D","b.nb","b.dist","dist")]
write.table(table1, file="Table_S1.txt", row.names=F, sep="\t")



#----------------------------------------------------------

# Phylogenetic signal

# Check if taxonomy matches
birdtree.sp <- read.csv("BLIOCPhyloMasterTax.csv")
tax_match <- read.csv("AVONET_BirdTree_match.csv")
spd <- left_join(niche.drivers[,c("species","Species1")],tax_match[,c("Species1","Species3")], by="Species1")
which(duplicated(spd$Species1)) # one duplicated species
spd <- spd[spd$Species3!="Dendroica petechia",] # keep the aestiva (North American subspecies)
all(spd$Species3 %in% birdtree.sp$Scientific)



# From https://github.com/nicholasjclark/BBS.occurrences/blob/master/Clark_etal_analysis/Appendix_S4_PhyloTraitData.Rmd

#Export the vector as a one-column .csv file 
#write.csv(spd$Species3, "sp.names.csv", row.names = FALSE)

# Once exported, copy the species names and paste into the `Select species` form 
# on the [Phylogeny Subsets](http://birdtree.org/subsets/) page at Birdtree.org. 
# We downloaded 100 trees from the *Ericsson All Species* dataset for our analyses. 
# Once processed, save the resulting .nex file and read in the multiphylo object 
# using functions in the `ape` package
setwd("/Users/viana/Dropbox/TFM_birds_BBS/Data/Phylo_100_BirdTree")
sp.trees <- ape::read.nexus("output.nex")
setwd("/Users/viana/Dropbox/TFM_birds_BBS/Data")
# Now check to make sure that all of the species names are represented in the tree. 
# Here, we have to include the underscore once again, as this is included in Birdtree.org phylogeny subsets. 
# This call should return `TRUE` if there are no unmatched names
sp.names.underscore <- gsub(' ', '_', spd$Species3)
all(sp.names.underscore %in% sp.trees[[1]]$tip.label) # Circus hudsonius is not in the tree!

# Add phylo species names
names(spd)[3] <- "phylo"
niche.drivers <- left_join(niche.drivers,spd[,c("Species1","phylo")], by="Species1")
niche.drivers$phylo <- gsub(' ', '_', niche.drivers$phylo)

# Phylogenetic covariance
#cons.tree <- sp.trees[[1]] # 1 phylogenetic tree
# From https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
# The phylo object contains information on the relationship between species. 
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).

# Build a consensus tree
library(phytools)
cons.tree <- consensus.edges(sp.trees,method="mean.edge")
#cons.tree <- di2multi(cons.tree)
A <- ape::vcv.phylo(cons.tree)

# Phylogenetic signal
library(pez)
trait.y<-setNames(niche.drivers$b.D,niche.drivers$phylo)
phylosig(cons.tree, trait.y, method = "lambda", test = TRUE)


#----------------------------------------------------------

# Run models
library(brms)
library(car)
# Based on
# https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html

niche.drivers2 <- niche.drivers
niche.drivers2$obs <- 1:nrow(niche.drivers2)
niche.drivers2$b.D <- niche.drivers2$b.D*100
niche.drivers2$se.D <- niche.drivers2$se.D*100
niche.drivers2$Hand.Wing.Index <- scale(niche.drivers2$Hand.Wing.Index)
niche.drivers2$Mass <- scale(niche.drivers2$Mass)
niche.drivers2$Habitat.Density <- scale(niche.drivers2$Habitat.Density)
niche.drivers2$Migration <- scale(niche.drivers2$Migration)
niche.drivers2$range.size <- scale(niche.drivers2$range.size)
niche.drivers2$Min.Latitude <- scale(niche.drivers2$Min.Latitude)
niche.drivers2$SSI <- scale(niche.drivers2$SSI)
niche.drivers2$Territoriality <- scale(niche.drivers2$Territoriality)

# Check for collinearity
vif(lm(b.D ~ Hand.Wing.Index + Mass + Habitat.Density +  
      Migration + range.size + Min.Latitude + SSI + Territoriality,
      data=niche.drivers2))

#------------------------------

# Models for species traits

# Niche overlap
m1.brm <- brm(b.D | se(se.D)  ~ Hand.Wing.Index + Mass + Habitat.Density +  
              Migration + range.size + Min.Latitude + SSI + Territoriality +
              (1|gr(phylo, cov = A)) + (1|obs),
              data=niche.drivers2, data2 = list(A = A),
              prior = c(prior(normal(0, 10), "Intercept"),
                        prior(student_t(3, 0, 10), "sd")),
              control = list(adapt_delta = .95), 
              chains=3, thin=10, iter=4000, warmup=1000, cores=3)


print(summary(m1.brm), digits = 3)
plot(m1.brm)
pp_check(m1.brm)


# Niche breadth
niche.drivers2$b.nb <- niche.drivers2$b.nb*100
niche.drivers2$se.nb <- niche.drivers2$se.nb*100

m2.brm <- brm(b.nb | se(se.nb)  ~ Hand.Wing.Index + Mass + Habitat.Density +  
                Migration + range.size + Min.Latitude + SSI + Territoriality +
                (1|gr(phylo, cov = A)) + (1|obs),
              data=niche.drivers2, data2 = list(A = A),
              prior = c(prior(normal(0, 10), "Intercept"),
                        prior(student_t(3, 0, 10), "sd")),
              control = list(adapt_delta = .95), 
              chains=3, thin=10, iter=4000, warmup=1000, cores=3)
print(summary(m2.brm), digits = 3)

# Niche position: distance trends
m3.brm <- brm(b.dist | se(se.dist)  ~ Hand.Wing.Index + Mass + Habitat.Density +  
                Migration + range.size + Min.Latitude + SSI + Territoriality +
                (1|gr(phylo, cov = A)) + (1|obs),
              data=niche.drivers2, data2 = list(A = A),
              prior = c(prior(normal(0, 10), "Intercept"),
                        prior(student_t(3, 0, 10), "sd")),
              control = list(adapt_delta = .95), 
              chains=3, thin=10, iter=4000, warmup=1000, cores=3)
print(summary(m3.brm), digits = 3)

# Niche position: net distance
hist(niche.drivers2$dist)

m4.brm <- brm(dist  ~ Hand.Wing.Index + Mass + Habitat.Density +  
                Migration + range.size + Min.Latitude + SSI + Territoriality +
                (1|gr(phylo, cov = A)),
              data=niche.drivers2, data2 = list(A = A),
              prior = c(prior(normal(0, 10), "Intercept"),
                        prior(student_t(3, 0, 10), "sd")),
              control = list(adapt_delta = .95), 
              chains=3, thin=10, iter=4000, warmup=1000, cores=3)
print(summary(m4.brm), digits = 3)

# Niche position: straightness
m5.brm <- brm(straightness  ~ Hand.Wing.Index + Mass + Habitat.Density +  
                Migration + range.size + Min.Latitude + SSI + Territoriality +
                (1|gr(phylo, cov = A)),
              data=niche.drivers2, data2 = list(A = A),
              prior = c(prior(normal(0, 10), "Intercept"),
                        prior(student_t(3, 0, 10), "sd")),
              control = list(adapt_delta = .97), 
              chains=3, thin=10, iter=4000, warmup=1000, cores=3)
print(summary(m5.brm), digits = 3)


#------------------------------

# Models for conservation status
# Cumulative (ordinal) models
# Conservation as a response variable


niche.drivers2$b.D.sc <- scale(niche.drivers2$b.D)
niche.drivers2$b.nb.sc <- scale(niche.drivers2$b.nb)
niche.drivers2$b.dist.sc <- scale(niche.drivers2$b.dist)

# Population trend
# with phylogenetic correlation
m21.brm <- brm(PT.c  ~ b.D.sc + b.nb.sc + b.dist.sc +
                 (1|gr(phylo, cov = A)), family= cumulative,
               data=niche.drivers2, data2 = list(A = A),
               control = list(adapt_delta = .97), 
               chains=4, thin=10, iter=5000, warmup=1000, cores=4)
print(summary(m21.brm), digits = 3)
# without phylogenetic correlation
m22.brm <- brm(PT.c  ~ b.D.sc + b.nb.sc + b.dist.sc, family= cumulative,
               data=niche.drivers2, 
               control = list(adapt_delta = .95), 
               chains=4, thin=10, iter=5000, warmup=1000, cores=4)
print(summary(m22.brm), digits = 3)

# Threats to breeding
# with phylogenetic correlation
m23.brm <- brm(TB.c  ~ b.D.sc + b.nb.sc + b.dist.sc +
                 (1|gr(phylo, cov = A)), family= cumulative,
               data=niche.drivers2, data2 = list(A = A),
               control = list(adapt_delta = .97), 
               chains=4, thin=10, iter=5000, warmup=1000, cores=4)
print(summary(m23.brm), digits = 3)
# without phylogenetic correlation
m24.brm <- brm(TB.c  ~ b.D.sc + b.nb.sc + b.dist.sc, family= cumulative,
               data=niche.drivers2, 
               control = list(adapt_delta = .97), 
               chains=4, thin=10, iter=5000, warmup=1000, cores=4)
print(summary(m24.brm), digits = 3)


# Threats to non-breeding
# with phylogenetic correlation
m25.brm <- brm(TN.c  ~ b.D.sc + b.nb.sc + b.dist.sc +
                 (1|gr(phylo, cov = A)), family= cumulative,
               data=niche.drivers2, data2 = list(A = A),
               control = list(adapt_delta = .97), 
               chains=4, thin=10, iter=5000, warmup=1000, cores=4)
print(summary(m25.brm), digits = 3)
# without phylogenetic correlation
m26.brm <- brm(TN.c  ~ b.D.sc + b.nb.sc + b.dist.sc +
                 (1|gr(phylo, cov = A)), family= cumulative,
               data=niche.drivers2, data2 = list(A = A),
               control = list(adapt_delta = .97), 
               chains=4, thin=10, iter=5000, warmup=1000, cores=4)
print(summary(m26.brm), digits = 3)




save(m1.brm, m2.brm, m3.brm, m4.brm, m5.brm, m21.brm, m22.brm, m23.brm, 
     m24.brm, m25.brm, m26.brm,file="brm_models.RData")




#############################################################################
#############################################################################

# Figures

library(ggplot2)
library(cowplot)
library(egg)
library(scales)
library(viridis)
library(patchwork)
theme_set(theme_classic())

#----------------------------------------------------------

# Figure 1

library(maps)
library("mapdata")
library(RColorBrewer)
pal <-  colorRampPalette(c("violetred4", "orange3","olivedrab3"))



library(factoextra)
quartz(height=4,width=4)
par(mar=c(0,0,0,0))
fviz_pca_biplot(pca.env, label="var", repel=TRUE, col.var = "black", col.ind = "grey", pointsize=0.3, title = "", xlab="PC1",ylab="PC2") +
  theme(panel.grid.major=element_blank())+
  theme(panel.border=element_rect(colour="black", fill=NA))+
  theme(axis.text = element_text(size=12, colour="black"))



# Prepare data
bbs.pc <- left_join(bbs.env, bbs[!duplicated(bbs$site_id),c("site_id","Longitude","Latitude")], by="site_id")
bbs.pc <- cbind(bbs.pc,scores.globclim)
bbs.pc <- bbs.pc %>% 
  group_by(site_id) %>%
  summarise(PC1=mean(Axis1),PC2=mean(Axis2),Longitude=unique(Longitude),Latitude=unique(Latitude)) %>% 
  data.frame()
bbs.pc<-SpatialPointsDataFrame(bbs.pc[,4:5],data=bbs.pc[,1:3],
                               proj4string=CRS("+proj=longlat +datum=WGS84"))

# Map sites
map('worldHires')
points(bbs.pc, col='red')
(e<-bbox(extent(bbs.pc)*1.2))

# Map PC1
quartz(height=4,width=7)
bbs.pc$order <- findInterval(bbs.pc$PC1, sort(bbs.pc$PC1))
legend_image <- as.raster(matrix(pal(nrow(bbs.pc)), nrow=1))
map('worldHires', xlim = e[1, ], ylim = e[2, ],mar=c(1.5,1,1.5,1))
#map.axes()
points(bbs.pc,col=pal(nrow(bbs.pc))[bbs.pc$order],pch=16,cex=2)
map('worldHires', xlim = e[1, ], ylim = e[2, ],add=T,lwd=2)
#polygon(x=c(-85,-60,-60,-85),y=c(28,28,31,31),col="white",border=NA)
rasterImage(legend_image, -77, 28, -58, 31)
text(x=seq(-77,-58,l=3), y=32 , labels = c(-7.7,2.7,13.2))
box(lwd=2)

# Map PC2
quartz(height=4,width=7)
bbs.pc$order <- findInterval(bbs.pc$PC2, sort(bbs.pc$PC2))
legend_image <- as.raster(matrix(pal(nrow(bbs.pc)), nrow=1))
map('worldHires', xlim = e[1, ], ylim = e[2, ],mar=c(1.5,1,1.5,1))
#map.axes()
points(bbs.pc,col=pal(nrow(bbs.pc))[bbs.pc$order],pch=16,cex=2)
map('worldHires', xlim = e[1, ], ylim = e[2, ],add=T,lwd=2)
#polygon(x=c(-85,-60,-60,-85),y=c(28,28,31,31),col="white",border=NA)
rasterImage(legend_image, -77, 28, -58, 31)
text(x=seq(-77,-58,l=3), y=32 , labels = c(-9.8,-3.2,3.4))
box(lwd=2)


#----------------------------------------------------------

# Figure 2

show_col(hue_pal()(2))

colf <- rep(0,nrow(niche.drivers))
colf[niche.drivers$p.D<0.05] <- "p<0.05"
colf[niche.drivers$p.D>=0.05] <- "p>=0.05"
niche.drivers$significance <- as.factor(colf)
gg.D <- ggplot(niche.drivers)  +
  geom_segment(aes(x = 1984, xend = 2018, y = a.D+b.D*1984, yend = a.D+b.D*2018, color=significance, alpha=I(0.3))) +
  geom_segment(aes(x = 1984, xend = 2018, y = mean(a.D)+mean(b.D)*1984, yend = mean(a.D)+mean(b.D)*2018), linewidth=2) +
  geom_segment(data=niche.drivers[niche.drivers$p.D<0.05,], 
               aes(x = 1984, xend = 2018, y = mean(a.D)+mean(b.D)*1984, yend = mean(a.D)+mean(b.D)*2018), linewidth=2, col="#F8766D") +
  xlim(1980, 2020) +
  ylim(0, 1) +
  xlab("Period") + ylab("Niche overlap (D)") +
  theme(legend.position = c(0.8, 0.9))

colf <- rep(0,nrow(niche.drivers))
colf[niche.drivers$p.nb<0.05] <- "p<0.05"
colf[niche.drivers$p.nb>=0.05] <- "p>=0.05"
niche.drivers$significance <- as.factor(colf)
gg.nb <- ggplot(niche.drivers)  +
  geom_segment(aes(x = 1981, xend = 2018, y = a.nb+b.nb*1981, yend = a.nb+b.nb*2018, color=significance, alpha=I(0.3))) +
  geom_segment(aes(x = 1981, xend = 2018, y = mean(a.nb)+mean(b.nb)*1981, yend = mean(a.nb)+mean(b.nb)*2018), linewidth=2) +
  geom_segment(data=niche.drivers[niche.drivers$p.nb<0.05,], 
               aes(x = 1981, xend = 2018, y = mean(a.nb)+mean(b.nb)*1981, yend = mean(a.nb)+mean(b.nb)*2018), linewidth=2, col="#F8766D") +
  xlim(1980, 2020) +
  ylim(0, 4) +
  xlab("Period") + ylab("Niche breadth") +
  theme(legend.position="none")

colf <- rep(0,nrow(niche.drivers))
colf[niche.drivers$p.dist<0.05] <- "p<0.05"
colf[niche.drivers$p.dist>=0.05] <- "p>=0.05"
niche.drivers$significance <- as.factor(colf)
gg.dist <- ggplot(niche.drivers)  +
  geom_segment(aes(x = 1984, xend = 2018, y = a.dist+b.dist*1984, yend = a.dist+b.dist*2018, color=significance, alpha=I(0.3))) +
  geom_segment(aes(x = 1984, xend = 2018, y = mean(a.dist)+mean(b.dist)*1984, yend = mean(a.dist)+mean(b.dist)*2018), linewidth=2) +
  geom_segment(data=niche.drivers[niche.drivers$p.dist<0.05,], 
               aes(x = 1984, xend = 2018, y = mean(a.dist)+mean(b.dist)*1984, yend = mean(a.dist)+mean(b.dist)*2018), linewidth=2, col="#F8766D") +
  xlim(1980, 2020) +
  ylim(0, 2) +
  xlab("Period") + ylab("Niche distance") +
  theme(legend.position="none")

# Figure
quartz(height=3.5,width=10)
ggarrange(gg.D, gg.nb, gg.dist, nrow=1, labels = c('(a)', '(b)', '(c)'))


#----------------------------------------------------------

# Figure 3

colf <- rep(0,nrow(niche.drivers))
colf[niche.drivers$p.pc1<0.05] <- "p<0.05"
colf[niche.drivers$p.pc1>=0.05] <- "p>=0.05"
niche.drivers$significance <- as.factor(colf)
gg.pc1 <- ggplot(niche.drivers)  +
  geom_segment(aes(x = 1981, xend = 2018, y = a.pc1+b.pc1*1981, yend = a.pc1+b.pc1*2018, color=significance, alpha=I(0.3))) +
  geom_segment(aes(x = 1981, xend = 2018, y = mean(a.pc1)+mean(b.pc1)*1981, yend = mean(a.pc1)+mean(b.pc1)*2018), linewidth=2) +
  geom_segment(data=niche.drivers[niche.drivers$p.pc1<0.05,],
               aes(x = 1981, xend = 2018, y = mean(a.pc1)+mean(b.pc1)*1981, yend = mean(a.pc1)+mean(b.pc1)*2018), linewidth=2, col="#F8766D") +
  xlim(1980, 2020) +
  ylim(-6, 4) +
  xlab("Period") + ylab("PC1") +
  theme(legend.position = c(0.85, 0.2))

colf <- rep(0,nrow(niche.drivers))
colf[niche.drivers$p.pc2<0.05] <- "p<0.05"
colf[niche.drivers$p.pc2>=0.05] <- "p>=0.05"
niche.drivers$significance <- as.factor(colf)
gg.pc2 <- ggplot(niche.drivers)  +
  geom_segment(aes(x = 1981, xend = 2018, y = a.pc2+b.pc2*1981, yend = a.pc2+b.pc2*2018, color=significance, alpha=I(0.3))) +
  geom_segment(aes(x = 1981, xend = 2018, y = mean(a.pc2)+mean(b.pc2)*1981, yend = mean(a.pc2)+mean(b.pc2)*2018), linewidth=2) +
  geom_segment(data=niche.drivers[niche.drivers$p.pc2<0.05,], 
               aes(x = 1981, xend = 2018, y = mean(a.pc2)+mean(b.pc2)*1981, yend = mean(a.pc2)+mean(b.pc2)*2018), linewidth=2, col="#F8766D") +
  xlim(1980, 2020) +
  ylim(-3, 3) +
  xlab("Period") + ylab("PC2") +
  theme(legend.position="none")

gg.d_ec <- ggplot(niche.drivers, aes(x = b.ec, y = b.D))  +
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Env. change in initial distribution") + ylab("Niche change")

gg.d_gc <- ggplot(niche.drivers, aes(x = b.gc, y = b.D))  +
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Env. change due to geographical change") + ylab("Niche change")



# Figure
quartz(height=6,width=7)
ggarrange(gg.pc1, gg.pc2, gg.d_ec, gg.d_gc, nrow=2, 
          labels = c('(a)','(b)','(c)','(d)'))


#----------------------------------------------------------

# Figure 4

gg.d_dgeo <- ggplot(niche.drivers, aes(x = b.Dgeo, y = b.D))  +
  geom_point() +
  geom_point(data=subset(niche.drivers, p.D<0.05), col="#F8766D") +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, col="black") +
  geom_smooth(data=subset(niche.drivers, p.D<0.05), method=lm, se=FALSE, fullrange=TRUE, col="#F8766D") +
  geom_vline(xintercept = 0, linetype="dashed") + 
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Geographical change") + ylab("Niche change (D trend)") 

gg.nb.occ <- ggplot(niche.drivers)  +
  geom_point(aes(x = b.occ, y = b.nb, col="All data")) +
  geom_point(data=subset(niche.drivers, p.nb<0.05), aes(x = b.occ, y = b.nb, col="p<0.05")) +
  geom_smooth(aes(x = b.occ, y = b.nb, col="p<0.05"), method=lm, se=FALSE, fullrange=TRUE) +
  geom_smooth(data=subset(niche.drivers, p.nb<0.05),aes(x = b.occ, y = b.nb, col="p<0.05"), method=lm, se=FALSE, fullrange=TRUE) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Occupancy change") + ylab("Niche breadth change") +
  scale_color_manual(name = "",
                     values = c("All data" = "black", "p<0.05" = "#F8766D"),
                     labels = c("All data", "p<0.05")) +
  theme(legend.position = c(0.8, 0.9))

gg.bdist.bgdist <- ggplot(niche.drivers, aes(x = b.gdist, y = b.dist))  +
  geom_point() +
  geom_point(data=subset(niche.drivers, p.dist<0.05), col="#F8766D") +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, col="black") +
  geom_smooth(data=subset(niche.drivers, p.dist<0.05), method=lm, se=FALSE, fullrange=TRUE, col="#F8766D") +
  geom_vline(xintercept = 0, linetype="dashed") + 
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Geographical distance trend") + ylab("Niche distance tend") 

gg.dist.gdist <- ggplot(niche.drivers, aes(x = log(dist.geo), y = dist))  +
  geom_point() +
  geom_point(data=subset(niche.drivers, p.D<0.05), col="#F8766D") +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, col="black") +
  geom_smooth(data=subset(niche.drivers, p.D<0.05), method=lm, se=FALSE, fullrange=TRUE, col="#F8766D") +
  xlab("Net geographical distance (log") + ylab("Net niche distance") 


quartz(height=6,width=7)
ggarrange(gg.d_dgeo, gg.nb.occ, gg.bdist.bgdist, gg.dist.gdist, nrow=2, 
          labels = c('(a)','(b)','(c)','(d)'),
          label.args=list(gp = grid::gpar(font = 4, cex = 1)))


#----------------------------------------------------------

# Figure 5
load("brm_models.RData")

brm.d <- summary(m1.brm)$fixed[-1,c(1,3,4)]
names(brm.d) <- c("Estimate","lCI","uCI")
brm.d$Predictor <- c("HWI","Mass","Habitat density","Migration","Range size","Min. Latitude","Habitat specialization","Territoriality")
gg.m.d <- ggplot(brm.d, aes(x=Estimate, y=Predictor)) + 
  geom_point()+
  geom_errorbar(aes(xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05)) +
  geom_vline(xintercept = 0)

brm.nb <- summary(m2.brm)$fixed[-1,c(1,3,4)]
names(brm.nb) <- c("Estimate","lCI","uCI")
brm.nb$Predictor <- c("HWI","Mass","Habitat density","Migration","Range size","Min. Latitude","Habitat specialization","Territoriality")
gg.m.nb <- ggplot(brm.nb, aes(x=Estimate, y=Predictor)) + 
  geom_point()+
  geom_errorbar(aes(xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05)) +
  geom_vline(xintercept = 0)

brm.bdist <- summary(m3.brm)$fixed[-1,c(1,3,4)]
names(brm.bdist) <- c("Estimate","lCI","uCI")
brm.bdist$Predictor <- c("HWI","Mass","Habitat density","Migration","Range size","Min. Latitude","Habitat specialization","Territoriality")
gg.m.bdist <- ggplot(brm.bdist, aes(x=Estimate, y=Predictor)) + 
  geom_point()+
  geom_errorbar(aes(xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05)) +
  geom_vline(xintercept = 0)

brm.dist <- summary(m4.brm)$fixed[-1,c(1,3,4)]
names(brm.dist) <- c("Estimate","lCI","uCI")
brm.dist$Predictor <- c("HWI","Mass","Habitat density","Migration","Range size","Min. Latitude","Habitat specialization","Territoriality")
gg.m.dist <- ggplot(brm.dist, aes(x=Estimate, y=Predictor)) + 
  geom_point()+
  geom_errorbar(aes(xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05)) +
  geom_vline(xintercept = 0)



levels(niche.drivers$fHabitat.Density) <- c("Dense","Semi-open","Open")
gg.d.hd <- ggplot(niche.drivers) +
  geom_boxplot(aes(x=fHabitat.Density, y=b.D), color="black", fill="lightgrey") +
  xlab("Habitat density") + ylab("Niche change (D trend)")

gg.d.lat <- ggplot(niche.drivers, aes(x = Min.Latitude, y = b.D))  +
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  xlab("Minimum latitude") + ylab("Niche change (D trend)")

gg.nb.ssi <- ggplot(niche.drivers, aes(x = SSI, y = b.nb))  +
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  xlab("Habitat specialisation") + ylab("Niche breadth trend")

niche.drivers$fTerritoriality <- factor(niche.drivers$Territoriality)
levels(niche.drivers$fTerritoriality) <- c("None","Weak","Strong")
gg.nb.ter <- ggplot(niche.drivers) +
  geom_boxplot(aes(x=fTerritoriality, y=b.nb), color="black", fill="lightgrey") +
  xlab("Territoriality") + ylab("Niche breadth trend")

gg.bdist.hd <- ggplot(niche.drivers) +
  geom_boxplot(aes(x=fHabitat.Density, y=b.dist), color="black", fill="lightgrey") +
  xlab("Habitat density") + ylab("Niche distance trend")

gg.dist.hd <- ggplot(niche.drivers) +
  geom_boxplot(aes(x=fHabitat.Density, y=dist), color="black", fill="lightgrey") +
  xlab("Habitat density") + ylab("Net niche distance")

gg.dist.rs <- ggplot(niche.drivers, aes(x = range.size, y = dist))  +
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  xlab("Range size (log)") + ylab("Net niche distance")


# Figure
quartz(height=9,width=8)
ggarrange(gg.m.d, gg.d.hd, gg.d.lat, gg.m.nb, gg.nb.ssi, gg.nb.ter,
          gg.m.bdist, gg.bdist.hd, ggplot() + theme_void(), gg.m.dist,
          gg.dist.hd, gg.dist.rs, nrow=4, 
          labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','','(i)','(j)','(k)'),
          label.args=list(gp = grid::gpar(font = 4, cex = 1)))

#----------------------------------------------------------

# Figure 6

brm.pt1 <- summary(m21.brm)$fixed[5:7,c(1,3,4)]
names(brm.pt1) <- c("Estimate_PT.c","lCI","uCI")
brm.pt1$Predictor <- c("Niche change (D trend)","Niche breadth trend","Niche distance trend")
brm.pt2 <- summary(m22.brm)$fixed[5:7,c(1,3,4)]
names(brm.pt2) <- c("Estimate_PT.c","lCI","uCI")
brm.pt2$Predictor <- c("Niche change (D trend)","Niche breadth trend","Niche distance trend")
quartz(height=3,width=4)
ggplot() + 
  geom_point(data=brm.pt1, aes(x=Estimate_PT.c, y=Predictor))+
  geom_errorbar(data=brm.pt1, aes(x=Estimate_PT.c, y=Predictor, xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05)) +
  geom_point(data=brm.pt2, aes(x=Estimate_PT.c, y=Predictor), col="#F8766D") +
  geom_errorbar(data=brm.pt2, aes(x=Estimate_PT.c, y=Predictor, xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05), col="#F8766D") +
  geom_vline(xintercept = 0)

brm.tbc1 <- summary(m23.brm)$fixed[4:6,c(1,3,4)]
names(brm.tbc1) <- c("Estimate_TB.c","lCI","uCI")
brm.tbc1$Predictor <- c("Niche change (D trend)","Niche breadth trend","Niche distance trend")
brm.tbc2 <- summary(m24.brm)$fixed[4:6,c(1,3,4)]
names(brm.tbc2) <- c("Estimate_TB.c","lCI","uCI")
brm.tbc2$Predictor <- c("Niche change (D trend)","Niche breadth trend","Niche distance trend")
quartz(height=3,width=4)
ggplot() + 
  geom_point(data=brm.tbc1, aes(x=Estimate_TB.c, y=Predictor))+
  geom_errorbar(data=brm.tbc1, aes(x=Estimate_TB.c, y=Predictor, xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05)) +
  geom_point(data=brm.tbc2, aes(x=Estimate_TB.c, y=Predictor), col="#F8766D") +
  geom_errorbar(data=brm.tbc2, aes(x=Estimate_TB.c, y=Predictor, xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05), col="#F8766D") +
  geom_vline(xintercept = 0)

brm.tbn1 <- summary(m25.brm)$fixed[4:6,c(1,3,4)]
names(brm.tbn1) <- c("Estimate_TN.c","lCI","uCI")
brm.tbn1$Predictor <- c("Niche change (D trend)","Niche breadth trend","Niche distance trend")
brm.tbn2 <- summary(m26.brm)$fixed[4:6,c(1,3,4)]
names(brm.tbn2) <- c("Estimate_TN.c","lCI","uCI")
brm.tbn2$Predictor <- c("Niche change (D trend)","Niche breadth trend","Niche distance trend")
quartz(height=3,width=4)
ggplot() + 
  geom_point(data=brm.tbn1, aes(x=Estimate_TN.c, y=Predictor))+
  geom_errorbar(data=brm.tbn1, aes(x=Estimate_TN.c, y=Predictor, xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05)) +
  geom_point(data=brm.tbn2, aes(x=Estimate_TN.c, y=Predictor), col="#F8766D") +
  geom_errorbar(data=brm.tbn2, aes(x=Estimate_TN.c, y=Predictor, xmin=lCI, xmax=uCI), width=.2, position=position_dodge(0.05), col="#F8766D") +
  geom_vline(xintercept = 0)










