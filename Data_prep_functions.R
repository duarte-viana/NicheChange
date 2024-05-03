##############################################################################
# Climate and land use change leading to niche expansion and shifts in birds #
# Pablo M. Avidad, Miguel Clavero, Duarte S. Viana                           #
##############################################################################


# Helper functions
# From https://rpubs.com/clanescher/processingBBSdata with some tweaks

unzipBBS <- function(path) {
  
  print("Unzipping States.zip")
  unzip(paste(path, "States.zip",sep = "/"), overwrite = T)
  
  base <- paste(path, "States",sep = "/")
  
  filesZip <- list.files(base, pattern = ".zip")
  
  print("Unzipping files")
  pb <- txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
  
  for (i in 1:length(filesZip)) {
    #print(paste0("Unzipping file ", i))
    unzip(paste0(base, "/", filesZip[i]), 
          exdir = base,
          overwrite = T)
    setTxtProgressBar(pb,i)
    
  }
  
  unzip(paste(path, "Routes.zip",sep = "/"), overwrite = T)
  unzip(paste(path, "Weather.zip",sep = "/"), overwrite = T)
}

readBBS <- function(path) {
  base <- paste(path, "States",sep = "/")
  
  files <- list.files(base, pattern = ".csv")
  
  print("Reading in breeding birds...")
  pb <- txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
  
  data <- read.csv(paste0(base, "/", files[1]), colClasses = "character")
  colnames(data) <- tolower(colnames(data))
  setTxtProgressBar(pb,1)
  
  for (j in 2:length(files)) {
    tmp <- read.csv(paste0(base, "/", files[j]), colClasses = "character")
    colnames(tmp) <- tolower(colnames(tmp))
    setTxtProgressBar(pb,j)
    
    data <- rbind(data, tmp)
  }
  
  return(data)
}

summarizeBBS <- function(raw) {
  raw$count <- raw$speciestotal
  raw$id <- paste0(raw$statenum, raw$route)
  raw1 <- raw[,c("id", "year", "aou", "count")]
  
  return(raw1)
}

aouToSp <- function(raw, path) {
  spList <- read.fwf(file = paste(path, "/SpeciesList.txt",sep = "/"),
                     widths = c(7, 5, 50, 50, 50, 
                                50, 50, 50, 50),
                     skip = 11,
                     col.names = c("Seq", "aou", "English",
                                   "French", "Spanish",
                                   "Order", "Family", "Genus", "Species"),
                     strip.white = T, fileEncoding="latin1")
  
  spList$gs <- paste(spList$Genus, spList$Species, sep = " ")
  spList <- spList[,c("aou", "gs", "English", "Spanish", "Order", "Family", "Genus")]
  #spList$aou <- as.character(spList$aou)
  
  bbsSp <- left_join(raw, spList, by = "aou")
  bbsSp <- bbsSp[,c("id", "year", "aou", "gs", "English", "Spanish", "Order", "Family", "Genus", "count")]
  colnames(bbsSp) <- c("id", "year", "aou", "species", "English", "Spanish", "Order", "Family", "Genus", "count")
  
  return(bbsSp)
}

# adapted from Harris et al. 2018 (https://doi.org/10.7717/peerj.4278)
combine_subspecies = function(df){
  
  species_table <- df %>% 
    dplyr::select(aou, Spanish) %>% distinct()
  
  # Subspecies have two spaces separated by non-spaces
  subspecies_names <-  species_table %>%
    pull(Spanish) %>%
    grep(" [^ ]+ ", ., value = TRUE)
  
  subspecies_ids = species_table %>%
    filter(Spanish %in% subspecies_names) %>%
    pull(aou)
  
  # Drop all but the first two words to get the root species name,
  # then find the AOU code
  new_subspecies_ids <-  species_table %>%
    slice(match(word(subspecies_names, 1,2),
                species_table$Spanish)) %>%
    pull(aou)
  
  # replace the full subspecies names with species-level names
  for (i in seq_along(subspecies_ids)) {
    df$aou[df$aou == subspecies_ids[i]] = new_subspecies_ids[i]
  }
  
  df %>%
    group_by(id, year, aou, Latitude, Longitude) %>%
    filter(count == sum(count)) %>%
    ungroup()
}


#################################################################
# From Harris 2018
#From raw monthly climate values calculate all the bioclim variables.
#You should not call this directly to load bioclim vars. Instead call
#get_bioclim_data(), which will 1st try to load the data from the sqlite
#db before processing it all from scratch.
#Columns present must be c('year','month','value','clim_var','site_id')
###################################################################
process_bioclim_data=function(monthly_climate_data){
  
  #Offset the year by 6 months so that the window for calculating bioclim variables
  #will be July 1 - June 30. See https://github.com/weecology/bbs-forecasting/issues/114
  monthly_climate_data$year = with(monthly_climate_data, ifelse(month %in% 7:12, year+1, year))
  
  #Spread out the climate variables ppt, tmean, etc into columns
  monthly_climate_data = monthly_climate_data %>%
    spread(clim_var, value)
  
  #Process the quarter ones first.
  quarter_info=data.frame(month=1:12, quarter=c(3,3,3,4,4,4,1,1,1,2,2,2))
  bioclim_quarter_data= monthly_climate_data %>%
    left_join(quarter_info, by='month') %>%
    group_by(site_id, year, quarter) %>%
    summarise(precip=sum(ppt), temp=mean(tmean)) %>%
    ungroup() %>%
    group_by(site_id,year) %>%
    summarise(bio8=max_min_combo(temp, precip, max=TRUE),
              bio9=max_min_combo(temp, precip, max=FALSE),
              bio10=max(temp),
              bio11=min(temp),
              bio16=max(precip),
              bio17=min(precip),
              bio18=max_min_combo(precip, temp, max=TRUE),
              bio19=max_min_combo(precip, temp, max=FALSE)) %>%
    ungroup()
  
  #Next the yearly ones, joining the quartely ones  back in at the end.
  bioclim_data=monthly_climate_data %>%
    group_by(site_id, year) %>%
    mutate(monthly_temp_diff=tmax-tmin) %>%
    summarise(bio1=mean(tmean),
              bio2=mean(monthly_temp_diff),
              bio4=sd(tmean)*100,
              bio5=max_min_combo(tmax,tmean,max=TRUE),
              bio6=max_min_combo(tmin,tmean,max=FALSE),
              bio12=sum(ppt),
              bio13=max(ppt),
              bio14=min(ppt),
              bio15=cv(ppt)) %>%
    ungroup() %>%
    mutate(bio7=bio5-bio6,
           bio3=(bio2/bio7)*100) %>%
    full_join(bioclim_quarter_data, by=c('site_id','year'))
  
  return(bioclim_data)
}

###################################################################
#Helper function to calculate some of the bioclim variables,
#like "precip in coldest month"
###################################################################
max_min_combo=function(vec1,vec2,max=TRUE){
  #Return the value in vec1 in the position where
  #vec2 is either highest or lowest. But 1st check for na
  #values.
  if(any(is.na(vec1)) | any(is.na(vec2))){
    return(NA)
  } else  if(max){
    return(vec1[which.max(vec2)])
  } else {
    return(vec1[which.min(vec2)])
  }
}
