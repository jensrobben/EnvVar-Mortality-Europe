##### 0) Preliminary settings  ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'abind',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'ecmwfr', 'httr',
              'stats', 'utils')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(sapply(packages, require, character.only = TRUE))

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source function
source('Functions/F. Diverse.R')
source('Functions/F. Population Weights.R')

# Country/countries of interest
ctry.spec <- c("ES")

# Minimum and maximum date
t_min <- as.Date('2013-01-01', format = '%Y-%m-%d')
t_max <- as.Date('2019-12-31', format = '%Y-%m-%d')

##### 1) Load shape file NUTS 3 European regions ----

# Read shape file - downloaded from Eurostat
shapefile <- read_sf('Shapefile Eurostat NUTS/NUTS_RG_20M_2021_3035.shp')

# Extract shapefile for NUTS 3 regions in the countries of interest 
nuts.spec <- 3
sf.ctry   <- shapefile[shapefile$CNTR_CODE %in% ctry.spec & 
                         shapefile$LEVL_CODE == nuts.spec,]

# Remove overseas areas
overseas  <- c('ES630', 'ES640', 'ES70')
ind.rm    <- unlist(sapply(overseas, function(x) 
  which(grepl(x, sf.ctry$NUTS_ID, fixed = TRUE))))
sf.ctry   <- sf.ctry[-c(ind.rm), ]

# Box long lat coordinates where country falls into
coords_ctry <- st_coordinates(st_transform(sf.ctry, 4326)) %>% 
  as.data.frame() %>% dplyr::select('X','Y')
long.range <- range(coords_ctry[,1])
lat.range  <- range(coords_ctry[,2])

##### 2) Download air pollution data (CAMS) ----

# Select air pollutants
var_list  <- c('no2', 'o3', 'pm10', 'pm2p5')
var_name  <- c('nitrogen_dioxide', 'ozone', 'particulate_matter_10um',
               'particulate_matter_2.5um')

# Login details - create account on CDS
wf_set_key(user    = "15882", 
           key     = "c75e9e6e-343f-4def-9b0a-7c419d4c1c48",
           service = "ads")

# Request information
request <- list(
  variable = "ozone",
  model = "ensemble",
  level = "0",
  type = "validated_reanalysis",
  year = "2013",
  month = c("01", "02", "03", "04", "05", "06"),
  format = "tgz",
  dataset_short_name = "cams-europe-air-quality-reanalyses",
  target = paste0('o3_2013_H1.tar.gz'))

# Download
wf_request(user = "15882", request = request, 
           path = paste0(getwd(),'/Data/CAMS'))

## Repeat the above downloading steps for every half year in 2013-2019 and for 
## every air pollutant

##### 3) Pre-process air quality data ----
read_cams <- function(tarfile, var){
  tmpdir  <- tempdir()
  tmpfile <- tempfile()
  dir.create(tmpfile)
  untar(tarfile, exdir = tmpfile)
  files <- list.files(tmpfile)

  nf <- length(files)
  nc_ds <- lapply(1:nf, function(s) 
    nc_open(filename = paste0(tmpfile, "\\", files[s])))
  
  # Extract longitude, latitude coordinates and time object
  long <- round(ncvar_get(nc_ds[[1]], "lon"),10) 
  lat  <- round(ncvar_get(nc_ds[[1]], "lat"),10)
  time <- lapply(1:nf, function(s) 
    ncvar_get(nc_ds[[s]], "time"))
  
  # Put as date object
  t_units <- lapply(1:nf, function(s) 
    ncatt_get(nc_ds[[s]], "time", "units"))
  t_ustr  <- lapply(1:nf, function(s) 
    strsplit(t_units[[s]]$value, " "))
  t_dstr  <- lapply(1:nf, function(s) 
    strsplit(unlist(t_ustr[[s]])[3], "-"))
  date    <- do.call('c', lapply(1:nf, function(s)
    ymd(t_dstr[[s]]) + dhours(time[[s]])))
  
  # Filter to save space
  # Indices to filter the spatial data-frames to the long/lat of interest
  id.long <- which(nc_ds[[1]]$dim$lon$vals >= long.range[1] & 
                     nc_ds[[1]]$dim$lon$vals <= long.range[2])
  id.lat  <- which(nc_ds[[1]]$dim$lat$vals >= lat.range[1] &
                     nc_ds[[1]]$dim$lat$vals <= lat.range[2])
  
  # Store the data in a 3-dimensional array (Kelvin - Celsius)
  arr.list <- lapply(1:nf, function(s) {
    obj     <- ncvar_get(nc_ds[[s]], var)
    obj[id.long,id.lat,]
  })
  arr.comb <- abind(arr.list, along = 3)
  gc()
  
  # Set dimnames
  dimnames(arr.comb) <- list(long[id.long], lat[id.lat], as.character(date))
  
  arr.comb
}

# Loop over different years
air.pol <- 'pm2p5'

for(y in year(t_min):year(t_max)){
  # Message
  message(paste0("Year ", y))
  
  ##### 1) Read air pollution files of year y ----
  
  # File names
  file0     <- paste0("Data/CAMS/")
  tfile1    <- paste0(file0, air.pol, '_', y, '_H1.tar.gz')
  tfile2    <- paste0(file0, air.pol, '_', y, '_H2.tar.gz')
  
  # Read files
  pol.yth1  <- read_cams(tarfile = tfile1, var = air.pol)
  gc()
  pol.yth2  <- read_cams(tarfile = tfile2, var = air.pol)
  gc()
  
  # Abind data
  pol.yt <- abind(pol.yth1, pol.yth2)
  rm(pol.yth1, pol.yth2)
  gc()
  
  # Extract longitude/latitude coordinates of grid
  long <- as.numeric(dimnames(pol.yt)[[1]])
  lat  <- as.numeric(dimnames(pol.yt)[[2]])
  time <- dimnames(pol.yt)[[3]]
  
  ##### 2) Intersection shape file with CAMS air quality interpolation grid ----
  
  # All grid points in CAMS air quality 
  dat <- expand.grid(long, lat)
  colnames(dat) <- c('Longitude', 'Latitude')
  
  # Put longitude latitude coordinates as sf object and on correct crs
  pnts_sf    <- st_as_sf(dat, coords = c('Longitude','Latitude'), crs = 4326) 
  pnts_trans <- pnts_sf                            
  tt1_trans  <- st_transform(sf.ctry, 4326)      
  
  # Make the intersection
  ints <- st_intersects(tt1_trans, pnts_trans, sparse = TRUE)
  for(k in 1:length(ints)){
    dat[ints[[k]],'region'] <- tt1_trans$NUTS_ID[k]
  }
  
  ##### 3) Add missing NUTS regions - closest grid point ----
  df <- dat %>% na.omit()
  
  ind.missing    <- which(! sf.ctry$NUTS_ID %in% unique(df$region))
  region.missing <- sf.ctry$NUTS_ID[ind.missing]
  
  if(length(ind.missing) > 0){ 
    centroid.region <- st_centroid(tt1_trans$geometry[ind.missing]) 
    PNT  <- st_coordinates(centroid.region) # Long-lat coordinates centroid
    GRID <- st_coordinates(pnts_trans) # Long-lat coordiantes grid points
    DIST <- lapply(1:nrow(PNT), function(s) sweep(GRID, 2, PNT[s,], `-`)) # Distances
    ind.min <- lapply(1:nrow(PNT), function(s) 
      which.min(sqrt(DIST[[s]][,1]^2 + DIST[[s]][,2]^2))) %>% unlist()
    
    df.add <- data.frame(st_coordinates(pnts_trans[ind.min,]), region.missing)
    colnames(df.add) <- c('Longitude', 'Latitude', 'region')
    df <- rbind(df, df.add)
  }
  
  ##### 4) Population weights ----
  
  # Weights - population weighted average
  wpop <- get_pop_weights(ctry.spec = ctry.spec, long.ext = long, lat.ext = lat)
  df   <- df %>% left_join(wpop, by = c('Longitude' = 'Long.ext', 
                                        'Latitude' = 'Lat.ext', 
                                        'region' = 'Region'))
  
  # If no point - weight of one
  if(length(ind.missing) > 0){ 
    index <- which(df$region %in% region.missing)
    df[index,-c(1:3)] <- 1
  }
  
  # If NA - zero population weight
  df[which(is.na(df), arr.ind = TRUE)] <- 0
  
  # Small check for population weights to sum up to 1 for each region
  test <- df %>% dplyr::group_by(region) %>%
    summarize(sum(WPC2000), sum(WPC2005), sum(WPC2010), sum(WPC2015), sum(WPC2020))
  
  # Split coordinates by region
  coord.region <- df %>% split(f = df$region)
  
  ##### 5) Average air pollutant over each region ----
  
  # Function to average over each region
  avg.AP.region <- function(df.sub){
    message(df.sub$region[1])
    lo.ind <- match(df.sub$Longitude, long)
    la.ind <- match(df.sub$Latitude, lat)
    obj <- sapply(1:nrow(df.sub), function(s) 
      pol.yt[lo.ind[s],la.ind[s],])
  
    objw <- obj %*% matrix(df.sub$WPC2015, ncol = 1) %>% as.vector()
    names(objw) <- dimnames(pol.yt)[[3]]
    objw
  }
  
  df.region.AP <- sapply(coord.region, function(sub) avg.AP.region(sub))
  
  ##### 6) Check if all times are present ----
  time1 <- ymd_hms(time)
  time2 <- seq(time1[1], tail(time1,1), by = 'hours')
  
  time.add <- time2[which(! time2 %in% time1)]
  df.region.AP.add <- matrix(NA, nrow = length(time.add), ncol = ncol(df.region.AP), 
                             dimnames = list(as.character(time.add), 
                                             colnames(df.region.AP)))
  df.region.AP <- rbind(df.region.AP, df.region.AP.add)
  df.region.AP <- df.region.AP[order(ymd_hms(rownames(df.region.AP))),]
  
  ##### 7) Daily average, minimum and maximum air pollutant level in each region ---- 
  mat.region.AVG.AP.daily <- sapply(1:ncol(df.region.AP), function(j)
    zoo::rollapply(unname(df.region.AP[,j]), 24, by = 24, mean))
  
  mat.region.MIN.AP.daily <- sapply(1:ncol(df.region.AP), function(j)
    zoo::rollapply(unname(df.region.AP[,j]), 24, by = 24, min))
  
  mat.region.MAX.AP.daily <- sapply(1:ncol(df.region.AP), function(j)
    zoo::rollapply(unname(df.region.AP[,j]), 24, by = 24, max))
  
  dates <- unique(as.Date(rownames(df.region.AP)))
  dimnames(mat.region.AVG.AP.daily) = dimnames(mat.region.MIN.AP.daily) = 
    dimnames(mat.region.MAX.AP.daily) <- list(as.character(dates), 
                                              colnames(df.region.AP))
  
  # Function to transform matrix to long data frame r (see ribiosUtils)
  matrix2longdf <- function(mat,
                            row.names, col.names,
                            longdf.colnames=c("row","column","value")) {
    if(missing(row.names)) row.names <- rownames(mat)
    if(missing(col.names)) col.names <- colnames(mat)
    
    if(is.null(row.names)) row.names <- 1:nrow(mat)
    if(is.null(col.names)) col.names <- 1:ncol(mat)
    
    value <- as.vector(mat)
    if(length(row.names)!=nrow(mat))
      warning("row.names is inconsistent with the matrix dim")
    if(length(col.names)!=ncol(mat))
      warning("col.names is inconsistent with the matrix dim")
    
    rn <- rep(row.names, ncol(mat))
    cn <- rep(col.names, each=nrow(mat))
    res <- data.frame(row=rn,
                      column=cn,
                      value=value)
    colnames(res) <- longdf.colnames
    return(res)
  }
  
  df.region.AVG.AP.daily  <- matrix2longdf(mat.region.AVG.AP.daily,
                                           longdf.colnames = c('Date', 'Region',
                                                               paste0(air.pol,'.avg.daily')))  
  df.region.MIN.AP.daily  <- matrix2longdf(mat.region.MIN.AP.daily,
                                           longdf.colnames = c('Date', 'Region',
                                                               paste0(air.pol,'.min.daily')))
  df.region.MAX.AP.daily  <- matrix2longdf(mat.region.MAX.AP.daily,
                                           longdf.colnames = c('Date', 'Region',
                                                               paste0(air.pol,'.max.daily')))
  
  df.region.AP.daily.temp <- df.region.AVG.AP.daily %>% 
    dplyr::full_join(df.region.MIN.AP.daily)  %>%  
    dplyr::full_join(df.region.MAX.AP.daily)
  
  df.region.AP.daily.temp$Date <- as.Date(df.region.AP.daily.temp$Date)
 
  # Rbind to previous result
  if(y > 2013){
    df.region.AP.daily <- rbind(df.region.AP.daily, df.region.AP.daily.temp) %>%
      arrange(Region, Date) } else if (y == 2013) {
        df.region.AP.daily <- df.region.AP.daily.temp
      }
  
  rm(pol.yt, df.region.AP.daily.temp, df.region.AVG.AP.daily, df, dat,
     mat.region.MIN.AP.daily, mat.region.AVG.AP.daily, mat.region.MAX.AP.daily,
     df.region.MAX.AP.daily, df.region.MIN.AP.daily, df.region.AP)
  gc()
}

file.save <- paste0("Data/CAMS/",paste0(ctry.spec, collapse = '_'),'_NUTS3_',
                    air.pol,'_daily.rds')
saveRDS(df.region.AP.daily, file = file.save)
