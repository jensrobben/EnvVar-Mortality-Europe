get_pop_weights <- function(ctry.spec, long.ext, lat.ext){
  
  ##### 1) Read Population gridded dataset ----
  
  # Open NetCDF Files 
  nc_ds <- nc_open(filename = "Data/SEDAC/gpw_v4_population_count_rev11_2pt5_min.nc")
  
  # Extract longitude, latitude coordinates 
  long <- ncvar_get(nc_ds, "longitude") 
  lat  <- ncvar_get(nc_ds, "latitude")

  # Store the data in a 3-dimensional array 
  pc.array <- ncvar_get(nc_ds) 
  dimnames(pc.array) <- list(long, lat, 1:20)
  
  ##### 2) Load shape file NUTS2/NUTS3 European countries (Eurostat) ----
  
  # Read shape file
  shapefile <- read_sf('Shapefile Eurostat NUTS/NUTS_RG_20M_2021_3035.shp')
  
  # Extract shapefile of NUTS levels for countries of interest 
  nuts.spec <- 3
  sf.ctry   <- shapefile[shapefile$CNTR_CODE %in% ctry.spec & 
                           shapefile$LEVL_CODE == nuts.spec,]
  
  # Remove overseas areas
  overseas  <- c('ES70', 'ES630', 'ES640')
  ind.rm    <- unlist(sapply(overseas, function(x) 
    which(grepl(x, sf.ctry$NUTS_ID, fixed = TRUE))))
  sf.ctry   <- sf.ctry[-c(ind.rm), ]
  
  ##### 3) Intersection shape file grid with population grid ----
  
  # All grid points 
  dat.full <- expand.grid(long, lat)
  colnames(dat.full) <- c('Longitude', 'Latitude')
  
  # Filter
  coords_ctry <- st_coordinates(st_transform(sf.ctry, 4326)) %>% 
    as.data.frame() %>% dplyr::select('X','Y')
  long.range <- range(coords_ctry[,1])
  lat.range  <- range(coords_ctry[,2])
  
  dat <- dat.full %>% dplyr::filter(Longitude >= long.range[1],
                                    Longitude <= long.range[2],
                                    Latitude  >= lat.range[1],
                                    Latitude  <= lat.range[2])
  
  # Filter extra for computational reasons
  dt_pnts  <- st_as_sf(dat, coords = c('Longitude','Latitude'), crs = 4326) 
  keep     <- st_intersects(dt_pnts, st_transform(st_union(sf.ctry), 4326))
  keep.ind <- which(! is.na(as.integer(keep)))
  dat      <- dat[keep.ind,]
  
  # Put on same CRS
  pop_sf  <- st_as_sf(dat, coords = c('Longitude','Latitude'), crs = 4326) 
  nuts_sf <- st_transform(sf.ctry, 4326)      
  
  # Make the intersection
  ints <- st_intersects(nuts_sf, pop_sf, sparse = TRUE)
  for(k in 1:length(ints)){
    dat[ints[[k]],'Region'] <- nuts_sf$NUTS_ID[k]
  }

  ##### 4) Population count per Region ----
  df <- dat %>% na.omit()
  
  # Split coordinates by Region
  coord.region <- df %>% split(f = df$Region)
  
  # Function to aggregate over each Region
  tot.pop.region <- function(df.sub){
    coords <- df.sub[,c(1,2)]
    
    obj <- sapply(1:nrow(coords), function(s) 
      pc.array[as.character(coords[s,1]), as.character(coords[s,2]),])
    
    rowSums(obj, na.rm = TRUE)[1:5]
  }
  
  df.region.Pop <- sapply(coord.region, function(sub) tot.pop.region(sub))
  rownames(df.region.Pop) <- c('TPC2000', 'TPC2005', 
                               'TPC2010', 'TPC2015', 
                               'TPC2020')

  df <- df %>% left_join(df.region.Pop %>% t() %>% as.data.frame() %>% 
                           rownames_to_column(var = 'Region'), by = c('Region'))

  ##### 5) Intersection with given grid ---- 
  dat.ext           <- expand.grid(long.ext, lat.ext)
  colnames(dat.ext) <- c('Longitude', 'Latitude')
  
  ext_sf <- st_as_sf(dat.ext, coords = c('Longitude','Latitude'), crs = 4326) 
  ints2  <- st_intersects(nuts_sf, ext_sf, sparse = TRUE)
  for(k in 1:length(ints2)){
    dat.ext[ints2[[k]],'Region'] <- nuts_sf$NUTS_ID[k]
  }
  df.ext <- dat.ext %>% na.omit()

  # closest point, in same region
  regions <- sort(sf.ctry$NUTS_ID)
  vars    <- c('Longitude', 'Latitude')
  sdf     <- df %>% split(f = df$Region)
  sdf.ext <- df.ext %>% split(f = df.ext$Region)
  regions <- unique(df.ext$Region)
  dst <- lapply(regions, function(r) {
    min.dist <- apply(proxy::dist(sdf.ext[[r]][,vars], sdf[[r]][,vars]), 2, which.min)
    ind      <- as.numeric(names(min.dist))
    cbind(sdf.ext[[r]][min.dist, c('Longitude','Latitude','Region')],
          'Ind' = ind)
  })
  newdf   <- do.call('rbind', dst)
  rownames(newdf) <- NULL
  
  df[,c('Long.ext', 'Lat.ext','Region.ext')] <- NA
  df[newdf$Ind,c('Long.ext', 'Lat.ext','Region.ext')] <- 
    newdf[,c('Longitude','Latitude','Region')]
  df$Ind <- as.numeric(rownames(df))

  df <- df %>% select(-c('Region.ext'))
  
  i1 <- match(df$Longitude, long)
  i2 <- match(df$Latitude, lat)
  df[,'PC2000'] <- pc.array[cbind(i1,i2,1)]
  df[,'PC2005'] <- pc.array[cbind(i1,i2,2)]
  df[,'PC2010'] <- pc.array[cbind(i1,i2,3)]
  df[,'PC2015'] <- pc.array[cbind(i1,i2,4)]
  df[,'PC2020'] <- pc.array[cbind(i1,i2,5)]
  
  # Aggregrate over extern grid
  dff <- df %>% group_by(Long.ext, Lat.ext) %>%
    na.omit() %>% summarize('Region'  = unique(Region),
                            'WPC2000' = sum(PC2000, na.rm = TRUE)/unique(TPC2000),
                            'WPC2005' = sum(PC2005, na.rm = TRUE)/unique(TPC2005),
                            'WPC2010' = sum(PC2010, na.rm = TRUE)/unique(TPC2010),
                            'WPC2015' = sum(PC2015, na.rm = TRUE)/unique(TPC2015),
                            'WPC2020' = sum(PC2020, na.rm = TRUE)/unique(TPC2020)) %>%
    arrange(Region, Long.ext, Lat.ext)
  
  dff
}
