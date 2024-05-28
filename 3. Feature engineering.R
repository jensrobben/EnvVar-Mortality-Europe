##### 0) Preliminary settings ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('lubridate','dplyr','ggplot2', 'robust', 'eurostat', 'ISOweek',
              'grDevices', 'RColorBrewer','sf','stats', 'viridis')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(sapply(packages, require, character.only = TRUE))
`%notin%` <- Negate(`%in%`)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# All times and ISO week/year information
t_min <- as.Date('2013-01-01', format = '%Y-%m-%d')
t_max <- as.Date('2019-12-31', format = '%Y-%m-%d')

time_frame <- seq(t_min, t_max, by="days")
df_itime   <- data.frame('Date' = time_frame, 'ISOWeek' = isoweek(time_frame),
                         'ISOYear' = isoyear(time_frame))

# Countries of interest
ctry.spec <- c("ES")

# Colour scales
ablue  <- '#56B4E9'
agreen <- '#A3D596'
ared   <- '#F8766D'

##### 1) Load shape file NUTS 3 European regions ----

# Read shape file
shapefile <- read_sf('Shapefile Eurostat NUTS/NUTS_RG_20M_2021_3035.shp')

# Extract shapefile of NUTS levels for countries of interest 
nuts.spec <- 3
sf.ctry   <- shapefile[shapefile$LEVL_CODE == nuts.spec,]
sf.ctry   <- sf.ctry[which(sf.ctry$CNTR_CODE %in% ctry.spec),] 

# Remove overseas areas
overseas  <- c('ES70', 'ES630', 'ES640')
ind.rm    <- unlist(sapply(overseas, function(x) 
  which(grepl(x, sf.ctry$NUTS_ID, fixed = TRUE))))
sf.ctry   <- sf.ctry[-c(ind.rm), ]

# Final regions
regions <- sort(unique(sf.ctry$NUTS_ID))

# Plot NUTS 3 regions
plot(sf.ctry[,c('NUTS_ID','geometry')])

##### 2) Load weekly deaths and preprocess ----

# Download weekly death counts from Eurostat
d.xtw.all <- get_eurostat(id = 'demo_r_mweek3', cache = FALSE,
                          compress_file = FALSE, time_format = 'raw')
colnames(d.xtw.all) <- c('Freq','Unit','Sex','Age','Region','Time', 'Deaths')

# Weekly deaths - countries of interest
d.xtw.c <- d.xtw.all %>% dplyr::filter(Region %in% sf.ctry$NUTS_ID &
                                         Age %notin% c('UNK', 'TOTAL')) %>%
  dplyr::select(-c('Freq','Unit')) %>%
  mutate(ISOYear = as.integer(substr(Time, 1, 4)), 
         ISOWeek = as.integer(substr(Time, 7, 9)),
         Age = factor(Age, levels = c('Y_LT5', paste0('Y',seq(5,85,5), '-', 
                                                      seq(9,89,5)), 'Y_GE90')),
         .before = 1) %>%
  select(-c('Time')) %>% arrange(Sex, ISOYear, ISOWeek, Age) %>%
  dplyr::filter(!grepl(':', Deaths),
                as.integer(Age) >= 8, ISOYear <= year(t_max), ISOYear >= year(t_min))

d.xtw.c$Age    <- droplevels(d.xtw.c$Age)

# Aggregate deaths 65+ per date, region, sex
d.xtw.c <- d.xtw.c %>% dplyr::filter(Age %in% c('Y65-69', 'Y70-74', 'Y75-79', 
                                                'Y80-84', 'Y85-89', 'Y_GE90') & 
                                       Sex == 'F') %>%
  group_by(ISOYear, ISOWeek, Region, Sex) %>%
  summarize(Deaths = sum(Deaths)) %>% ungroup()

# Add date
wdate   <- ISOweek2date(sprintf("%d-W%02d-%d", d.xtw.c$ISOYear, d.xtw.c$ISOWeek, 1))
d.xtw.c <- d.xtw.c %>% mutate('Date' = wdate, .before = 1)

# Remove Sex column
df.d <- d.xtw.c %>% select(-c(Sex))

# Plot of death counts to show seasonality
region <- c('ES300','ES521', 'ES243')
dfsubr <- df.d %>% dplyr::filter(Region %in% region)
pvis <- ggplot(dfsubr) + 
  geom_line(aes(x = Date, y = Deaths, group = Region, col = Region),
            linewidth = 0.8) + 
  xlab('Date') + ylab('Death counts') + 
  theme_bw(base_size = 15) + 
  theme(legend.position = 'bottom') + 
  scale_color_manual(name = 'Region', 
                     values = c('ES243' = ared, 'ES300' = agreen,
                                'ES521' = ablue))
pvis

##### 3) Load daily weather data + weekly transformation: tg, tx, tn ----

# File average, maximum, and minimum temperature
file.tg  <- paste0("Data/E-OBS/", paste0(ctry.spec, collapse = '_'),
                   "_NUTS3_tg_daily.rds")
file.tx <- paste0("Data/E-OBS/", paste0(ctry.spec, collapse = '_'),
                  "_NUTS3_tx_daily.rds")
file.tn <- paste0("Data/E-OBS/", paste0(ctry.spec, collapse = '_'),
                  "_NUTS3_tn_daily.rds")

# Read temperature data files and join with ISO-date 
df.tn <- df_itime %>% left_join(readRDS(file = file.tn), by = 'Date')
df.tg <- df_itime %>% left_join(readRDS(file = file.tg), by = 'Date')
df.tx <- df_itime %>% left_join(readRDS(file = file.tx), by = 'Date')

df.t  <- df.tn %>% left_join(df.tg) %>% left_join(df.tx) %>%
  mutate(Day = yday(Date), .before = 'ISOWeek') 

# Other temperature related features 
excess_temp <- function(r){
  message(r)
  
  # Data from region r
  dfsub.daily <- df.t %>% dplyr::filter(Region == r)
  
  # Baseline fits
  fitr.tg <- robust::lmRob(tg ~ sin(2*pi*Day/365.25) + cos(2*pi*Day/365.25),
                 data = dfsub.daily)
  fitr.tx <- robust::lmRob(tx ~ sin(2*pi*Day/365.25) + cos(2*pi*Day/365.25), 
                 data = dfsub.daily)
  fitr.tn <- robust::lmRob(tn ~ sin(2*pi*Day/365.25) + cos(2*pi*Day/365.25), 
                 data = dfsub.daily)
  
  # Predictions
  pred.tg <- predict(fitr.tg, dfsub.daily)
  pred.tx <- predict(fitr.tx, dfsub.daily)
  pred.tn <- predict(fitr.tn, dfsub.daily)
  
  # Quantiles
  qtg <- quantile(dfsub.daily$tg, c(0.05, 0.95), na.rm = TRUE)
  qtx <- quantile(dfsub.daily$tx, c(0.05, 0.95), na.rm = TRUE)
  qtn <- quantile(dfsub.daily$tn, c(0.05, 0.95), na.rm = TRUE)
  
  # Information
  dd <- data.frame(Date = dfsub.daily$Date, Region = dfsub.daily$Region,
                   tg_anom = dfsub.daily$tg - pred.tg,
                   tx_anom = dfsub.daily$tx - pred.tx, 
                   tn_anom = dfsub.daily$tn - pred.tn,
                   Tind95 = (dfsub.daily$tg > qtg[2])*1 + 
                     (dfsub.daily$tx > qtx[2])*1 + (dfsub.daily$tn > qtn[2])*1,
                   Tind5 = (dfsub.daily$tg < qtg[1])*1 + 
                     (dfsub.daily$tx < qtx[1])*1 + (dfsub.daily$tn < qtn[1])*1)
  dd
}

list.exc.t <- lapply(intersect(sort(unique(df.t$Region)),regions), excess_temp)
df.exc.t   <- do.call('rbind', list.exc.t)
df.t       <- df.t %>% left_join(df.exc.t)

# Transform to weekly - on ISO week basis
df.t.w <- df.t %>% group_by(ISOYear, ISOWeek, Region) %>%
  summarize(w_avg_tg  = mean(tg), w_avg_tx  = max(tx), w_avg_tn  = min(tn),
            w_avg_tg_anom = mean(tg_anom), w_avg_tx_anom = mean(tx_anom), 
            w_avg_tn_anom = mean(tn_anom), w_avg_Tind95 = mean(Tind95), 
            w_avg_Tind5 = mean(Tind5)) %>%
  ungroup() 

# Plot on world map
isoy <- 2018
isow <- 30

temp.wrldmp <- df.t.w %>% 
  dplyr::filter(ISOYear == isoy & ISOWeek == isow) %>%
  full_join(sf.ctry, by = c('Region' = 'NUTS_ID'))

# Weekly average of maximum temperature anomalies
pebtx <- ggplot() + 
  theme_bw(base_size = 10) + 
  geom_sf(data = temp.wrldmp$geometry, aes(fill = temp.wrldmp$w_avg_tx_anom,
                                           group = temp.wrldmp$CNTR_CODE,
                                           col = temp.wrldmp$w_avg_tx_anom)) + 
  scale_fill_viridis(name = 'Tmax', option = 'rocket', direction = -1,
                     na.value = 'white') + 
  scale_color_viridis(name = 'Tmax', option = 'rocket', direction = -1,
                      na.value = 'white') + 
  ggtitle(paste0('w_avg_Tmax_anom: ', isoy, '-',isow)) + 
  theme(legend.position = 'bottom',
        legend.title = element_text(margin = margin(t = -1.5, unit = "lines")))

pTind <- ggplot() + 
  theme_bw(base_size = 10) + 
  geom_sf(data = temp.wrldmp$geometry, aes(fill = temp.wrldmp$w_avg_Tind95, 
                                           group = temp.wrldmp$CNTR_CODE,
                                           col = temp.wrldmp$w_avg_Tind95)) + 
  scale_fill_viridis(name = 'Tind95', option = 'rocket', direction = -1,
                     na.value = 'white') +
  scale_color_viridis(name = 'Tind95', option = 'rocket', direction = -1,
                      na.value = 'white') + 
  ggtitle(paste0('w_avg_Tind95: ', isoy, '-',isow)) +
  theme(legend.position = 'bottom',
        legend.title = element_text(margin = margin(t = -1.5, unit = "lines")))

##### 4) Load daily weather data + weekly transformation: hum, wind sp, rain fall ----

# File names humidity, rain fall, wind speed
file.hu  <- paste0("Data/E-OBS/", paste0(ctry.spec, collapse = '_'), 
                   "_NUTS3_hu_daily.rds")
file.rr <- paste0("Data/E-OBS/", paste0(ctry.spec, collapse = '_'),
                  "_NUTS3_rr_daily.rds")
file.fg <- paste0("Data/E-OBS/", paste0(ctry.spec, collapse = '_'), 
                  "_NUTS3_fg_daily.rds")

# Read climate data files and join with isodate information
df.hu = df.rr = df.fg <- df_itime

df.hu <- df.hu %>% left_join(readRDS(file = file.hu), by = 'Date')
df.rr <- df.rr %>% left_join(readRDS(file = file.rr), by = 'Date')
df.fg <- df.fg %>% left_join(readRDS(file = file.fg), by = 'Date')

df.hrf  <- df.hu %>% left_join(df.rr) %>% left_join(df.fg) %>%
  mutate(Day = yday(Date), .before = ISOWeek) 

# Other climate related features of hu, fg, rr 

excess_cli <- function(r){
  message(r)
  
  # Data from region r
  dfsub.daily <- df.hrf %>% dplyr::filter(Region == r)
  
  # Baseline fits
  fitr.hu <- robust::lmRob(hu ~ sin(2*pi*Day/365.25) + cos(2*pi*Day/365.25),
                       data = dfsub.daily, control = robust::lmRob.control(mxr = 500))
  fitr.fg <- robust::lmRob(fg ~ sin(2*pi*Day/365.25) + cos(2*pi*Day/365.25),
                          data = dfsub.daily, control = robust::lmRob.control(mxr = 500))
  
  # Predictions
  pred.hu <- predict(fitr.hu, dfsub.daily)
  pred.fg <- predict(fitr.fg, dfsub.daily)

  # Quantiles
  qhu <- quantile(dfsub.daily$hu, c(0.05, 0.95), na.rm = TRUE)
  qfg <- quantile(dfsub.daily$fg, c(0.05, 0.95), na.rm = TRUE)
  qrr <- quantile(dfsub.daily$rr, c(0.05, 0.95), na.rm = TRUE)
  
  # Data frame
  data.frame(Date = dfsub.daily$Date, Region = dfsub.daily$Region,
             hu_anom  = dfsub.daily$hu - pred.hu,
             rr_anom  = dfsub.daily$rr - 0,
             fg_anom  = dfsub.daily$fg - pred.fg,
             hu.ind95 = (dfsub.daily$hu >= qhu[2])*1,
             hu.ind5  = (dfsub.daily$hu <= qhu[1])*1,
             fg.ind95 = (dfsub.daily$fg >= qfg[2])*1,
             fg.ind5  = (dfsub.daily$fg <= qfg[1])*1,
             rr.ind95 = (dfsub.daily$rr >= qrr[2])*1,
             rr.ind5  = (dfsub.daily$rr <= qrr[1])*1)
}

list.exc.cli <- lapply(intersect(sort(unique(df.hrf$Region)),regions), excess_cli)
df.exc.cli   <- do.call('rbind', list.exc.cli)

df.hrf <- df.hrf %>% left_join(df.exc.cli)

# Transform to weekly - on ISO week basis
df.hrf.w <- df.hrf %>% group_by(ISOYear, ISOWeek, Region) %>%
  summarize(w_avg_hu = mean(hu), w_avg_rr = mean(rr), w_avg_fg = mean(fg),
            w_avg_hu_anom  = mean(hu_anom), w_avg_rr_anom = mean(rr_anom),
            w_avg_fg_anom  = mean(fg_anom), 
            w_avg_hu.ind95 = mean(hu.ind95), w_avg_hu.ind5 = mean(hu.ind5),
            w_avg_rr.ind95 = mean(rr.ind95), w_avg_rr.ind5 = mean(rr.ind5),
            w_avg_fg.ind95 = mean(fg.ind95), w_avg_fg.ind5 = mean(fg.ind5)) %>%
  ungroup() 


# Plot weekly average of daily relative humidity anomalies
isoy <- 2018
isow <- 30

hrf.wrldmp <- df.hrf.w %>% 
  dplyr::filter(ISOYear == isoy & ISOWeek == isow) %>%
  full_join(sf.ctry, by = c('Region' = 'NUTS_ID'))

pal <- colorRampPalette(brewer.pal(9,"Blues"))(50)

phum <- ggplot() + 
  theme_bw(base_size = 10) + 
  geom_sf(data = hrf.wrldmp$geometry, aes(fill = hrf.wrldmp$w_avg_hu_anom, 
                                          group = hrf.wrldmp$CNTR_CODE,
                                          col = hrf.wrldmp$w_avg_hu_anom)) + 
  scale_fill_gradientn(colours = pal, name = 'Hum', na.value = 'white') + 
  scale_color_gradientn(colours = pal, name = 'Hum', na.value = 'white') +   
  ggtitle(paste0('w_avg_Hum_anom: ', isoy, '-',isow)) + 
  theme(legend.position = 'bottom',
        legend.title = element_text(margin = margin(t = -1.5, unit = "lines"))) 
  
phum



##### 5) Load daily air pollution data + weekly transformation ----

# Merge daily data
airpols <- c('o3','pm10', 'pm2p5', 'no2')
file.airpols <- paste0('Data/CAMS/', rep(paste0(ctry.spec, collapse = '_'), 
                                         length(airpols)),'_NUTS3_',
                       airpols, '_daily.rds')
list.ap <- rep(list(df_itime), length(airpols))
names(list.ap) <- airpols

for(a in 1:length(airpols)){
    list.ap[[a]] <- list.ap[[a]] %>% 
      left_join(readRDS(file.airpols[a]), by = 'Date')
}

df.ap <- list.ap[[1]]
for(j in 2:length(airpols)){
  df.ap <- df.ap %>% left_join(list.ap[[j]])
}

df.ap$Day <- yday(df.ap$Date)

# Other air pollution related features 
excess_ap <- function(r){
  message(r)
  
  # Data from region r
  dfsub.daily     <- df.ap %>% dplyr::filter(Region == r)
  
  # Fits
  fits <- list()
  for(a in airpols){
    form <- paste0(a, ".avg.daily ~ sin(2*pi*Day/365.25) + cos(2*pi*Day/365.25)")
    fits[[which(a ==  airpols)]] <- robust::lmRob(as.formula(form), 
                                                  data = dfsub.daily,
                                                  control = robust::lmRob.control(mxr = 500))
  }
  names(fits) <- airpols
  
  # Put in data frame 
  df <- data.frame(Date = dfsub.daily$Date, Region = dfsub.daily$Region)
  
  for(a in airpols){
    df[,paste0(a,'_anom')] <- dfsub.daily[,paste0(a,'.avg.daily')] - 
      fits[[a]]$fitted.values
  }
  
  # Quantiles 
  qo3    <- quantile(dfsub.daily$o3.avg.daily, c(0.05, 0.95), na.rm = TRUE)
  qno2   <- quantile(dfsub.daily$no2.avg.daily, c(0.05, 0.95), na.rm = TRUE)
  qpm10  <- quantile(dfsub.daily$pm10.avg.daily, c(0.05, 0.95), na.rm = TRUE)
  qpm2p5 <- quantile(dfsub.daily$pm2p5.avg.daily, c(0.05, 0.95), na.rm = TRUE)
  
  for(a in airpols){
    df[,paste0(a,'.ind95')] <- (dfsub.daily[,paste0(a,'.avg.daily')] >= 
                                    get(paste0('q',a))[2])*1
    df[,paste0(a,'.ind5')] <- (dfsub.daily[,paste0(a,'.avg.daily')] <= 
                                    get(paste0('q',a))[1])*1
  }
  
  df
}

list.exc.ap <- lapply(regions, excess_ap)
df.exc.ap   <- do.call('rbind', list.exc.ap)
df.ap       <- df.ap %>% left_join(as.data.frame(df.exc.ap))

# Transform to weekly - on ISO week basis
list.ap.w <- list()
for(j in 1:length(airpols)){
  list.ap.w[[j]] <- df.ap %>% group_by(ISOYear, ISOWeek, Region) %>%
    summarize(!!paste0('w_avg_',airpols[j]) := mean(!!sym(paste0(airpols[j],'.avg.daily'))),
              !!paste0('w_avg_',airpols[j],'_anom')  := mean(!!sym(paste0(airpols[j],'_anom'))),
              !!paste0('w_avg_',airpols[j],'.ind95') := mean(!!sym(paste0(airpols[j],'.ind95'))),
              !!paste0('w_avg_',airpols[j],'.ind5')  := mean(!!sym(paste0(airpols[j],'.ind5')))) %>%
    ungroup() 
}

df.ap.w <- list.ap.w[[1]]
for(j in 2:length(airpols)){
  df.ap.w <- df.ap.w %>% left_join(list.ap.w[[j]])
}

# Plot on world map
isoy <- 2018
isow <- 30

ap.wrldmp <- df.ap.w %>% 
  dplyr::filter(ISOYear == isoy & ISOWeek == isow) %>%
  full_join(sf.ctry, by = c('Region' = 'NUTS_ID'))

pal <- colorRampPalette(brewer.pal(9,"Greys"))(100)

po3 <- ggplot() + 
  theme_bw(base_size = 10) + 
  geom_sf(data = ap.wrldmp$geometry, aes(fill = ap.wrldmp$w_avg_o3_anom, 
                                         group = ap.wrldmp$CNTR_CODE,
                                         col = ap.wrldmp$w_avg_o3_anom)) +
  scale_fill_gradientn(colours = pal, name = 'O3', na.value = 'white') +
  scale_colour_gradientn(colours = pal, name = 'O3', na.value = 'white') + 
  ggtitle(paste0("w_avg_O3_anom: ", isoy, '-',isow)) + 
  theme(legend.position = 'bottom',
        legend.title = element_text(margin = margin(t = -1.5, unit = "lines")))


##### 6) Combine all data sets ----

# Combine
df <- df.d %>% 
  left_join(df.t.w,   by = c('ISOYear', 'ISOWeek', 'Region')) %>%
  left_join(df.hrf.w, by = c('ISOYear', 'ISOWeek', 'Region')) %>%
  left_join(df.ap.w,  by = c('ISOYear', 'ISOWeek', 'Region')) %>%
  na.omit()

# Save
file <- paste0('Results/Data/df_NUTS', nuts.spec, '_',
               paste0(ctry.spec, collapse = '-'), '.rds')

saveRDS(df, file)
 