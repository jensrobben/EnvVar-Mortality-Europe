##### 0) Preliminary settings ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('dplyr','sf','ncdf4','lubridate','tidyverse','tibble', 'spdep',
              'Polychrome','tidyr','ISOweek', 'ggplot2', 'data.table', 'xgboost',
              'parallel', 'doParallel')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(sapply(packages, require, character.only = TRUE))

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
ablue   <- '#56B4E9'
agreen  <- '#A3D596'
ared    <- '#F8766D'
aorange <- '#FFB347'
apurple <- '#9E5B9B'
ayellow <- '#FFD92F'
apink   <- '#FF6EB4'
ateal   <- '#008080'

# Source functions
source('Functions/F. Diverse.R')

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

##### 2) Data processing ----

# Training data
file    <- paste0('Results/Data/df_NUTS', nuts.spec, '_',
                  paste0(ctry.spec, collapse = '-'), '.rds')
df      <- readRDS(file)

# Add baseline deaths 
fitb <- readRDS(paste0('Results/Data/bDeaths_NUTS', nuts.spec, '_',
                       paste0(ctry.spec, collapse = '-'), '.rds'))
df   <- df %>% left_join(fitb[,c('Date','Region','bDeaths')], 
                         by = c('Date','Region')) %>% na.omit()

# Season
season    <- getSeason(df$Date)
df$Summer <- (season == 'Summer')*1
df$Winter <- (season == 'Winter')*1
df$Spring <- (season == 'Spring')*1
df$Fall   <- (season == 'Fall')*1

# Long/Lat coordinates
ctrd <- st_centroid(st_transform(sf.ctry$geometry, 4326)) %>% st_coordinates()
colnames(ctrd) <- c('long', 'lat')
sf.ctry <- cbind(sf.ctry, ctrd)
df <- df %>% left_join(sf.ctry %>% as.data.frame() %>% select(c('long','lat','NUTS_ID')), 
                       by = c('Region' = 'NUTS_ID'))

# Create lagged variables
var.t.lag <- c('w_avg_tx_anom', 'w_avg_tn_anom', 'w_avg_Tind95', 'w_avg_Tind5', 
               'w_avg_hu_anom', 'w_avg_rr_anom', 'w_avg_fg_anom', 
               'w_avg_hu.ind95', 'w_avg_hu.ind5', 'w_avg_rr.ind95', 'w_avg_rr.ind5',
               'w_avg_fg.ind95', 'w_avg_fg.ind5', 'w_avg_o3_anom', 'w_avg_o3.ind95',
               'w_avg_o3.ind5', 'w_avg_pm10_anom', 'w_avg_pm10.ind95',
               'w_avg_pm10.ind5', 'w_avg_pm2p5_anom', 'w_avg_pm2p5.ind95',
               'w_avg_pm2p5.ind5', 'w_avg_no2_anom', 'w_avg_no2.ind95',
               'w_avg_no2.ind5')
list.lagdf <- list()
for (v in var.t.lag){
  list.lagdf[[which(var.t.lag == v)]] <- df %>% group_by(Region) %>% 
    reframe(Date, !!paste0(v,'_l1') := lag(!!sym(v), 1))
}
df.lag <- plyr::join_all(list.lagdf, by = c('Region', 'Date'))
df     <- df %>% left_join(df.lag, by = c('Region', 'Date'))

# Input features XGBoost
vars.nlag <- var.t.lag
vars.lag  <- paste0(vars.nlag, '_l1')
vars.loc  <- c('long','lat')
vars.time <- c('Summer', 'Winter', 'Spring', 'Fall')
vars      <- c(vars.nlag, vars.lag, vars.loc, vars.time)

# Remove NA's (first week for each region due to lagged features)
df <- df %>% na.omit()

# Select features in df
df <- df %>% select(c('Date', 'ISOYear', 'ISOWeek', 'Deaths', 'bDeaths', 
                      'Region', all_of(vars))) 
df$Date <- as.Date(df$Date)

##### 3) Tuning the machine learning model -----

# Split the dataset into 7 different folds according to the year
df <- df %>% mutate(Fold = (ISOYear - year(t_min)) + 1)

# Save final training data set
file    <- paste0('Results/Data/df_tune_NUTS3_ES.rds')
saveRDS(df, file)

# Number of cross validations
ncv <- length(unique(df$Fold))

# Poisson deviance
poisson_deviance <- function(data, lev = NULL, model = NULL) {
  ypred <- data$pred
  yobs  <- data$obs
  out   <- 2 * sum(yobs * log((yobs + 10^(-50))/ypred) - (yobs - ypred))
  names(out) <- 'DEV'
  out
}

# Custom objective function - Poisson
myobjective <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  offset <- attr(dtrain, "offset")
  grad <- (exp(preds + offset) - labels)
  hess <- exp(preds + offset)
  return(list(grad = grad, hess = hess))
}

# Custom Metric - Poisson
evalerror <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  offset <- attr(dtrain, "offset")
  err    <- mean(exp(preds + offset) - labels*(preds + offset))
  return(list(metric = "MyError", value = err))
}

# Tuning grid - add fold to grid (can be extended)
grid <- expand.grid(nrounds = c(seq(10, 1000, 10)),
                    max_depth = c(1,3,5,7,9),
                    fold = 1:max(df$Fold))

# Reducing tuning grid
gridr <- unique(grid[,c('max_depth', 'fold')])

# Function to fit XGBoost on specific training fold and tuning parameter combination
cv.xgb <- function(df.train, df.val, vars, nrounds.max, nrounds.vec, max_depth){
  
  # XGBoost Matrix
  dtrain <- xgb.DMatrix(data  = as.matrix(df.train %>% select(vars)),
                        label = as.matrix(df.train %>% select(Deaths)))
  dval   <- xgb.DMatrix(data  = as.matrix(df.val %>% select(vars)),
                        label = as.matrix(df.val %>% select(Deaths)))
  
  vfold  <- df.val$Fold[1]
  attr(dtrain, 'offset') <- log(df.train[['bDeaths']])
  attr(dval, 'offset')   <- log(df.val[['bDeaths']])
  
  # Watchlist
  watchlist <- list(eval1 = dtrain, eval2 = dval)
  
  set.seed(1996)
  xgb1 <- xgb.train(params = list(booster = 'gbtree',
                                  eta = 0.01,
                                  max_depth = max_depth,
                                  min_child_weight = 100,
                                  subsample = 0.75,
                                  colsample_bytree = 0.75,
                                  base_score = 0,
                                  lambda = 0,
                                  objective = myobjective,
                                  eval_metric = evalerror),
                    watchlist = watchlist,
                    early_stopping_rounds = NULL,
                    data = dtrain,
                    maximize = FALSE,
                    nrounds = nrounds.max,
                    verbose = 1,
                    nthread = 1)
  
  # Predictions validation set
  pred <- lapply(nrounds.vec, function(ntr) 
    exp(predict(object = xgb1, dval, iterationrange = c(1,ntr)) + 
          attr(dval, 'offset')))
  
  # Poisson deviance
  pois.dev.base <- poisson_deviance(data.frame('obs' = df.val$Deaths, 
                                               'pred' = df.val$bDeaths))
  pois.dev.gbm  <- sapply(1:length(nrounds.vec), function(x)
    poisson_deviance(data.frame('obs' = df.val$Deaths, 
                                'pred' = unlist(pred[[x]]))))
  
  names(pois.dev.gbm) <- nrounds.vec
  
  pois.dev.gbm
}

# Tuning by cross-validation
nC <- 12
cl <- snow::makeCluster(nC)
doSNOW::registerDoSNOW(cl)
pb <- txtProgressBar(max = nrow(gridr), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
results <- foreach(k = 1:nrow(gridr),
                   .packages = c('dplyr', 'xgboost'), .options.snow = opts,
                   .noexport = c('df', 'file', 'wd_data')) %dopar% {
                     # Load data set
                     wd_data <- "C:/Users/u0131219/OneDrive - KU Leuven/Gridded datasets/"
                     file    <- paste0(wd_data, 'Weekly combined/df_tune.rds')
                     df      <- readRDS(file)
                     
                     # Number of iterations
                     nrounds.max <- max(grid$nrounds)
                     nrounds.vec <- unique(grid$nrounds)
                     
                     # Parameter values
                     max_depth        <- gridr[k,'max_depth']
                     vfold            <- gridr[k,'fold']
                     
                     # Training and validation set
                     df.train <- df %>% dplyr::filter(Fold != vfold)
                     df.val   <- df %>% dplyr::filter(Fold == vfold)
                     
                     # Fit GBM for each parameter combination in tuning grid
                     obj <- cv.xgb(df.train, df.val, vars, nrounds.max,
                                   nrounds.vec, max_depth)
                     
                     obj
                   }
close(pb)
parallel::stopCluster(cl)

# Set names 
names(results) <- sapply(1:nrow(gridr), function(j) paste0('maxdepth_', gridr[j,'max_depth'],
                                                           '/vfold_',gridr[j,'fold']))
names(results) 

# Save the results
file.xgb <-  'Results/XGB/xgbCV_results.rds'
# saveRDS(results, file.xgb)

# Read the results
results.list <- readRDS(file = file.xgb)
results      <- do.call('rbind', results.list) %>% as.data.frame()

# Manipulate data frame
info              <- str_extract_all(rownames(results),"\\(?[0-9,.]+\\)?")
results$max_depth <- as.numeric(unlist(sapply(info,  `[`, 1)))
results$vfold     <- as.numeric(unlist(sapply(info,  `[`, 2)))
rownames(results) <- NULL

# Sum accross folds
results <- results %>% 
  group_by(max_depth) %>%
  summarise(across(everything(), sum)) %>%
  dplyr::select(-c('vfold')) %>% ungroup()

# Find minimum
trees.vec <- as.numeric(colnames(results)[-c(1)])
min.val   <- apply(results[-c(1)], 1, min)
ind.min   <- apply(results[-c(1)], 1, which.min)
results   <- results[,c(1)] %>%
  mutate('ntree_opt' = trees.vec[ind.min],
         'minCV' = min.val)
tune.opt <- results[which.min(results$minCV),]

###### Fit xgb on entire training data set with optimal parameters

# Full data
file    <- 'Results/Data/df_tune_NUTS3_ES.rds'
df      <- readRDS(file)

# XGBoost Matrix
xgb.df <- xgb.DMatrix(data  = as.matrix(df %>% select(all_of(vars))),
                      label = as.matrix(df %>% select(Deaths)))

attr(xgb.df, 'offset') <- log(df[['bDeaths']])

# Watchlist
watchlist <- list(eval1 = xgb.df)

# Fit
set.seed(1996)
xgb.opt <- xgb.train(params = list(booster = 'gbtree',
                                   eta = 0.01,
                                   max_depth = tune.opt$max_depth,
                                   min_child_weight = 100,
                                   subsample = 0.75,
                                   colsample_bytree = 0.75,
                                   base_score = 0,
                                   lambda = 0,
                                   objective = myobjective,
                                   eval_metric = evalerror),
                  watchlist = watchlist,
                  early_stopping_rounds = NULL,
                  data = xgb.df,
                  maximize = FALSE,
                  nrounds = tune.opt$ntree_opt,
                  verbose = 1)
# Save
filexgb <- 'Results/XGB/xgb_fit.rds'
# saveRDS(xgb.opt, file = filexgb)

