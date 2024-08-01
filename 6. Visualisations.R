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

# Colour scales
ablue   <- '#56B4E9'
agreen  <- '#A3D596'
ared    <- '#F8766D'
aorange <- '#FFB347'
apurple <- '#9E5B9B'
ayellow <- '#FFD92F'
apink   <- '#FF6EB4'
ateal   <- '#008080'

##### 1) XGBoost predictions ----

# Full data
file    <- 'Results/Data/df_tune_NUTS3_ES.rds'
df      <- readRDS(file)

# All considered variables
var.t.lag <- c('w_avg_tx_anom', 'w_avg_tn_anom', 'w_avg_Tind95', 'w_avg_Tind5', 
               'w_avg_hu_anom', 'w_avg_rr_anom', 'w_avg_fg_anom', 
               'w_avg_hu.ind95', 'w_avg_hu.ind5', 'w_avg_rr.ind95', 'w_avg_rr.ind5',
               'w_avg_fg.ind95', 'w_avg_fg.ind5', 'w_avg_o3_anom', 'w_avg_o3.ind95',
               'w_avg_o3.ind5', 'w_avg_pm10_anom', 'w_avg_pm10.ind95',
               'w_avg_pm10.ind5', 'w_avg_pm2p5_anom', 'w_avg_pm2p5.ind95',
               'w_avg_pm2p5.ind5', 'w_avg_no2_anom', 'w_avg_no2.ind95',
               'w_avg_no2.ind5')
vars.nlag <- var.t.lag
vars.lag  <- paste0(vars.nlag, '_l1')
vars.loc  <- c('long','lat')
vars.time <- c('Summer', 'Winter', 'Spring', 'Fall')
vars      <- c(vars.nlag, vars.lag, vars.loc, vars.time)

# Read XGBoost fit
filexgb <- 'Results/XGB/xgb_fit.rds'
xgb.fit <- readRDS(file = filexgb)

# Predict function xgboost
predict.xgb <- function(model, newdata){
  exp(predict(object = model, 
              newdata =  xgb.DMatrix(data = newdata %>% as.matrix())))
}  

# Predictions
df$xgbDeaths <- predict.xgb(model = xgb.fit, 
                            newdata = df %>% select(all_of(vars))) * df$bDeaths


##### 2) Visualisation ----

# Select region
region <- 'ES511'

# Filter data set
dfr <- df %>% dplyr::filter(Region == region)

# Plot
ggplot(dfr) + 
  theme_bw(base_size = 15) + 
  geom_line(aes(x = Date, y = Deaths), col = 'gray80') + 
  geom_line(aes(x = Date, y = bDeaths), col = ared, linewidth = 0.75) +
  geom_line(aes(x = Date, y = xgbDeaths), col = ablue, linewidth = 0.75) 
  
