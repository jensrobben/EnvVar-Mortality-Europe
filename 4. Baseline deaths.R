##### 0) Preliminary settings ----

# Clear environment
rm(list = ls())
gc()

# Required packages
packages <- c('lubridate', 'dplyr', 'rstudioapi', 'eurostat', 'ggplot2', 'utils',
              'ISOweek', 'mgcv', 'sf', 'spdep', 'stats', 'stringr', 'tidyr')
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


##### 3) Load population exposures and preprocess ----

# Download annual population exposures in broad age categories (January 1)
P.xt.all <- get_eurostat(id = 'demo_r_pjanaggr3', cache = FALSE,
                         compress_file = FALSE, time_format = 'raw')
colnames(P.xt.all) <- c('Freq','Unit','Sex','Age','Region','Time', 'Pop')

# Annual exposures- countries of interest
P.xt.c <- P.xt.all %>% dplyr::filter(Region %in% regions &
                                       Age %notin% c('UNK', 'TOTAL')) %>%
  dplyr::select(-c('Freq','Unit')) %>%
  mutate(ISOYear = as.integer(substr(Time, 1, 4)),
         Age = factor(Age, levels = c('Y_LT15', 'Y15-64','Y_GE65')),
         .before = 1) %>%
  select(-c('Time')) %>% arrange(Sex, ISOYear, Age) %>%
  dplyr::filter(!grepl(':', Pop), ISOYear <= year(t_max) + 1, ISOYear >= year(t_min))

# Focus on females and age group 65+
P.xt.c <- P.xt.c %>% dplyr::filter(Age == 'Y_GE65', Sex == "F") %>%
  select(-c(Age, Sex))

# Create weekly exposures from population count
E.xt.c <- P.xt.c %>% group_by(Region) %>% 
  arrange(ISOYear) %>% reframe(ISOYear, 
                               'Expo' = c((Pop[-length(ISOYear)] + Pop[-1])/2, 
                                          Pop[length(ISOYear)])/52.1775) %>%
  dplyr::filter(ISOYear <= year(t_max))

# Add to death data set
df.d <- df.d %>% left_join(E.xt.c, by = c('ISOYear', 'Region')) 

##### 4) Fit spatially smoothed Poisson GLM ----

# Adjacency matrix
adjmat <- nb2mat(poly2nb(sf.ctry[which(sf.ctry$NUTS_ID %in% regions),'geometry']), 
                 zero.policy = TRUE)
adjmat[adjmat > 0] <- -1
dimnames(adjmat) <- list(sf.ctry[which(sf.ctry$NUTS_ID %in% regions),]$NUTS_ID, 
                         sf.ctry[which(sf.ctry$NUTS_ID %in% regions),]$NUTS_ID)
diag(adjmat) <- sapply(1:nrow(adjmat), function(x) sum(adjmat[,x] != 0))
nr <- nrow(adjmat)

# Data
df <- df.d %>%
  mutate('fsin52' = sin(2*pi*ISOWeek/52.1775),
         'fcos52' = cos(2*pi*ISOWeek/52.1775),
         'fsin26' = sin(4*pi*ISOWeek/52.1775),
         'fcos26' = cos(4*pi*ISOWeek/52.1775))
df$Region <- factor(df$Region, levels = sort(unique(df$Region)))

# Save Exposure information
dfExpo <- df %>% dplyr::select('Date', 'ISOYear', 'ISOWeek', 'Region', 'Expo')
file   <- paste0('Results/Data/Expo_NUTS3_ES.rds')
saveRDS(dfExpo, file)

# Formulas
formula0 <- Deaths ~ -1 + Region
formula1 <- Deaths ~ -1 + ISOYear:Region
formula2 <- Deaths ~ -1 + fsin52:Region
formula3 <- Deaths ~ -1 + fcos52:Region
formula4 <- Deaths ~ -1 + fsin26:Region
formula5 <- Deaths ~ -1 + fcos26:Region

# Model matrices
X0 <- model.matrix(formula0, data = df)
X1 <- model.matrix(formula1, data = df)
X2 <- model.matrix(formula2, data = df)
X3 <- model.matrix(formula3, data = df)
X4 <- model.matrix(formula4, data = df)
X5 <- model.matrix(formula5, data = df)

# Add model matrices to data frame
df$M0 <- X0
df$M1 <- X1
df$M2 <- X2
df$M3 <- X3
df$M4 <- X4
df$M5 <- X5

# Formula
formula <- Deaths ~ -1 + M0 + M1 + M2 + M3 + M4 + M5

# Remove to save memory
rm(X0,X1,X2,X3,X4,X5)
gc()

# Fit
fit <- gam(formula, offset = log(df$Expo), data = df,
           family = poisson(link = 'log'),
           control = list(epsilon = 10^(-6), trace = TRUE,
                          newton = list(conv.tol = 1e-5), nthreads = 8),
           drop.intercept = TRUE,
           paraPen = list(M0 = list(adjmat, rank = nr, sp = -1),
                          M1 = list(adjmat, rank = nr, sp = -1),
                          M2 = list(adjmat, rank = nr, sp = -1),
                          M3 = list(adjmat, rank = nr, sp = -1),
                          M4 = list(adjmat, rank = nr, sp = -1),
                          M5 = list(adjmat, rank = nr, sp = -1)))

# Save model fit
saveRDS(fit, paste0('Results/GAM/GAM.rds'))

# Save baseline deaths
df.fitb <- data.frame(df[,c('Date', 'ISOYear', 'ISOWeek', 'Region', 'Deaths')],
                      'bDeaths' = fit$fitted.values)
saveRDS(df.fitb, 'Results/Data/bDeaths_NUTS3_ES.rds')

##### 5) Visualizations ----

# Read fit
fit <- readRDS('Results/GAM/GAM.rds')

### Parameter coefficients
coef.new    <- fit$coefficients
df.coef.new <- data.frame('Region' = sub(".*Region", "",names(coef.new)), 
                          'Param' = coef.new, 
                          'Type' = str_extract(names(coef.new), 
                                               "(ISOYear|cos\\d+|sin\\d+)"))
rownames(df.coef.new) <- NULL
df.coef.new[which(is.na(df.coef.new$Type)),'Type'] <- 'INTERCEPT'
df.coef.new  <- df.coef.new %>% pivot_wider(names_from = 'Type', values_from = 'Param')
df.coef.new$NewInt <- df.coef.new$INTERCEPT + df.coef.new$ISOYear*2013
sf.ctry.plot <- df.coef.new %>% full_join(sf.ctry, by = c('Region' = 'NUTS_ID'))

# Plots
p0 <- ggplot() + ggtitle(bquote(beta[0]))  +
  theme_bw(base_size = 12.5) + 
  geom_sf(data = sf.ctry.plot$geometry, aes(fill = sf.ctry.plot$INTERCEPT,
                                            col  = sf.ctry.plot$INTERCEPT,
                                            group = sf.ctry.plot$CNTR_CODE)) + 
  scale_fill_viridis_c(option = "magma", name = bquote(beta[0]~' '), 
                       na.value = 'white', direction = -1) +
  scale_color_viridis_c(option = "magma", name = bquote(beta[0]~' '), 
                        na.value = 'white', direction = -1) + 
  theme(legend.position = 'bottom',
        legend.key.width=unit(1,"cm"),
        legend.title = element_text(margin = margin(t = -1.25, unit = "lines"), 
                                    size = 14),
        plot.title = element_text(size=20))

p1 <- ggplot() + ggtitle(bquote(beta[1]))  +
  theme_bw(base_size = 12.5) + 
  geom_sf(data = sf.ctry.plot$geometry, aes(fill = sf.ctry.plot$ISOYear,
                                            col  = sf.ctry.plot$ISOYear,
                                            group = sf.ctry.plot$CNTR_CODE)) + 
  scale_fill_viridis_c(option = "magma", name = bquote(beta[1]~' '), 
                       na.value = 'white', direction = -1) +
  scale_color_viridis_c(option = "magma", name = bquote(beta[1]~' '), 
                        na.value = 'white', direction = -1) + 
  theme(legend.position = 'bottom',
        legend.key.width=unit(1,"cm"),
        legend.title = element_text(margin = margin(t = -1.25, unit = "lines"), 
                                    size = 14),
        plot.title = element_text(size=20))

p2 <- ggplot() + ggtitle(bquote(beta[2]))  +
  theme_bw(base_size = 12.5) + 
  geom_sf(data = sf.ctry.plot$geometry, aes(fill = sf.ctry.plot$sin52,
                                            col  = sf.ctry.plot$sin52,
                                            group = sf.ctry.plot$CNTR_CODE)) + 
  scale_fill_viridis_c(option = "magma", name = bquote(beta[2]~' '), 
                       na.value = 'white', direction = -1) +
  scale_color_viridis_c(option = "magma", name = bquote(beta[2]~' '), 
                        na.value = 'white', direction = -1) + 
  theme(legend.position = 'bottom',
        legend.key.width=unit(1,"cm"),
        legend.title = element_text(margin = margin(t = -1.25, unit = "lines"), 
                                    size = 14),
        plot.title = element_text(size=20))

p3 <- ggplot() + ggtitle(bquote(beta[3]))  +
  theme_bw(base_size = 12.5) + 
  geom_sf(data = sf.ctry.plot$geometry, aes(fill = sf.ctry.plot$cos52,
                                            col  = sf.ctry.plot$cos52,
                                            group = sf.ctry.plot$CNTR_CODE)) + 
  scale_fill_viridis_c(option = "magma", name = bquote(beta[3]~' '), 
                       na.value = 'white', direction = -1) +
  scale_color_viridis_c(option = "magma", name = bquote(beta[3]~' '), 
                        na.value = 'white', direction = -1) + 
  theme(legend.position = 'bottom',
        legend.key.width=unit(1,"cm"),
        legend.title = element_text(margin = margin(t = -1.25, unit = "lines"),
                                    size = 14),
        plot.title = element_text(size=20)) 

p4 <- ggplot() + ggtitle(bquote(beta[4]))  +
  theme_bw(base_size = 12.5) + 
  geom_sf(data = sf.ctry.plot$geometry, aes(fill = sf.ctry.plot$sin26,
                                            col  = sf.ctry.plot$sin26,
                                            group = sf.ctry.plot$CNTR_CODE)) + 
  scale_fill_viridis_c(option = "magma", name = bquote(beta[4]~' '), 
                       na.value = 'white', direction = -1, n.breaks = 4) +
  scale_color_viridis_c(option = "magma", name = bquote(beta[4]~' '), 
                        na.value = 'white', direction = -1, n.breaks = 4) + 
  theme(legend.position = 'bottom',
        legend.key.width=unit(1,"cm"),
        legend.title = element_text(margin = margin(t = -1.25, unit = "lines"), 
                                    size = 14),
        plot.title = element_text(size=20))

p5 <- ggplot() + ggtitle(bquote(beta[5]))  +
  theme_bw(base_size = 12.5) + 
  geom_sf(data = sf.ctry.plot$geometry, aes(fill = sf.ctry.plot$cos26,
                                            col  = sf.ctry.plot$cos26,
                                            group = sf.ctry.plot$CNTR_CODE)) + 
  scale_fill_viridis_c(option = "magma", name = bquote(beta[5]~' '), 
                       na.value = 'white', direction = -1) +
  scale_color_viridis_c(option = "magma", name = bquote(beta[5]~' '), 
                        na.value = 'white', direction = -1) + 
  theme(legend.position = 'bottom',
        legend.key.width=unit(1,"cm"),
        legend.title = element_text(margin = margin(t = -1.25, unit = "lines"), 
                                    size = 14),
        plot.title = element_text(size=20))

### Fitted death rates
mxtfit <- predict(fit, newdata = df, type = 'response')
dfpred <- df %>% mutate('qtr.obs'  = 1 - exp(-Deaths/Expo),
                        'qtr.fit.new'  = 1 - exp(-mxtfit))

# Regions
region1 <- 'ES111'
region2 <- 'ES300'
region3 <- 'ES511'

# Range
range <- dfpred %>% dplyr::filter(Region %in% c(region1, region2, region3)) %>%
  pull(qtr.obs) %>% range()

# Plots
p1 <- ggplot(dfpred %>% dplyr::filter(Region == region1)) + 
  ggtitle(paste0(region1,': ', sf.ctry$NUTS_NAME[sf.ctry$NUTS_ID == region1])) + 
  theme_bw(base_size = 15) + ylab(bquote(q['t,w']^'(r)')) + 
  geom_line(aes(x = Date, y = qtr.obs), col = 'gray70', alpha = 0.5) + 
  geom_line(aes(x = Date, y = qtr.fit.new), col = ared, linewidth = 1) +
  coord_cartesian(ylim = range)
p1

p2 <- ggplot(dfpred %>% dplyr::filter(Region == region2)) + 
  ggtitle(paste0(region2,': ', sf.ctry$NUTS_NAME[sf.ctry$NUTS_ID == region2])) + 
  theme_bw(base_size = 15) + ylab(bquote(q['t,w']^'(r)')) + 
  geom_line(aes(x = Date, y = qtr.obs), col = 'gray70', alpha = 0.5) + 
  geom_line(aes(x = Date, y = qtr.fit.new), col = ared, linewidth = 1) +
  coord_cartesian(ylim = range)
p2

p3 <- ggplot(dfpred %>% dplyr::filter(Region == region3)) + 
  ggtitle(paste0(region3,': ', sf.ctry$NUTS_NAME[sf.ctry$NUTS_ID == region3])) + 
  theme_bw(base_size = 15) + ylab(bquote(q['t,w']^'(r)')) + 
  geom_line(aes(x = Date, y = qtr.obs), col = 'gray70', alpha = 0.5) + 
  geom_line(aes(x = Date, y = qtr.fit.new), col = ared, linewidth = 1) + 
  coord_cartesian(ylim = range)
p3
