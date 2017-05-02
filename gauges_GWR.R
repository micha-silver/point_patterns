library('spatstat')
library('dplyr')
library('sp')
library("raster")
library('RColorBrewer')
library('spgwr')
library('dismo')

# ---------------------------------
# Load data,
# ---------------------------------
# Load data and metadata
setwd(paste0('/home/micha/Studies/Courses',
             '/Geostatistics-Tal/Final Project/'))

# WHich data interval?
data_interval <- 'daily'

if (data_interval == 'daily') {
  date_str <- '20170112'
  data_file <- paste0('gauge_data/','gauge_data_daily.csv')
  data_cols <- c('station_id','date_time','quality',
                 'obs_precip', 'precip_ind','snow','eor')
} else {
  date_str <- '2017011201'
  data_file <- paste0('gauge_data/','gauge_data_hourly.csv')  
  data_cols <- c('station_id','date_time','quality',
                'precip_ind', 'obs_precip','form','eor')
}
gauge_data <- read.csv(data_file, col.names=data_cols)

meta_cols <- c('station_id','from_date','to_date','elevation',
               'latitude','longitude',
               'stn_name','province')
gauge_metadata <- read.csv('gauge_data/gauge_metadata.csv', col.names=meta_cols)
head(gauge_data)
head(gauge_metadata)

# Get one day (or hour), and attach metadata
gauge_data_filtered <- filter(gauge_data, date_time==date_str)
head(gauge_data_filtered)
gauges <- merge(gauge_data_filtered, gauge_metadata, by='station_id', all.y=TRUE)
# Make sure to clean out NA or < 0 (unknown values)
gauges <- na.omit(gauges)
gauges <- filter(gauges, obs_precip>=0)


# ---------------------------------
# Prepare SpatialPointsDataFrame
# ---------------------------------
xy_coords <- cbind(gauges$longitude, gauges$latitude) 
proj_wgs84 = CRS("+init=epsg:4326")
coordinates(gauges) <- xy_coords
proj4string(gauges) <- proj_wgs84
# Now gauges should be a SpatialPointsDataFrame
class(gauges)

# Show gauges on a map
gadm_1 <- readRDS('GIS/DEU_adm1.rds')
#gadm_2 <- readRDS('GIS/DEU_adm2.rds')
png('bavaria_gauges.png', width=600, height=600)
bg <- gmap('Bavaria')
# Reproject gadm borders and gauges to EPSG:3857
gadm_1_reproj <- spTransform(gadm_1, bg@crs)
gauges_reproj <- spTransform(gauges, bg@crs)
plot(bg, width=3, height=4)
plot(gadm_1_reproj, fill=F, border='black', lwd=2, add=T)
plot(gauges_reproj, pch=18, col="blue", cex=gauges$obs_precip/5, add=T)
title(main="Rain Gauges - Bavaria", sub="11/01/2017")
dev.off()

# Read in aspect raster
aspect_bavaria <- raster('GIS/aspect_bavaria.tif')
# Extract aspect from raster at gauge locations, 
# and add values to new data column
gauges$aspect <- extract(aspect_bavaria, gauges)
# Now we have three possible predictors: elevation, latitude, aspect
# dependant variable is obs_precip
head(gauges$aspect)

# ---------------------------------
# OLS
# ---------------------------------
# Try three predictor variables
gauges_ols <- lm(obs_precip ~ elevation + aspect + latitude, data=gauges)
summary(gauges_ols)
# Print confidence intervals
print(confint(gauges_ols))

# Now remove latitude, use only two 
gauges_ols <- lm(obs_precip ~ elevation + aspect, data=gauges)
summary(gauges_ols)
plot(gauges_ols, which=3)
print(confint(gauges_ols))

# Scatterplot of observations vs fitted
plot(gauges$obs_precip, gauges_ols$fitted.values, pch=16, col="Blue",
     xlab="Observed", ylab="Model", main="Observed vs Modeled")
abline(lm(gauges$obs_precip~gauges_ols$fitted.values), col="Red")

# Check for spatial correlation of residuals
gauges_resids<-residuals(gauges_ols)
clrs <- brewer.pal(5, "RdBu")
gauges_resids_sp <- SpatialPointsDataFrame(data.frame(gauges_resids), 
                                           coords=xy_coords, 
                                           proj4string = proj_wgs84)
plot(gauges_resids_sp, cuts=cut(gauges_resids_sp$gauges_resids, breaks=5), 
     col=clrs, pch=18, cex=2)
plot(gadm_1, fill=F, border='black', lwd=2, add=T)

# Check for normal distribution of residuals
hist(gauges_resids_sp$gauges_resids, col="orange", 
     breaks=12, main="Frequency Distribution of Residuals")

# ---------------------------------
# GWR
# ---------------------------------
# Use only two predictors
# First determine the bandwidth
bw <- gwr.sel(obs_precip ~ elevation + aspect, 
              data = gauges, longlat=T)
print(bw)
gauges_gwr <- gwr(obs_precip ~ elevation + aspect, 
                  data = gauges, longlat=T, bandwidth=bw,
                  se.fit=T, hatmatrix=T)
str(gauges_gwr$SDF)
# Show R squared values
print(gauges_gwr$SDF@data$localR2)

# Analyze results
gauge_results <- as.data.frame(gauges_gwr$SDF)
gauge_results <- SpatialPointsDataFrame(gauge_results, 
                                       coords=xy_coords, 
                                       proj4string = proj_wgs84)
# Plot predicted precipitation rate
plot(gauge_results, pch=16, col="Blue", cex=gauge_results$pred/4)
plot(gadm_1, fill=F, border='black', lwd=2, add=T)

# Plot overall standard error
plot(gauge_results, pch=5, col="blue", cex=gauge_results$gwr.e, lwd=2)
plot(gadm_1, fill=F, border='black', lwd=2, add=T)

# Plot coefficients of each of the predictors
opar <- par(mfrow=c(1,2))
plot(gauge_results, pch=17, col="brown", cex=gauge_results$elevation*100,
     main="Elevation Coefficient")
plot(gadm_1, fill=F, border='black', lwd=2, add=T)
plot(gauge_results, pch=15, col="purple", cex=gauge_results$aspect*200,
     main="Aspect coefficient")
plot(gadm_1, fill=F, border='black', lwd=2, add=T)

# Plot t-values for both coefficients
# t-value = coefficient / standard error
# If it's around 0, highly significant, > +- 4 not significant
elev_t <- gauge_results$elevation / gauge_results$elevation_se
aspect_t <- gauge_results$aspect / gauge_results$aspect_se
gauge_signif <- cbind(gauge_results, elev_t, aspect_t)
#gauge_signif <- SpatialPointsDataFrame(gauge_signif, 
#                                        coords=xy_coords, 
#                                        proj4string = proj_wgs84)
opar <- par(mfrow=c(1,2))
clrs2=c("green","red","purple","red","green") 
breaks=c(min(gauge_signif$elev_t),-4,0,4,max(gauge_signif$elev_t)) 
plot(gauge_signif, pch=17, col="brown", cex=gauge_results$elevation*100,
     main="Elevation Coefficient")
plot(gadm_1, fill=F, border='black', lwd=2, add=T)
plot(gauge_signif, pch=16, col=clrs2, cex=2,
     main="Elevation t-value")
plot(gadm_1, fill=F, border='black', lwd=2, add=T)

opar <- par(mfrow=c(1,2))
breaks=c(min(gauge_signif$aspect_t),-4,0,4,max(gauge_signif$aspect_t)) 
plot(gauge_signif, pch=15, col="purple", cex=gauge_results$aspect*200,
     main="Aspect Coefficient")
plot(gadm_1, fill=F, border='black', lwd=2, add=T)
plot(gauge_signif, pch=16, col=clrs2, cex=2,
     main="Aspect t-value")
plot(gadm_1, fill=F, border='black', lwd=2, add=T)
