library('dplyr')
library('spdep')
library('geoR')
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

# Get one day (or hour), and attach metadata
gauge_data_filtered <- filter(gauge_data, date_time==date_str)
head(gauge_data_filtered)
gauges <- merge(gauge_data_filtered, gauge_metadata, by='station_id', all.y=TRUE)
# Make sure to clean out NA or < 0 (unknown values)
gauges <- na.omit(gauges)
gauges <- filter(gauges, obs_precip>=0)

# Weights matrix
distances <- as.matrix(dist(gauges))
inv_dist = 1/distances
diag(inv_dist) <- 0
#inv_dist_squared <- 1/distances^2
print(inv_dist[1:5,1:5])
wts <- mat2listw(inv_dist)
#wts_squared <- mat2listw(inv_dist_squared)

# Prepare SPDF
xy_coords <- cbind(gauges$longitude, gauges$latitude) 
proj_wgs84 = CRS("+init=epsg:4326")
coordinates(gauges) <- xy_coords
proj4string(gauges) <- proj_wgs84

# Moran's I
I <- moran.test(gauges$obs_precip, wts)
print(I)

I_mc <- moran.mc(gauges$obs_precip, wts, nsim=9999)
# Get value of 1.96 Stdev of the Monte Carlo distribution
sd2 <- 1.96 * sd(I_mc$res)
hist(I_mc$res, xlab="Monte Carlo runs of Moran's I", 
     breaks=50, col="lightgrey")
abline(v=I$estimate[1], col="red", lwd=3)
abline(v=I$estimate[2], col="orange", lwd=2)
abline(v=sd2, col="grey", lwd=2, lty="dashed")
abline(v=-sd2, col="steelblue4", lwd=2, lty="dashed")

# Geary's C
C <- geary.test(gauges$obs_precip, wts)

# Variogram and Envelope
coord_matrix <- cbind(gauges$longitude, gauges$latitude)
vg <- variog(coords=coord_matrix, 
             data=gauges$obs_precip)
plot(vg, pch=16, col="blue")
vg_envel <- variog.mc.env(coords=coord_matrix, 
                          data=gauges$obs_precip, 
                          obj.var=vg, 
                          nsim=999)
plot(vg, envelope=vg_envel, pch=16, col="blue")