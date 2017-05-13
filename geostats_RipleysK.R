library('spatstat')
library('dplyr')

# ---------------------------------
# Point Density Analysis
# ---------------------------------
# Load data and metadata
setwd(paste0('/home/micha/Studies/Courses',
             '/Geostatistics-Tal/Project/'))


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

# Convert to ppp
gauges_owin <- ripras(gauges$longitude, gauges$latitude, shape="rectangle")
gauges_marks <- gauges$obs_precip
gauges_ppp <- ppp(gauges$longitude, gauges$latitude, 
                  window=gauges_owin,
                  marks=gauges_marks)
#str(gauges_ppp)
class(gauges_ppp)

# Check overall intensity
print(summary(gauges_ppp)$intensity)
# Check quadrat count
plot(quadratcount(gauges_ppp, nx=12, ny=9))

# Get point density and do Kolmogorov-Smirnov
gauges_dens <- density(gauges_ppp)
plot(gauges_dens, main='Kernel Density of Gauges')
gauges_ks <- kstest(gauges_ppp, gauges_dens)
print(gauges_ks)
plot(gauges_ks, main="Kolmogorov-Smirnoff")

# Ripley's K test
gauges_fv <- Kest(gauges_ppp, correction=c("best"), var.approx=TRUE)
plot(gauges_fv, main="Ripley's K Test")
plot(envelope(gauges_ppp, Kest, 49), main="Ripley's K with Envelope")
# L test
gauges_l <- Lest(gauges_ppp, correction="Ripley")
#plot(gauges_l, main="Ripley's L function")
plot(envelope(gauges_ppp, Lest, 49), main="Ripley's L with Envelope")

# L test with varying Dmax, i.e. global envelopes
# 5% confidence throughout 
plot(envelope(gauges_ppp, Lest, 19, global=T), main="Ripley's L with Envelope")

