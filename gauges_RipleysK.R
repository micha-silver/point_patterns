library('spatstat')
library('dplyr')

# Load data and metadata
setwd(paste0('/home/micha/Studies/Courses',
             '/Geostatistics-Tal/Final Project/Point pattern statistics/'))
data_cols <- c('station_id','date_time','relialbility',
               'rain_bool','obs_precip','novalue','eor')
gauge_data <- read.csv('gauge_data.csv', col.names=data_cols)
meta_cols <- c('station_id','height','latitude','longitude',
               'from_date','stn_name','province')
gauge_metadata <- read.csv('gauge_metadata.csv', col.names=meta_cols)
head(gauge_data)
head(gauge_metadata)

# Get one hour, and attach metadata
gauge_data_2016110703 <- filter(gauge_data, date_time=='2016110703')
head(gauge_data_2016110703)
gauges <- merge(gauge_data_2016110703, gauge_metadata, by='station_id', all.y=TRUE)
gauges <- na.omit(gauges)


# Convert to ppp
gauges_owin <- ripras(gauges$longitude, gauges$latitude, shape="rectangle")
gauges_marks <- gauges$obs_precip
gauges_ppp <- ppp(gauges$longitude, gauges$latitude, 
                  window=gauges_owin,
                  marks=gauges_marks)
str(gauges_ppp)
class(gauges_ppp)

# Check overall intensity
print(summary(gauges_ppp)$intensity)
# Check quadrat count
plot(quadratcount(gauges_ppp, nx=9, ny=5))

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
gauges_l <- Lest(gauges_ppp, correction="Ripley")
plot(gauges_l, main="Ripley's L function")
plot(envelope(gauges_ppp, Lest, 49), main="Ripley's L with Envelope")