# Install hdf5 package for reading NDF5 formatted CML data:
#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")

library(rhdf5)  # Library to read HDF5 formatted CML data
library(dplyr)
library(spdep)
library(gstat)    # library with krige functions
library(raster)
library(RColorBrewer)
library(maptools)
library(automap)
library(rgdal)

# ---------------------------------
# Initial setup
# ---------------------------------
setwd('/home/micha/Studies/Courses/Geostatistics-Tal/Project/')
cml_file  <- "CML/cml_bavaria_v1.h5"
gauge_file  <- 'gauge_data/gauge_data_daily.csv'
gauge_meta_file <- 'gauge_data/gauge_metadata.csv'
# Load administrative boundaries for plotting
gadm_1 <- readRDS('GIS/DEU_adm1.rds')
cities <- shapefile('GIS/cities_de.shp')
# Filter out only Bavaria
cities <- subset(cities, cities$adm1_code==2)
# List of date strings
datestrs = c('2016-07-01', '2016-07-02',
             '2016-07-03', '2016-07-04',
             '2016-07-05', '2016-07-06',
             '2016-07-07', '2016-07-08',
             '2016-07-09', '2016-07-10',
             '2016-07-11', '2016-07-12',
             '2016-07-13', '2016-07-14',
             '2016-07-15', '2016-07-16',
             '2016-07-17', '2016-07-18',
             '2016-07-19', '2016-07-20')

radar_dir <- '/home/micha/Studies/Research/IMAP/Data/DE/RADOLAN/daily/'

# Coordinate reference Systems
# Original data are in Long/Lat
# Transform to ETRS LAEA
proj_wgs84 = CRS("+init=epsg:4326")
proj_laea <- CRS("+init=epsg:3035")

# source functions file
source('Point_Patterns/geostats_Final_functions.R')

# ---------------------------------
# Load CML data
# ---------------------------------
# Open HDF5 file of CML data, and get list of links (groups in HDF5)
grp_list <- filter(h5ls(cml_file), grepl("product_0", name))
dt_col <- rep('time', length(grp_list$group))
rr_col <- rep('rain rate', length(grp_list$group))
dt_list <- paste(grp_list$group, grp_list$name, dt_col, sep="/")
rr_list <- paste(grp_list$group, grp_list$name, rr_col, sep="/")

num_links <- length(rr_list)
link_data <- list()
link_meta <- list()

# Loop thru all links, 
# get meta_data and precip data for each 
for (i in 1:num_links) {
  rr_dset <- h5read(cml_file, rr_list[i])
  dt_dset <- h5read(cml_file, dt_list[i])
  attrs <- h5readAttributes(cml_file, grp_list[i,]$group)
  # Create list of metadata for all cmls
  link_id <- attrs$cml_id
  # Set link point at *middle* between antennas
  lat <- (attrs$site_a_latitude + attrs$site_b_latitude)/2
  lon <- (attrs$site_a_longitude + attrs$site_b_longitude)/2
  meta <- data.frame(link_id, lon, lat)
  colnames(meta) <- c('link_id','lon','lat')
  link_meta[[i]] <- meta 
  
  # Make a list of all daily data for this cml
  dset <- as.data.frame(cbind(rr_dset, dt_dset))
  colnames(dset) <- c('rain_rate', 'date_time')
  dset$date_time <- as.POSIXct(dset$date_time, 
                               origin="1970-01-01", tz='GMT')
  cml_daily <- dset %>% 
    group_by(date_time = cut(date_time, breaks="1 day")) %>%
    summarize(rain_rate = mean(rain_rate)*24)
  
  link_id_list <- rep(link_id, length(cml_daily$date_time))
  cml_daily['link_id'] <- link_id_list
  link_data[[i]] <-  cml_daily
} 

# Convert lists to DF
link_meta <- do.call(rbind, link_meta)
link_data <- do.call(rbind, link_data)
link_data <- as.data.frame(link_data)
link_data <- merge(link_data, link_meta, all.x=FALSE)
link_data <- na.omit(link_data)

H5close()
print(paste("Found: ",length(link_meta$link_id), "links"))
print(paste("Aggregated: ", length(link_data$link_id),
            "rows of daily CML data"))

# ---------------------------------
# Load gauge data
# ---------------------------------
# Read CSV files of gauge data and metadata 
data_cols <- c('station_id','date_time','quality',
               'obs_precip', 'precip_ind','snow','eor')
gauge_data <- read.csv(gauge_file, col.names=data_cols)
meta_cols <- c('station_id','from_date','to_date','elevation',
               'latitude','longitude',
               'stn_name','province')
gauge_meta <- read.csv(gauge_meta_file, col.names=meta_cols)

# Join spatial data with gauge data
gauge_data <- merge(gauge_data, gauge_meta, 
                    by='station_id', all.x=FALSE)
# Make sure to clean out NA or < 0 (unknown values)
gauge_data <- na.omit(gauge_data)
gauge_data <- filter(gauge_data, obs_precip>=0)
print(paste("Found: ",length(gauge_meta$station_id), "gauges"))
print(paste("Found: ", length(gauge_data$station_id),
            "rows of daily gauge data"))

# ---------------------------------
# Make SPDF from data frames
# ---------------------------------
# First SpatialPointsDataFrame from the Link data
link_coords <- cbind(link_data$lon, link_data$lat)
coordinates(link_data) <- link_coords
proj4string(link_data) <- proj_wgs84
links_bbox <- link_data@bbox

# also spatial points data frame of gauges
gauge_coords <- cbind(gauge_data$longitude, gauge_data$latitude)
coordinates(gauge_data) <- gauge_coords
proj4string(gauge_data) <- proj_wgs84
# and crop to links BBOX
gauge_data <- crop(gauge_data, extent(links_bbox))

# Transform all spatial data to European LAEA projection, EPSG:3035
cities <- spTransform(cities, proj_laea)
gadm_1 <- spTransform(gadm_1, proj_laea)
gauge_data <- spTransform(gauge_data, proj_laea)
link_data <- spTransform(link_data, proj_laea)

# Plot inital map
plot_map(cities, gadm_1, gauge_data, link_data)

# Setup grid for kriging results
grd <- prepare_krige_grid(link_data)

# ---------------------------------
# Begin process - one day
# ---------------------------------
for (i in seq(12,13)) {
  datestr = datestrs[i]
  # Extract data for one day
  data_1day <- slice_data(link_data, gauge_data, datestr)
  links_1day <- data_1day$links_1day
  gauges_1day <- data_1day$gauges_1day
  # Moran's I for autocorrelation
  Moran_I <- calculate_moran(links_1day, datestr)
  
  # Create a variogram and fitted model
  vg_fit <- prepare_variogram(links_1day, datestr)  

  # Do Ordinary Kriging
  perform_ok(links_1day, vg_fit, gauges_1day, grd, datestr)
  
  # Do Kriging with External Drift
  radar_file <- paste0(radar_dir,'rw_',gsub('-','',datestr), '.asc')
  perform_ked(links_1day, gauges_1day, radar_file, grd, datestr)
}