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
cat("\n============================\n")
print("Initial Setup")
setwd('/home/micha/Studies/Courses/Geostatistics-Tal/Project/')
out_dir <- paste0(getwd(),'/Report/Output/')
cml_file  <- "CML/cml_bavaria_v1.h5"
# Load administrative boundaries for plotting
gadm_1 <- readRDS('GIS/DEU_adm1.rds')
cities <- shapefile('GIS/cities_de.shp')
# Filter out only Bavaria
cities <- subset(cities, cities$adm1_code==2)
# List of date strings
datestrs = c('2016-06-21', '2016-06-22', '2016-06-23', '2016-06-24',
             '2016-06-25', '2016-06-26', '2016-06-27', '2016-06-28',
             '2016-06-29', '2016-06-30', '2016-07-01', '2016-07-02',
             '2016-07-03', '2016-07-04', '2016-07-05', '2016-07-06',
             '2016-07-07', '2016-07-08', '2016-07-09', '2016-07-10',
             '2016-07-11', '2016-07-12', '2016-07-13', '2016-07-14',
             '2016-07-15', '2016-07-16', '2016-07-17', '2016-07-18',
             '2016-07-19', '2016-07-20', '2016-07-21', '2016-07-22',
             '2016-07-23', '2016-07-24', '2016-07-25', '2016-07-26',
             '2016-07-27', '2016-07-28', '2016-07-29', '2016-07-20')

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
# Make SPDF from data frame
# ---------------------------------
# First SpatialPointsDataFrame from the Link data
link_coords <- cbind(link_data$lon, link_data$lat)
coordinates(link_data) <- link_coords
proj4string(link_data) <- proj_wgs84
links_bbox <- link_data@bbox

# Transform all spatial data to European LAEA projection, EPSG:3035
cities <- spTransform(cities, proj_laea)
gadm_1 <- spTransform(gadm_1, proj_laea)
link_data <- spTransform(link_data, proj_laea)

# Plot inital map
plot_map(cities, gadm_1, link_data)

# Setup grid for kriging results
grid_result <- prepare_krige_grid(link_data)
grd <- grid_result$grid
validation_points <- grid_result$validation


# ---------------------------------
# Begin process - one day
# ---------------------------------
for (i in c(5,6,22,23,24,35)) {
  datestr = datestrs[i]
  cat("\n============================\n")
  print(paste("Starting Process for date:", datestr))
  # Extract data for one day
  links_1day <- slice_data(link_data, datestr)
  # Get radar for extracting validation points
  radar <- load_radar(datestr, grd)
  
  # Moran's I for autocorrelation
  Moran_I <- calculate_moran(links_1day, datestr)
  
  # Create a variogram and fitted model
  vg_fit <- prepare_variogram(links_1day, datestr)  

  # Do Ordinary Kriging
  OK_result <- perform_ok(links_1day, vg_fit, grd, 
                           validation_points, radar, datestr)
  # Produce plots
  OK_scatter(radar, OK_result )
  
  # Do Kriging with External Drift
  KED_results <- perform_ked(links_1day, grd, radar,
                             validation_points, datestr)
  KED_scatter(KED_results, datestr)
  check_KED(links_1day, radar)
  
}
