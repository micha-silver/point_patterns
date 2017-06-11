library(rhdf5)
library(dplyr)
setwd('/home/micha/Studies/Courses/Geostatistics-Tal/Project/')
cml_file <- "CML/cml_bavaria_v1.h5"

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
  # Set link point at middle between antennas
  lat <- (attrs$site_a_latitude + attrs$site_b_latitude)/2
  lon <- (attrs$site_a_longitude + attrs$site_b_longitude)/2
  meta <- as.data.frame(cbind(link_id, lon, lat))
  colnames(meta) <- c('link_id','lon','lat')
  link_meta[[i]] <- meta 
  
  # Make a list of all daily data for this cml
  dset <- as.data.frame(cbind(rr_dset, dt_dset))
  colnames(dset) <- c('rain_rate', 'date_time')
  dset$date_time <- as.POSIXct(dset$date_time, origin="1970-01-01", tz='GMT')
  cml_daily <- dset %>% 
      group_by(date_time = cut(date_time, breaks="1 day")) %>%
      summarize(rain_rate = mean(rain_rate)*24)
  
  link_id <- rep(link_id, length(cml_daily$date_time))
  cml_daily['link_id'] <- link_id
  link_data[[i]] <-  cml_daily
} 

# Convert lists to DF
link_meta <- do.call(rbind, link_meta)
link_data <- do.call(rbind, link_data)
H5close()
print(paste("Found: ",length(link_meta$link_id), "links"))
print(paste("Aggregated: ", length(link_data$link_id), "rows of daily data"))