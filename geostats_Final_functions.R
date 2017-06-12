require(dplyr)
require(spdep)
require(gstat)
require(raster)
require(RColorBrewer)
require(maptools)
require(automap)
require(rgdal)

# ---------------------------------
# Plot initial map
# ---------------------------------
plot_map <- function(cities, gadm_1, gauge_data, link_data) {
  # Clip gauges to region covered by links BBOX
  # first create spatial points data frame of links
 
  # Plot links and gauges
  outpng <- 'Output/gauges_links_map.png'
  png(outpng, width=800, height=600)
  orig_par <- par()
  plot(link_data, col="red", main="Bavaria - Links and Gauges",
       pch=18, cex=1.6)
  plot(gadm_1, border='black', add=T)
  points(gauge_data, col='blue', pch=7, cex=1.0)
  points(cities, pch=16, col='grey30', cex=1.2)
  pointLabel(coordinates(cities), cities$asciiname, cex=0.9)
  lclrs <- c('red','blue','grey30')
  lgnd <- c('CML Links','Gauges','Cities')
  legend('bottomright', col=lclrs, legend=lgnd, 
         lwd=2, cex=0.9, pch=c(18, 7, 16), lty=NA)
  scalebar(40000, type="bar", xy=c(4450000,2700000),
           divs=4, cex=0.9)
  par(orig_par)
  dev.off()
  
}

# ---------------------------------
# Subset one day of data
# ---------------------------------
# Function to slice out one day of data
slice_data <- function(link_data, gauge_data, datestr) {
  # Slice out link data for one day
  links_1day <- link_data[link_data$date_time == datestr,]
  # Remove all NaN values
  links_1day <- na.omit(links_1day)

  # Slice out gauge data for one day
  datestr_clean <- gsub("-","",datestr)
  gauges_1day <- gauge_data[gauge_data$date_time==datestr_clean,]
  # Crop out only the gauges that are in the links region
  gauges_1day <- crop(gauges_1day, extent(links_1day@bbox))
  
  retval <- list("gauges_1day"=gauges_1day, "links_1day"=links_1day)
  return(retval)
}


# ---------------------------------
# Check Point Pattern Dispersion
# Moran's I 
# ---------------------------------
# Function to calculate Moran's I, do Monte Carlo simulations and plot
calculate_moran <- function(links_1day, datestr) {
 
  # Weights matrix as inverse distances between all link points
  distances <- as.matrix(dist(as.data.frame(links_1day)))
  inv_dist = 1/distances
  diag(inv_dist) <- 0
  #print(inv_dist[1:5,1:5])
  wts <- mat2listw(inv_dist)
 
  # Moran's I
  Moran_I <- moran.test(links_1day$rain_rate, wts)
  print(Moran_I)
  I_stat <- Moran_I$estimate[1]
  I_expect <- Moran_I$estimate[2]
  
  # Monte Carlo simulations
  Moran_I_mc <- moran.mc(links_1day$rain_rate, wts, nsim=9999)
  # Get value of 1.96 Stdev of the Monte Carlo distribution
  sd2 <- 1.96 * sd(Moran_I_mc$res)
  
  # Create Histogram plot
  png('Output/Morans_histogram.png', width=800, height=800)
  title <- paste("Moran's I - Distribution of Monte Carlo Simulations", datestr)
  hist(Moran_I_mc$res, xlab="Monte Carlo runs of Moran's I", 
       breaks=50, col="lightgrey", main=title)
  mtext(datestr)
  abline(v=I_stat, col="red", lwd=3)
  abline(v=I_expect, col="orange", lwd=3)
  abline(v=sd2+I_expect, col="steelblue4", lwd=2, lty="dashed")
  abline(v=-sd2+I_expect, col="steelblue4", lwd=2, lty="dashed")
  lclrs <- c('red','orange','steelblue4')
  lgnd <- c('I statistic','Expected Value','2 Std Dev')
  legend('topright', col=lclrs, legend=lgnd, 
         lwd=2, cex=0.9)
  dev.off()
  return(I_stat)
}


# ---------------------------------
# Ordinary Kriging of CML data
# ---------------------------------
# Function to prepare variogram and examine envelope
prepare_variogram <- function(links_1day, datestr) {
  
  vg <- variogram(rain_rate~1, data=links_1day)
  plot(vg)
  # Guess at range and sill
  vg_fit_man <- fit.variogram(vg, 
                              vgm(psill=0.5,c('Exp','Sph','Gau'),range=20000))
  # vg_fit_sph <- fit.variogram(vg, 
  #                             vgm(psill=100,'Sph',range=1.0))
  # vg_fit_gau <- fit.variogram(vg, 
  #                             vgm(psill=100,'Gau',range=1.0))
  # Run autofit function from automap
  vg_fit_auto <- autofitVariogram(rain_rate~1, 
                             input_data=links_1day,
                             model=c("Exp","Sph","Gau"),
                             verbose=TRUE)
  
  # Plot variogram and models
  outpng <- paste0("Output/variogram_", datestr, ".png")
  png(outpng, width=600, height=400)
  y_upper <- max(vg$gamma)*1.1
  x_upper <- max(vg$dist)*1.1
  title <- paste("Variogram_",datestr)
  orig_par <- par()
  plot(gamma~dist, vg, col="blue", main=title,
       xlab="Distance", ylab="Semivariance", 
       ylim=c(0, y_upper), xlim=c(0,x_upper))
  lines(variogramLine(vg_fit_man, x_upper), col="red")
  lines(variogramLine(vg_fit_auto$var_model, x_upper), col="green", lwd=2, lty=3)

  lclrs <- c('red','green')
  lgnd <- c('Manual fit','Auto')
  legend('bottomright', col=lclrs, legend=lgnd, lwd=2, lty=c(1,4,1,2), cex=1.2)
  par(orig_par)
  dev.off()
  # Return the fitted model 
  return(vg_fit_man)
}

prepare_krige_grid <- function(link_data) {
  # Setup grid for output raster, resolution 1000 m.
  # links_1day is already tranformed to LAEA, 
  # be sure to use @coords (NOT lon, lat columns)
  # Extend grid 6000 m. beyond the links region (6 pixels)
  mnx <- min(link_data@coords[,1])-6000
  mxx <- max(link_data@coords[,1])+6000
  mny <- min(link_data@coords[,2])-6000
  mxy <- max(link_data@coords[,2])+6000
  grd <- expand.grid(lon=seq(mnx,mxx,1000), lat=seq(mny,mxy,1000))
  xy_coords <- cbind(grd$lon, grd$lat)
  coordinates(grd) <- xy_coords
  proj4string(grd) <- proj_laea
  gridded(grd) <- TRUE
  return(grd)
}

# Function to preform Ordinary Kriging
perform_ok <- function(links_1day, vg_fit, gauges_1day, grd, datestr) {
 
  # Kriging on grid
  # First: Remove duplicated locations!
  zdist <- zerodist(links_1day)[,1]
  if (length(zdist)>0) {
    # Remove rows identified by zerodist
    links_1day <- links_1day[-zdist,]
  }
  OK_grid <- krige(rain_rate~1, 
                    locations=links_1day, newdata=grd, model=vg_fit)
  
  # Plot krige raster
  clrs <- brewer.pal(8,'BuPu')
  title <- paste("Ordinary Kriging Interpolation ", datestr)
  outpng <- paste0("Output/OK_interp_", links_1day$date_time[1], ".png")
  png(outpng, width=800, height=600)
  orig_par <- par()
  plot(OK_grid, col=clrs, main=title)
  plot(links_1day, col="red", pch=18, cex=1.2, add=T)
  plot(gadm_1, border='black', lwd=1, add=T)
  par(orig_par)
  dev.off()
  
  # Kriging at gauge locations
  OK_gauges <- krige(rain_rate~1, 
                           locations=links_1day, newdata=gauges_1day, model=vg_fit)
  
  # Join results with gauges to get the station_id, and observed precip 
  gauges_obs <- data.frame("station_id"=gauges_1day$station_id, 
                           "obs_precip"=gauges_1day$obs_precip,
                          "x_coord"=gauges_1day@coords[,1], 
                          "y_coord"=gauges_1day@coords[,2])
  links_pred <- data.frame("pred"=OK_gauges$var1.pred,
                           "x_coord"=OK_gauges@coords[,1], 
                           "y_coord"=OK_gauges@coords[,2])
  OK_results <- merge(gauges_obs, links_pred, by=c('x_coord', 'y_coord'))
  
  return(OK_results)
}


OK_scatter <- function(OK_results) {
  corrcoef <- cor(OK_results$pred, OK_results$obs_precip)
  print(paste("Correlation Coefficient Gauge Observations vs OK: ", corrcoef))
  
  outpng <- paste0("Output/OK_Scatter_", links_1day$date_time[1], ".png")
  title <- paste("OK vs Observed Precipitation ", datestr)
  png(outpng, width=800, height=600)
  plot(OK_results$pred~OK_results$obs_precip, main=title, 
       xlab="Observed", ylab="OK prediction", pch=16, col="blue")
  abline(lm(OK_results$pred~OK_results$obs_precip), col="red")
  legend("topleft", paste("Correlation: ", round(corrcoef,3)),
         cex=1.1, bty="n")
  dev.off()
}
  
  
load_radar <- function(datestr, grd) {
  # Load hourly radar (AAI format), crop to extent of grid and plot
  radar_file <- paste0(radar_dir,'rw_',gsub('-','',datestr), '.asc')
  print(paste("Loading RADOLAN file:", radar_file))
  radar <- raster(radar_file)
  radar <- projectRaster(radar, res=1000, crs=proj_laea)
  radar <- crop(radar, grd@bbox)
  
  # Plot radar
  outpng <- paste0("Output/Radar_", links_1day$date_time[1], ".png")
  title <- paste("Radar rain rate ", datestr)
  png(outpng, width=800, height=600)
  clrs <- brewer.pal(8,'YlGnBu')
  plot(radar, col=clrs, main=title)
  plot(links_1day, col="red", pch=18, cex=0.8, add=T)
  plot(gadm_1, border='black', lwd=1, add=T)
  dev.off()
  return(radar)   
}

perform_ked <- function(links_1day, gauges_1day, radar, grd, datestr) {
  # Kriging with External Drift
  # use radar values as external 
  # First: Remove duplicated locations!
  zdist <- zerodist(links_1day)[,1]
  if (length(zdist)>0) {
    # Remove rows identified by zerodist
    links_1day <- links_1day[-zdist,]
  }
  # Prepare new variogram and model with independant variables
  links_1day$radar_precip <- extract(radar, links_1day)
  vg2 <- variogram(rain_rate~radar_precip, data=links_1day)
  vg2_fit <- fit.variogram(vg2, 
                          vgm(psill=0.5,c('Exp','Sph','Gau'),range=20000))
  
  # Add to the grd object an attibute column 
  # for the independant variable "radar_precip"
  grd$radar_precip <- 0.0
  KED_grid <- krige(rain_rate ~ radar_precip, 
                    links_1day, grd, model=vg2_fit)
  # Plot KED 
  clrs <- brewer.pal(8,'BuPu')
  outpng <- paste0("Output/KED_interp_", links_1day$date_time[1], ".png")
  title <- paste("Kriging with External Drift Interpolation ", datestr)
  png(outpng, width=800, height=600)
  plot(KED_grid, col=clrs, main=title)
  plot(links_1day, col="red", pch=18, cex=1.2, add=T)
  plot(gadm_1, border='black', lwd=1, add=T)
  dev.off()
  
  # Kriging at gauge locations
  # Add attribute column "radar_precip"  to the gauges SPDF
  gauges_1day$radar_precip <- 0.0
  KED_gauges <- krige(rain_rate ~ radar_precip, 
                      links_1day, gauges_1day, model=vg_fit)
  
  # Join results with gauges to get the station_id, and observed precip 
  gauges_obs <- data.frame("station_id"=gauges_1day$station_id, 
                           "obs_precip"=gauges_1day$obs_precip,
                          "x_coord"=gauges_1day@coords[,1], 
                          "y_coord"=gauges_1day@coords[,2])
  links_pred <- data.frame("pred"=KED_gauges$var1.pred,
                           "x_coord"=KED_gauges@coords[,1], 
                           "y_coord"=KED_gauges@coords[,2])
  KED_results <- merge(gauges_obs, links_pred, by=c('x_coord', 'y_coord'))
  return(KED_results)
}

KED_scatter <- function(KED_results, datestr) {
  # Scatter plot of KED vs gauge observations
  corrcoef <- cor(KED_results$pred, KED_results$obs_precip)
  print(paste("Correlation Coefficient Gauge Observations vs KED: ", corrcoef))
  title <- paste("KED vs Observed Precipitation", datestr)
  outpng <- paste0("Output/KED_Scatter_", datestr, ".png")
  png(outpng, width=800, height=600)
  plot(KED_results$pred ~ KED_results$obs_precip, 
       main=title, 
       xlab="Observed", ylab="KED prediction", pch=16, col="blue")
  abline(lm(KED_results$pred ~ KED_results$obs_precip), col="red")
  legend("topleft", paste("Correlation: ", round(corrcoef,3)),
         cex=1.1, bty="n")
  dev.off()
}

gauge_radar_scatter <- function(gauges_1day, radar, datestr) {
  # Second scatterplot comparing gauge observations with radar values
  gauges_1day$radar_precip <- extract(radar, gauges_1day)
  corrcoef <- cor(gauges_1day$obs_precip, gauges_1day$radar_precip)
  print(paste("Correlation Coefficient Gauge Observations vs Radar: ", corrcoef))
  
  outpng <- paste0("Output/Gauge_Radar_Scatter_", datestr, ".png")
  png(outpng, width=800, height=600)
  title <- paste("Radar vs Observed Precipitation", datestr)
  plot(gauges_1day$obs_precip, gauges_1day$radar_precip, 
       main=title, 
       xlab="Observed", ylab="Radar precip", pch=16, col="blue")
  abline(lm(gauges_1day$obs_precip ~ gauges_1day$radar_precip), col="red")
  legend("topleft", paste("Correlation: ", round(corrcoef,3)), 
         cex=1.1, bty="n")
  dev.off() 
}