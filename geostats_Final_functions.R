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
plot_map <- function(cities, gadm_1, link_data) {
  # Clip gauges to region covered by links BBOX
  # first create spatial points data frame of links
 
  # Plot links and gauges
  outpng <- paste0(out_dir,'links_map.png')
  png(outpng, width=800, height=600)
  orig_par <- par()
  plot(link_data, col="red", main="Bavaria - CML Links",
       pch=18, cex=1.6)
  plot(gadm_1, border='black', add=T)
  points(cities, pch=16, col='grey30', cex=1.2)
  pointLabel(coordinates(cities), cities$asciiname, cex=0.9)
  lclrs <- c('red','grey30')
  lgnd <- c('CML Links','Cities')
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
slice_data <- function(link_data, datestr) {
  # Slice out link data for one day
  links_1day <- link_data[link_data$date_time == datestr,]
  # Remove all NaN values
  links_1day <- na.omit(links_1day)

  return(links_1day)
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
  outpng <- paste0(out_dir, 'Morans_histogram.png')
  png(outpng, width=800, height=800)
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
  outpng <- paste0(out_dir,"OK_variogram_", datestr, ".png")
  png(outpng, width=600, height=400)
  y_upper <- max(vg$gamma)*1.1
  x_upper <- max(vg$dist)*1.1
  title <- paste("OK Variogram_",datestr)
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
  # Also create set of 400 validation locations
  mnx <- min(link_data@coords[,1])-6000
  mxx <- max(link_data@coords[,1])+6000
  mny <- min(link_data@coords[,2])-6000
  mxy <- max(link_data@coords[,2])+6000
  grd <- expand.grid(lon=seq(mnx,mxx,1000), lat=seq(mny,mxy,1000))
  xy_coords <- cbind(grd$lon, grd$lat)
  coordinates(grd) <- xy_coords
  proj4string(grd) <- proj_laea
  gridded(grd) <- TRUE
  
  # And choose 400 random validation points
  validation_points <- spsample(grd, 400, 'random')
  return(c("grid"=grd, "validation"=validation_points))
}

load_radar <- function(datestr, grd) {
  # Load hourly radar (AAI format), crop to extent of grid and plot
  radar_file <- paste0(radar_dir,'rw_',gsub('-','',datestr), '.asc')
  print(paste("Loading RADOLAN file:", radar_file))
  radar <- raster(radar_file)
  radar <- projectRaster(radar, res=1000, crs=proj_laea)
  radar <- crop(radar, grd@bbox)
  
  # Plot radar
  outpng <- paste0(out_dir, "Radar_", links_1day$date_time[1], ".png")
  title <- paste("Radar rain rate ", datestr)
  png(outpng, width=800, height=600)
  clrs <- brewer.pal(8,'YlGnBu')
  plot(radar, col=clrs, main=title)
  points(links_1day, col="red", pch=18, cex=0.8)
  plot(gadm_1, border='black', lwd=1, add=T)
  dev.off()
  
  return(radar)   
}


perform_ok <- function(links_1day, vg_fit, grd, 
                       validation_points, radar, datestr) {
  # Function to preform Ordinary Kriging
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
  outpng <- paste0(out_dir, "OK_interp_", links_1day$date_time[1], ".png")
  png(outpng, width=800, height=600)
  orig_par <- par()
  
  plot(OK_grid, col=clrs, main=title)
  plot(links_1day, col="red", pch=18, cex=1.2, add=T)
  plot(gadm_1, border='black', lwd=1, add=T)
  par(orig_par)
  dev.off()
  
  # Kriging at validation locations
  OK_validation <- krige(rain_rate~1, 
                     locations=links_1day, 
                     newdata=validation_points, model=vg_fit)
  
  # Extract radar values at validation points and merge with OK result
  validation_precip <- extract(radar, validation_points)
  OK_results <- cbind(OK_validation$var1.pred, validation_precip)
  colnames(OK_results) <- c('OK_pred', 'radar_precip')
  OK_results <- as.data.frame(OK_results)
  return(OK_results)
}


OK_scatter <- function(OK_results) {
  corrcoef <- cor(OK_results$OK_pred, 
                  OK_results$radar_precip, use='complete')
  print(paste("Correlation Coefficient Radar vs OK: ", corrcoef))

  outpng <- paste0(out_dir, "OK_Scatter_", links_1day$date_time[1], ".png")
  title <- paste("OK vs Radar Precipitation ", datestr)
  png(outpng, width=800, height=600)
  plot(OK_results$OK_pred~OK_results$radar_precip, main=title, 
       xlab="Radar precip", ylab="OK prediction", pch=16, col="blue")
  abline(lm(OK_results$OK_pred~OK_results$radar_precip), col="red")
  legend("topleft", paste("Correlation: ", round(corrcoef,3)),
         cex=1.1, bty="n")
  dev.off()
}
  

perform_ked <- function(links_1day, grd, radar, 
                        validation_points, datestr) {
  # Kriging with External Drift
  # use radar values as external 
  # First: Remove duplicated locations!
  zdist <- zerodist(links_1day)[,1]
  if (length(zdist)>0) {
    # Remove rows identified by zerodist
    links_1day <- links_1day[-zdist,]
  }
  # Prepare new variogram and model with independant variables
  links_1day$radar_values <- extract(radar, links_1day)
  
  vg2 <- variogram(rain_rate~radar_values, data=links_1day)
  
  vg2_fit <- fit.variogram(vg2, 
                          vgm(psill=50,c('Exp','Sph','Gau'),range=100000))
  x_upper <- max(vg2$dist)*1.1
  y_upper <- max(vg2$gamma)*1.1
  title <- paste("KED Variogram", datestr)
  outpng <- paste0(out_dir,"KED_variogram_", datestr, ".png")
  png(outpng, width=600, height=400)
  
  plot(gamma~dist, vg2, col="blue", main=title,
       xlab="Distance", ylab="Semivariance", 
       ylim=c(0, y_upper), xlim=c(0,x_upper))
  lines(variogramLine(vg2_fit, x_upper), col='red')
  dev.off()
  print(vg2_fit)
  
  # Add to the grd object an attibute column 
  # for the independant variable "radar_values"
  grd$radar_values <- 0.0
  KED_grid <- krige(rain_rate ~ radar_values, 
                    locations=links_1day, 
                    newdata=grd, model=vg2_fit)
  # Plot KED 
  clrs <- brewer.pal(8,'BuPu')
  outpng <- paste0(out_dir, "KED_interp_", links_1day$date_time[1], ".png")
  title <- paste("Kriging with External Drift Interpolation ", datestr)
  png(outpng, width=800, height=600)
  plot(KED_grid, col=clrs, main=title)
  plot(links_1day, col="red", pch=18, cex=1.2, add=T)
  plot(gadm_1, border='black', lwd=1, add=T)
  dev.off()
  
  # Kriging at validation locations
  # First prepare an attrib column in validation points 
  # to accept the radar values
  radar_values <- rep(0.0, length(validation_points))
  validation_points <- SpatialPointsDataFrame(validation_points, 
                                              as.data.frame(radar_values))
  # Now KED at the validation locations
  KED_validation <- krige(rain_rate ~ radar_values, 
                      locations=links_1day, 
                      newdata=validation_points, model=vg2_fit)

  # Extract actual radar values at validation points and merge with OK result
  validation_precip <- extract(radar, validation_points)
  KED_results <- cbind(KED_validation$var1.pred, validation_precip)
  colnames(KED_results) <- c('KED_pred', 'radar_precip')
  KED_results <- as.data.frame(KED_results)
  pred_radar_correlation <- cor(KED_results$KED_pred, 
                                KED_results$radar_precip,
                                use='complete')
  print(pred_radar_correlation)
  return(KED_results)
}

KED_scatter <- function(KED_results, datestr) {
  # Scatter plot of KED vs gauge observations
  corrcoef <- cor(KED_results$KED_pred, 
                  KED_results$radar_precip, 
                  use="complete")
  print(paste("Correlation Coefficient Radar vs KED: ", corrcoef))
  title <- paste("KED vs Radar Precipitation", datestr)
  outpng <- paste0(out_dir, "KED_Scatter_", datestr, ".png")
  png(outpng, width=800, height=600)
  plot(KED_results$KED_pred ~ KED_results$radar_precip, 
       main=title, 
       xlab="Observed", ylab="KED prediction", pch=16, col="blue")
  abline(lm(KED_results$KED_pred ~ KED_results$radar_precip), col="red")
  legend("topleft", paste("Correlation: ", round(corrcoef,3)),
         cex=1.1, bty="n")
  dev.off()
}

check_KED <- function(links_1day, radar) {
 
  links_1day$radar_values <- extract(radar, links_1day)
  # Check correlation of radar (independant) to links
  # Pearson NOT suitable since the distributions of rainrall are not Gaussian
  par(mfrow=c(1,2))
  outpng <- paste0(out_dir,"Histograms_", links_1day$date_time[1],".png")
#  png(outpng, width=800, height=600)
  hist(links_1day$radar_values, main="Histogram of Radar Precipitation", 
       breaks=30, col="green", xlab="radar_precip")  
  hist(links_1day$rain_rate, main="Histogram of Link Rain Rates", 
       breaks=30, col="blue", xlab="link rain_rates") 
#  dev.off()
  
  #link_radar_pearson <- cor.test(links_1day$rain_rate, 
  #                                  links_1day$radar_values, 
  #                                  use="complet", alternative = 't')
  
  link_radar_spearman <- cor.test(links_1day$rain_rate, 
                                  links_1day$radar_values,
                                  method='spearman', exact=FALSE)
  print(link_radar_spearman)
  
}