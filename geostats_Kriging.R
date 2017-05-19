library('dplyr')
library('spatstat')
library('spdep')
library('geoR')
library('raster')
library('RColorBrewer')
#library('gstat')

# Setup
#----------------------------------------
setwd(paste0('/home/micha/Studies/Courses',
             '/Geostatistics-Tal/Project/'))
date_str <- '20170112'
data_file <- paste0('gauge_data/','gauge_data_daily.csv')
data_cols <- c('station_id','date_time','quality',
               'obs_precip', 'precip_ind','snow','eor')

gauge_data <- read.csv(data_file, col.names=data_cols)
meta_cols <- c('station_id','from_date','to_date','elevation',
               'latitude','longitude',
               'stn_name','province')
gauge_metadata <- read.csv('gauge_data/gauge_metadata.csv', col.names=meta_cols)
# Get one day (or hour), and attach metadata
gauge_data_filtered <- filter(gauge_data, date_time==date_str)
gauges <- merge(gauge_data_filtered, gauge_metadata, by='station_id', all.y=TRUE)
# Make sure to clean out NA or < 0 (unknown values)
gauges <- na.omit(gauges)
gauges <- filter(gauges, obs_precip>=0)

gadm_1 <- readRDS('GIS/DEU_adm1.rds')

# Prepare SpatialPointsDataFrame
# ---------------------------------
xy_coords <- cbind(gauges$longitude, gauges$latitude) 
proj_wgs84 = CRS("+init=epsg:4326")
coordinates(gauges) <- xy_coords
proj4string(gauges) <- proj_wgs84
# Now gauges should be a SpatialPointsDataFrame
class(gauges)

# IDW Interpolation
# ---------------------------------
# Convert to ppp
gauges_owin <- ripras(gauges$longitude, gauges$latitude, shape="rectangle")
gauges_marks <- gauges$obs_precip
gauges_ppp <- ppp(gauges$longitude, gauges$latitude, 
                  window=gauges_owin,
                  marks=gauges_marks)
#class(gauges_ppp)

# Set dimensions
minx <- min(gauges$longitude)
maxx <- max(gauges$longitude)
miny <- min(gauges$latitude)
maxy <- max(gauges$latitude)
nsdim <- ceiling((maxy - miny)/0.01)
ewdim <- ceiling((maxx - minx)/0.01)

gauges_idw <- idw(gauges_ppp, dimyx=c(nsdim, ewdim))

clrs <- brewer.pal(8,'BuPu')
opar <- par(pch=3, col="grey30", cex=0.5, lwd=1)
plot(gauges_idw, col=clrs)
plot(gadm_1, border='black', add=T)
plot(gauges, add=T)

# Get IDW values at original points for comparison
gauges_idw_pts <- extract(raster(gauges_idw), gauges)
rmse_idw <- sqrt(mean((gauges$obs_precip - gauges_idw_pts)^2))
print(rmse_idw)

# Variogram and Envelope
#----------------------------------------
coord_matrix <- cbind(gauges$longitude, gauges$latitude)
vg <- variog(coords=coord_matrix, 
             data=gauges$obs_precip, max.dist=3)
# Run 9999 Monte Carlo simulations
vg_envel <- variog.mc.env(coords=coord_matrix, 
                          data=gauges$obs_precip, 
                          obj.var=vg, 
                          nsim=999)

options(repr.plot.width=5, repr.plot.height=5)
plot(vg, envelope=vg_envel, pch=16, col="blue")
# Guess at nugget, sill and range
krige_params = c(11, 1.5)
vg_fit_exp <- variofit(vg, ini.cov.pars=krige_params, 
                       cov.model = "exponential",
                       fix.nugget = F, nugget = 3)
vg_fit_sph <- variofit(vg, ini.cov.pars=krige_params, 
                       cov.model = "spherical",
                       fix.nugget = F, nugget = 3)
vg_fit_gau <- variofit(vg, ini.cov.pars=krige_params, 
                       cov.model = "gaussian",
                       fix.nugget = F, nugget = 3)
lines.variomodel(vg_fit_exp, col="red")
lines.variomodel(vg_fit_sph, col="orange")
lines.variomodel(vg_fit_gau, col="magenta")

grd <- expand.grid(seq(minx,maxx,0.01), seq(miny,maxy,0.01))
gauges_ok <- krige.conv(coords=coord_matrix, 
                        data=gauges$obs_precip, 
                        krige = krige.control(cov.pars=krige_params), 
                        locations = grd)
# Plot krige raster
image(gauges_ok, col=clrs)
plot(gadm_1, fill=F, border='black', lwd=1, add=T)

# Check RMSE with original points
gauges_ok_rast <- raster(xmn=minx, xmx=maxx, 
                         ymn=miny, ymx=maxy, 
                         nrows=nsdim, ncols=ewdim, 
                         vals=gauges_ok$predict)
gauges_ok_pts <- extract(gauges_ok_rast, gauges)
rmse_ok <- sqrt(mean((gauges$obs_precip - gauges_ok_pts)^2))
print(rmse_ok)