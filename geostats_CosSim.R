library('dplyr')
library('sp')
library('reshape2')

# ---------------------------------
# Load data and metadata
# ---------------------------------
setwd(paste0('/home/micha/Studies/Courses',
             '/Geostatistics-Tal/Project/'))

date_str_from <- '20170101'
date_str_to <- '20170131'
data_file <- paste0('gauge_data/','gauge_data_daily.csv')
data_cols <- c('station_id','date_time','quality',
                 'obs_precip', 'precip_ind','snow','eor')

gauge_data <- read.csv(data_file, col.names=data_cols)

meta_cols <- c('station_id','from_date','to_date','elevation',
               'latitude','longitude',
               'stn_name','province')
gauge_metadata <- read.csv('gauge_data/gauge_metadata_bavaria.csv', col.names=meta_cols)

# Get one day (or hour), and attach metadata, for plotting
gauges <- filter(gauge_data,
                 date_time == date_str_from)
# Data.frame for all of January
gauges_jan17 <- filter(gauge_data, 
                              date_time >= date_str_from & date_time <= date_str_to)

gauges <- merge(gauges, gauge_metadata, by='station_id', all.y=TRUE)
gauges_jan17 <- merge(gauges_jan17, gauge_metadata, by='station_id', all.y=TRUE)
# Make sure to clean out NA or < 0 (unknown values)
gauges_jan17 <- na.omit(gauges_jan17)
gauges_jan17 <- filter(gauges_jan17, obs_precip>=0)

# Create "Wide" matrix with one row for each day, 270 columns of precip data 
gauges_precip = as.data.frame(acast(gauges_jan17, date_time~station_id, value.var="obs_precip"))

# ---------------------------------
# Show gauges on a map
# ---------------------------------
library('ggplot2')
library('ggmap')

bavaria_map <- get_map(location = c(11.5,49), zoom=7)

# Select test gauge
# Subset out the two test gauges, to plot in red
test_gauge1 <- '1161'     # Center of region
test_gauge2 <- '4617'   # In the mountains
select_gauges <- subset(gauges, station_id %in% c(test_gauge1,test_gauge2))
ggmap(bavaria_map) + 
  geom_point(data=gauges, aes(x=longitude, y=latitude), 
             col="blue", size=1.5) +
  geom_point(data=select_gauges, aes(x=longitude, y=latitude),
             col="red", size=3) +
  geom_text(aes(label=station_id, x=longitude, y=latitude, vjust=1.5), 
            size=2, data=gauges)

# ---------------------------------
# Write function for cosine similarity index
# ---------------------------------
cosine_similarity <- function(x,y) {
  # x, y are vectors
  sum_xy <- sum(x*y)
  sum_x2 <- sum(x^2)
  sum_y2 <- sum(y^2)
  csi <- (sum_xy)/(sqrt(sum_x2)*sqrt(sum_y2))
  return(csi)
}

# ---------------------------------
# Loop thru gauges_precip 
# and calculate correlation coeff and cosine similarity
# ---------------------------------
corr_coef1 <- as.data.frame(cor(gauges_precip, gauges_precip[test_gauge1]))
corr_coef2 <- as.data.frame(cor(gauges_precip, gauges_precip[test_gauge2]))

cos_sim1 <- as.data.frame(matrix(ncol=1, nrow=dim(gauges_precip)[2]),
                       row.names = colnames(gauges_precip))
cos_sim2 <- as.data.frame(matrix(ncol=1, nrow=dim(gauges_precip)[2]), 
                       row.names = colnames(gauges_precip))

#Use rownames to add column for station_id
corr_coef1$station_id <- rownames(corr_coef1)
corr_coef2$station_id <- rownames(corr_coef2)
cos_sim1$station_id <- rownames(cos_sim1)
cos_sim2$station_id <- rownames(cos_sim2)
# Setup column names
colnames(corr_coef1) <- c('Rsquared1','station_id')
colnames(corr_coef2) <- c('Rsquared2','station_id')
colnames(cos_sim1) <- c('CosSim1','station_id')
colnames(cos_sim2) <- c('CosSim2','station_id')

for (i in colnames(gauges_precip)) {
  cos_sim1[i,]$CosSim1 <- cosine_similarity(gauges_precip[i], gauges_precip[test_gauge1])
  cos_sim2[i,]$CosSim2 <- cosine_similarity(gauges_precip[i], gauges_precip[test_gauge2])
}

# ---------------------------------
# Now plot correlation coefficients and cosine similarities
# ---------------------------------
library(gridExtra)
library(RColorBrewer)
clrs <- brewer.pal(name='RdYlBu', n=6)

gauges <- merge(gauges, corr_coef1, by='station_id', all.y=TRUE)
gauges <- merge(gauges, corr_coef2, by='station_id', all.y=TRUE)
gauges <- merge(gauges, cos_sim1, by='station_id', all.y=TRUE)
gauges <- merge(gauges, cos_sim2, by='station_id', all.y=TRUE)

#Pair of maps for first test point
png('Rsq_CosSim1.png', width=1200, height=800)
map1 <- ggmap(bavaria_map) + 
  geom_point(data=gauges, aes(x=longitude, y=latitude, colour=Rsquared1), 
              size=2) +
  scale_color_gradientn(colours = clrs)
map2 <- ggmap(bavaria_map) + 
  geom_point(data=gauges, aes(x=longitude, y=latitude, colour=CosSim1), 
             size=2) +
  scale_color_gradientn(colours = clrs)
grid.arrange(map1, map2, nrow=1)
dev.off()

# Pair of maps for second test point
png('Rsq_CosSim2.png', width=1200, height=800)
map3 <- ggmap(bavaria_map) + 
  geom_point(data=gauges, aes(x=longitude, y=latitude, colour=Rsquared2), 
             size=2) +
  scale_color_gradientn(colours = clrs)
map4 <- ggmap(bavaria_map) + 
  geom_point(data=gauges, aes(x=longitude, y=latitude, colour=CosSim2), 
             size=2) +
  scale_color_gradientn(colours = clrs)
grid.arrange(map3, map4, nrow=1)
dev.off()

# And scatter plot
ggplot(data=gauges) + geom_point(aes(x=Rsquared1, y=CosSim1), col='Blue')
