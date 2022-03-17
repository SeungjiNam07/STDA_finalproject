## spatio-temporal data
library(sf)
library(classInt)
library(maps)
library(sp)
library(classInt)
library(gstat)
library(fields)
library(stringr)
library(mapdata)
library(mapproj)
library(batchmeans)

library(rgdal)
library(nimble);library(mvtnorm)
library(spacetime)
load('/Users/seungji/Library/Mobile Documents/com~apple~CloudDocs/STDA/Project/modeling/working.RData')
imac = '/Users/namseungji/Library/Mobile Documents/com~apple~CloudDocs/STDA/Project/'
macair = '/Users/seungji/Library/Mobile Documents/com~apple~CloudDocs/STDA/Project/'
# censor = st_read('/Users/namseungji/Library/Mobile Documents/com~apple~CloudDocs/STDA/Project/data/censor.shp')
sensor = readOGR(paste0(macair,'/data/censor.shp'))
seoul = readOGR(paste0(macair,'data/LARD_ADM_SECT_SGG_서울/LARD_ADM_SECT_SGG_11.shp'))
sensor@data$date = as.Date(sensor@data$date)
sensor@data = sensor@data[sensor@data$date>='2021-01-07',]
sensor@data <- sensor@data[order(sensor@data$date),]
sensor@data$building_c = as.numeric(sensor@data$building_c)
sensor@data = sensor@data[, -c(4,5)]
lat = sensor@data$lat
lon = sensor@data$lon
head(sensor@data)

sp = SpatialPoints(unique(sensor@coords))

# library(xts)
# library(PerformanceAnalytics)
# sensor2 = xts(sensor@data[,-2], sensor@data$date)
# sensor2 = PerformanceAnalytics::Return.calculate(sensor2)
# sensor@data = sensor2
# sensor_st = STFDF(sp, unique(sensor@data$date), sensor@data[,c(3,8,9,10,11,12)])
sensor_st = STFDF(sp, unique(sensor@data$date), sensor@data[,c(3,4)])
tm=unique(sensor_st@time)
# subway+bus+building_c+building_a+pop,
vv = variogramST(fp ~1 + subway,
          data = sensor_st,
          # width = 0.5,
          # cutoff = 0.01,
          tlags = 0.01:7.01)
sepVgm <- vgmST(stModel = 'separable',
                space = vgm(10,'Mat', 12, nugget = 0.5),
                time = vgm(10, 'Mat',100, nugget = 0.5),
                sill = 20)
sepVgm <- fit.StVariogram(vv, sepVgm)

metricVgm <- vgmST(stModel = "metric",
                   joint = vgm(100, "Exp", 400, nugget = 0.1),
                   sill = 10,
                   stAni = 100)
metricVgm <- fit.StVariogram(vv, metricVgm)

metricMSE <- attr(metricVgm, "optim")$value
sepMSE <- attr(sepVgm, "optim")$value

# plot(vv, list(sepVgm, metricVgm), main = "Semi-variance")
plot(vv, sepVgm, main = "Semi-variance")
plot(vv, metricVgm, main = "Semi-variance")

############################################################################
# sensor_st = STFDF(sp, unique(sensor@data$date), sensor@data[,c(3,8,9,10,11,12)])
sensor_st = STFDF(sp, unique(sensor@data$date), sensor@data[,c(3,6,7,8,9,10)])
sensor_df = data.frame(sensor_st)

col = c('subway','bus','building_c','building_a','pop')
for (c in col){
  sensor_st@data[,c] =sensor_df[,c]/(max(sensor_df[,c])-min(sensor_df[,c]))
  
}
tm=unique(sensor_st@time)

# subway+bus+building_c+building_a+pop,
vv = variogramST(fp ~1 + subway+bus,
                 data = sensor_st,
                 # width = 0.5,
                 # cutoff = 0.01,
                 tlags = 0.01:7.01)
sepVgm <- vgmST(stModel = 'separable',
                space = vgm(10,'Exp', 12, nugget = 0.5),
                time = vgm(10, 'Exp',75, nugget = 0.1),
                sill = 20)
sepVgm <- fit.StVariogram(vv, sepVgm)


# vv = variogramST(fp ~1 + subway+bus,
#                  data = sensor3,
#                  # width = 0.5,
#                  # cutoff = 0.01,
#                  tlags = 0.01:7.01)
sepVgm <- vgmST(stModel = 'separable',
                space = vgm(10,'Exp', 12, nugget = 0.5),
                time = vgm(10, 'Exp',75, nugget = 0.1),
                sill = 20)
sepVgm <- fit.StVariogram(vv, sepVgm)

metricVgm <- vgmST(stModel = "metric",
                   joint = vgm(100, "Exp", 400, nugget = 0.1),
                   sill = 10,
                   stAni = 100)
metricVgm <- fit.StVariogram(vv, metricVgm)

metricMSE <- attr(metricVgm, "optim")$value
sepMSE <- attr(sepVgm, "optim")$value

# plot(vv, list(sepVgm, metricVgm), main = "Semi-variance")
plot(vv, sepVgm, main = "Semi-variance")
plot(vv, metricVgm, main = "Semi-variance")
############################################################################
############################################################################
############################################################################

# spatial grid
min_lat = min(lat); max_lat = max(lat)
min_lon = min(lon); max_lon = max(lon)

spat_pred_grid <- expand.grid(
  lon = seq(min_lon-0.01, max_lon+0.01, length = 2),
  lat = seq(min_lat-0.01, max_lat+0.01, length = 2)) %>%
  SpatialPoints(proj4string = CRS(proj4string(sensor_st)))
gridded(spat_pred_grid) <- TRUE

# temporal grid 
temp_pred_grid <- as.Date("2021-01-07") + seq(0, 210, length = 2)

DE_pred <- STF(sp = spat_pred_grid, # spatial part
               time = temp_pred_grid) # temporal part

sensor_st <- as(sensor_st, "STIDF") # convert to STIDF
# sensor3 <- subset(sensor3, !is.na(sensor3$fp)) # remove missing data

pred_kriged <- krigeST(fp ~ 1 , # latitude trend
                       data = sensor_st, # data set w/o 14 July
                       newdata = DE_pred, # prediction grid
                       modelList = sepVgm, # semivariogram
                       computeVar = TRUE) # compute variances

color_pal <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(16))

stplot(pred_kriged,
       main = "Predictions (degrees Fahrenheit)",
       layout = c(3, 2),
       col.regions = color_pal)





############################################################################3
install_github("andrewzm/STRbook")
library(STRbook)
data("STObj3", package = "STRbook")

