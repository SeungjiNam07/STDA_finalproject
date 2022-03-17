#### STDA Final Project
#### Seung ji Nam

# load library
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
library(rgdal)
library(nlme)
library(batchmeans)
library(nimble);library(mvtnorm)
library(spacetime)

# set your directory
directory = '/Users/seungji/Library/Mobile Documents/com~apple~CloudDocs/STDA/Project/'
sensor = readOGR(paste0(directory,'processeddata/censor.shp'))
sensor_oct = readOGR(paste0(directory, 'processeddata/censor_oct.shp'))
sensor_march = readOGR(paste0(directory, 'processeddata/censor_march.shp'))
sensor_mt = readOGR(paste0(directory,'processeddata/sensor_mt.shp'))
sensor_wk = readOGR(paste0(directory,'processeddata/sensor_wk.shp'))
################################################################################

## 1. Spatial correlation
# 1) linear model 

# October
# ols1 = lm(fp~ bus+subway, data = sensor_oct)
ols1 = lm(fp~ bus+subway+building_c+ building_a+pop, data = sensor_oct)
summary(ols1)
fitted1 = predict(ols1, sensor_oct[,c('bus','subway','building_c','building_a','pop')])
ehat1 = sensor_oct$fp-fitted1
ols1$coefficients

# March
ols2 = lm(fp~ bus+subway+building_c+ building_a+pop, data = sensor_march)
summary(ols2)
fitted2 = predict(ols2, sensor_march[,c('bus','subway','building_c','building_a','pop')])
ehat2 = sensor_march$fp-fitted2
ols2$coefficients


# draw plot
plot.point.ref <- function ( spatialdata , vals ) {
  
  pal <- tim.colors (10)
  ints <- classIntervals (vals , n = 8, style = "pretty") # Determine breakpoints
  # also see style options " quantile " and " fisher "
  intcols <- findColours (ints , pal ) # vector of colors
  # if pal doesn 't have the same length as # classes , findColours will interpolate
  
  par ( mar = rep (0, 4))
  plot ( spatialdata , col = intcols , pch = 19)
  points ( spatialdata , pch = 1)
  legend ("topleft", fill = attr ( intcols , "palette"),
          legend = names ( attr ( intcols , "table")), bty = "n")
}
plot.point.ref(sensor_oct, ehat1)
plot.point.ref(sensor_march, ehat1)

# 2) variogram
# October
sensor_oct$ehat = ehat1
vg1 <- variogram(ehat1~1, data = sensor_oct, width = 0.5)
par(mfrow=c(1,2))
# plot(vg, xlab = 'Distance',ylab = 'semi-variogram estimate', width =0.5)

ggplot(data = vg1)+geom_point(aes(x = dist,y = gamma))+
  labs( x = 'Distance',y = 'semi-variogram estimate')+
  coord_cartesian( ylim = c(0,max(vg1$gamma)+10000))+
  theme_bw()+theme_light()+ theme(legend.position = "none")

vgangle1 <- variogram(ehat ~ 1, data = sensor_oct, width = 1,alpha = c(0, 45, 90, 135))
# plot(vgangle, xlab = "Distance", ylab = "Semi-variogram estimate")
ggplot(data = vgangle1)+geom_point(aes(x = dist,y = gamma))+
  labs( x = 'Distance',y = 'semi-variogram estimate')+
  coord_cartesian( ylim = c(0,max(vgangle1$gamma)+10000))+
  facet_wrap(dir.hor~.)+ scale_x_discrete(limits=seq(0,12,2))+
  theme_bw()+theme_light()+ theme(legend.position = "none")

# March
sensor_march$ehat = ehat2
vg2 <- variogram(ehat2~1, data = sensor_march, width = 1)
par(mfrow=c(1,2))
# plot(vg, xlab = 'Distance',ylab = 'semi-variogram estimate', width =0.5)

ggplot(data = vg2)+geom_point(aes(x = dist,y = gamma))+
  labs( x = 'Distance',y = 'semi-variogram estimate')+
  coord_cartesian( ylim = c(0,max(vg2$gamma)+10000))+
  theme_bw()+theme_light()+ theme(legend.position = "none")

vgangle2 <- variogram(ehat ~ 1, data = sensor_march, width = 1,alpha = c(0, 45, 90, 135))
# plot(vgangle, xlab = "Distance", ylab = "Semi-variogram estimate")
ggplot(data = vgangle2)+geom_point(aes(x = dist,y = gamma))+
  labs( x = 'Distance',y = 'semi-variogram estimate')+
  coord_cartesian( ylim = c(0,max(vgangle2$gamma)+10000))+
  facet_wrap(dir.hor~.)+ scale_x_discrete(limits=seq(0,12,2))+
  theme_bw()+theme_light()+ theme(legend.position = "none")


# fit covariance function
vg1 <- variogram(ehat~1, data = sensor_oct, width = 0.5)
fitvg1 = fit.variogram(vg1, vgm(1,'Exp',12, 0.01))
# fitvg = fit.variogram(vg, vgm(1,'Mat',12, 0.01))
print(fitvg1)

tau2.hat = fitvg1$psill[1]
s2.hat = fitvg1$psill[2]
rho.hat = fitvg1$range[2]

# 3) gls
x = sensor_oct$lat
y = sensor_oct$lon

gls.fit <- gls(fp~ bus+subway+building_c+ building_a+pop, data = sensor_oct,
               corSpher(value = c(range = rho.hat, nugget = tau2.hat/(tau2.hat+s2.hat)),
                        nugget = TRUE, form=~x+y, fixed = TRUE))
summary(gls.fit)

# March
vg2 <- variogram(ehat~1, data = sensor_march, width = 0.5)
fitvg2 = fit.variogram(vg2, vgm(1,'Exp',12, 0.01))
# fitvg = fit.variogram(vg, vgm(1,'Mat',12, 0.01))
print(fitvg2)

tau2.hat = fitvg2$psill[1]
s2.hat = fitvg2$psill[2]
rho.hat = fitvg2$range[2]

# 3) gls
x = sensor_march$lat
y = sensor_march$lon

gls.fit <- gls(fp~ bus+subway+building_c+ building_a+pop, data = sensor_march,
               corSpher(value = c(range = rho.hat, nugget = tau2.hat/(tau2.hat+s2.hat)),
                        nugget = TRUE, form=~x+y, fixed = TRUE))
summary(gls.fit)

################################################################################


# 2. Temporal correlation
# * Daily
# 1) variogram 
sensor@data$date = as.Date(sensor@data$date)
sensor@data = sensor@data[sensor@data$date>='2021-01-07',]
sensor@data <- sensor@data[order(sensor@data$date),]
sensor@data$building_c = as.numeric(sensor@data$building_c)
sensor@data = sensor@data[, -c(4,5)]
lat = sensor@data$lat
lon = sensor@data$lon

sp = SpatialPoints(unique(sensor@coords))
sensor_st = STFDF(sp, unique(sensor@data$date), sensor@data[,c(3,4,5,6,7,8)])
tm=unique(sensor_st@time)

# variogram

# inclue only subway and bus 
vv1 = variogramST(fp ~1 + subway+bus,
                     data = sensor_st, tlags = 0.01:7.01)
sepVgm <- vgmST(stModel = 'separable',
                space = vgm(10,'Exp', 12, nugget = 0.5),
                time = vgm(10, 'Exp',75, nugget = 0.1),
                sill = 20)
sepVgm1 <- fit.StVariogram(vv1, sepVgm)

metricVgm <- vgmST(stModel = "metric",
                   joint = vgm(100, "Exp", 400, nugget = 0.1),
                   sill = 10,
                   stAni = 100)
metricVgm1 <- fit.StVariogram(vv1, metricVgm)

metricMSE <- attr(metricVgm, "optim")$value
sepMSE <- attr(sepVgm, "optim")$value
plot(vv1, list(sepVgm1,metricVgm1))


# 2) ACF plot
df = data.frame(sensor@data)
serial = sample(unique(df$serial),8,replace =F)
par(mfrow=c(4,2),mar=c(2,2,2,2))
for (s in serial){
  df1 = df[df$serial ==s,]
  # acf(df1$fp, main = 'floating population')
  model = lm(fp~bus+subway, data = df1)
  acf(fitted(model)-df1$fp, main = 'residuals')
}


## * monthly data
sensor_mt = readOGR(paste0(directory,'processeddata/sensor_mt.shp'))
sensor_mt@data$month = as.Date(as.character(sensor_mt@data$month))
# unique(sensor_mt@data$month)
sensor_mt@data <- sensor_mt@data[order(sensor_mt@data$month),]
sensor_mt@data$building_c = as.numeric(sensor_mt@data$building_c)
sensor_mt@data = sensor_mt@data[, -c(4,5)]
sp = SpatialPoints(unique(sensor_mt@coords))
sensor_mtst = STFDF(sp, unique(sensor_mt@data$month), sensor_mt@data[,c(3,4,5,6,7,8)])


# 1) Variogram
# include all variables
vv_all = variogramST(fp ~1 + subway+bus,
                     data = sensor_mtst, tlags = 0.01:1.01)
sepVgm_all <- vgmST(stModel = 'separable',
                    space = vgm(10,'Mat', 12, nugget = 0.5),
                    time = vgm(10, 'Mat',100, nugget = 0.5),
                    sill = 20)
sepVgm_all <- fit.StVariogram(vv_all, sepVgm_all)

metricVgm_all <- vgmST(stModel = "metric",
                       joint = vgm(100, "Exp", 400, nugget = 0.1),
                       sill = 10,
                       stAni = 100)
metricVgm_all <- fit.StVariogram(vv_all, metricVgm_all)
plot(vv_all, list(sepVgm_all,metricVgm_all))


# 2) acf
df = data.frame(sensor_mt@data)
serial = sample(unique(df$serial),8,replace =F)
par(mfrow=c(4,2),mar=c(2,2,2,2))
for (s in serial){
  df1 = df[df$serial ==s,]
  plot(acf,  main = 'floating population')
  model = lm(fp~bus+subway, data = df1)
  acf(fitted(model)-df1$fp, main = 'residuals')
}

## * weekly data
sensor_wk = readOGR(paste0(directory,'processeddata/sensor_wk.shp'))
sensor_wk@data$week = as.Date(as.character(sensor_wk@data$week2))
# unique(sensor_mt@data$month)
sensor_wk@data <- sensor_wk@data[order(sensor_wk@data$week),]
sensor_wk@data$building_c = as.numeric(sensor_wk@data$building_c)
sensor_wk@data = sensor_wk@data[, -c(4,5)]
sp = SpatialPoints(unique(sensor_wk@coords))
sensor_wkst = STFDF(sp, unique(sensor_wk@data$week), sensor_wk@data[,c(3,4,5,6,7,8)])


# 1) Variogram
# include all variables
vv_all = variogramST(fp ~1 + subway+bus,
                     data = sensor_wkst, tlags = 0.01:4.01)
sepVgm_all <- vgmST(stModel = 'separable',
                    space = vgm(10,'Exp', 12, nugget = 0.5),
                    time = vgm(10, 'Exp',10, nugget = 0.5),
                    sill = 20)
sepVgm_all <- fit.StVariogram(vv_all, sepVgm_all)

metricVgm_all <- vgmST(stModel = "metric",
                       joint = vgm(100, "Exp", 400, nugget = 0.1),
                       sill = 10,
                       stAni = 100)
metricVgm_all <- fit.StVariogram(vv_all, metricVgm_all)
plot(vv_all, list(sepVgm_all,metricVgm_all))


# 2) acf
df = data.frame(sensor_mt@data)
serial = sample(unique(df$serial),8,replace =F)
par(mfrow=c(4,2),mar=c(2,2,2,2))
for (s in serial){
  df1 = df[df$serial ==s,]
  acf(df1$fp, xlab = 'floating population')
  model = lm(fp~bus+subway, data = df1)
  acf(fitted(model)-df1$fp, main = 'residuals')
}

dev.off()

