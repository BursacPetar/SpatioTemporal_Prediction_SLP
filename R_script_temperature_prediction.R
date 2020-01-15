library(rmarkdown)
library(knitr)
library(kableExtra)
library(sp)
library(spacetime)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mapview)
library(gstat)
library(RColorBrewer)
library(rworldmap)
library(zoo)
library(xts)
library(gridExtra)
library(ggsflabel)
library(stringr)
library(magrittr)
library(STRbook)
library(ggpubr)
library(maps)
library(ggspatial)
library(raster)
library(stars)
library(viridis)
library(doParallel)
library(leafsync)
library(meteo)
library(plotGoogleMaps)
library(lubridate)

load("d:/Doktorske_akademske_studije_Geodezija_i_geoinformatika/Semestar_2/Prostorno_vremenska_statisitika/Temperature [Sale]/folder/ogimet_temp_stfdf.rda")
stfdf_temp

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Vremenska serija podataka
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

df_temp <- as.data.frame(stfdf_temp)

# Remove NAs
df_temp <- df_temp[!is.na(df_temp$tmax), ]
df_temp <- df_temp[!is.na(df_temp$tmin), ]
df_temp <- df_temp[!is.na(df_temp$tmean), ]


df_temp$time <- str_sub(df_temp$time, end=-1)
obs_temp <- df_temp
sta_temp <- as.data.frame(stfdf_temp@sp)

obs_temp$date_1 <- obs_temp$time
obs_temp %<>% separate(date_1, "-", into = c("y", "m", "d"))

tempv <- obs_temp
tempv %<>% mutate(y = as.numeric(y),
                 m = as.numeric(m),
                 d = as.numeric(d))

toJulian <- function(y, m, d){
  mm<- m
  xx<- 0
  if( mm<=2) {xx<- 1}
  mm<- (12*xx)+mm
  yy<- y-xx
  nc<- floor(0.01*yy)
  jd<- floor(365.25*yy)+floor(30.6001*(1+mm))+d+1720995+(2-(nc-floor(0.25*nc)))
  return(jd)
}

tempv$julian <- toJulian(tempv$y, tempv$m, tempv$d)
tempv$t <- tempv$julian - 2454832

# TMEAN
tmean_dan <- tempv %>% 
  dplyr::group_by(y, m, d) %>%
  summarise(tmean_dan = mean(tmean), 
            julian = mean(julian))
tmean_dan$vreme <- paste(tmean_dan$y, "-", tmean_dan$m, "-", tmean_dan$d)

# TMAX
tmax_dan <- tempv %>% 
  dplyr::group_by(y, m, d) %>%
  summarise(tmax_dan = mean(tmax), 
            julian = mean(julian))
tmax_dan$vreme <- paste(tmax_dan$y, "-", tmax_dan$m, "-", tmax_dan$d)

# TMIN
tmin_dan <- tempv %>% 
  dplyr::group_by(y, m, d) %>%
  summarise(tmin_dan = mean(tmin), 
            julian = mean(julian))
tmin_dan$vreme <- paste(tmin_dan$y, "-", tmin_dan$m, "-", tmin_dan$d)

j_temp <- ggplot()+
  geom_line(data = tmean_dan, aes(x = julian, y = tmean_dan, group = y, colour = "tmean"))+
  geom_line(data = tmax_dan, aes(x = julian, y = tmax_dan, group = y, colour = "tmax"))+
  geom_line(data = tmin_dan, aes(x = julian, y = tmin_dan, group = y, colour = "tmin"))+
  xlab("Julian [2455000 = 2010 godina]") + 
  ylab("Temperatura [°C]") +
  ggtitle("Vremenska linija opazanja Temperature vazduha [°C]")+
  theme_bw()

tm_y <- ggplot()+
  geom_boxplot(data = tmean_dan, aes(x = y, y = tmean_dan, group = y), colour = "#00BA38")+
    xlab("Julian [2455000 = 2010 godina]") + 
  ylab("Temperatura [°C]") +
  ggtitle("Prosecna temperatura po godinama [°C]")+
  theme_bw()
tmax_y <- ggplot()+
  geom_boxplot(data = tmax_dan, aes(x = y, y = tmax_dan, group = y), colour = "#F8766D")+
  xlab("Julian [2455000 = 2010 godina]") + 
  ylab("Temperatura [°C]") +
  ggtitle("Maksimalne vrednosti temperature po godinama [°C]")+
  theme_bw()
tmin_y <- ggplot()+
  geom_boxplot(data = tmin_dan, aes(x = y, y = tmin_dan, group = y), colour = "#619CFF")+
  xlab("Julian [2455000 = 2010 godina]") + 
  ylab("Temperatura [°C]") +
  ggtitle("Minimalne vrednosti temperature po godinama [°C]")+
  theme_bw()

j_temp
grid.arrange(tm_y, tmax_y, tmin_y, ncol = 1, nrow = 3)


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Histogram
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

df_t <- na.omit(stfdf_temp@data)

#TMEAN

hist_tmean <- ggplot(data = df_t, aes(tmean)) + 
  geom_histogram(binwidth = 1,
                 col="red", 
                 aes(y = ..density..,
                     fill = ..count..)) +
  scale_fill_gradient("Count", low = "blue", high = "red")+
  stat_function(fun = dnorm, 
                color = "black",
                size = 1.5,
                args = list(mean = mean(df_t$tmean), sd = sd(df_t$tmean)))+
  labs(title="Histogram prosecne vrednosti Temperature vazduha") +
  labs(x="Prosecna temperatura vazduha [°C]", y="Density") + 
  annotate("text", x = -10, y=0.05, label = paste("Mean:",round(mean(df_t$tmean),digits = 2), "°C"))+
  annotate("text", x = -10, y=0.045, label = paste("SD:",round(sd(df_t$tmean),digits = 2), "°C"))+
  scale_x_continuous(limits = c(-25, 50)) + 
  scale_y_continuous(limits = c(0, 0.065))+
  theme_bw()

#TMAX

hist_tmax <- ggplot(data = df_t, aes(tmax)) + 
  geom_histogram(binwidth = 1,
                 col="red", 
                 aes(y = ..density..,
                     fill = ..count..)) +
  scale_fill_gradient("Count", low = "blue", high = "red")+
  stat_function(fun = dnorm, 
                color = "black",
                size = 1.5,
                args = list(mean = mean(df_t$tmax), sd = sd(df_t$tmax)))+
  labs(title="Histogram maksimalne vrednosti Temperature vazduha") +
  labs(x="Maksimalna temperatura vazduha [°C]", y="Density") + 
  annotate("text", x = -10, y=0.05, label = paste("Mean:",round(mean(df_t$tmax),digits = 2), "°C"))+
  annotate("text", x = -10, y=0.045, label = paste("SD:",round(sd(df_t$tmax),digits = 2), "°C"))+
  scale_x_continuous(limits = c(-25, 50)) + 
  scale_y_continuous(limits = c(0, 0.065))+
  theme_bw()

#TMIN

hist_tmin <- ggplot(data = df_t, aes(tmin)) + 
  geom_histogram(binwidth = 1,
                 col="red", 
                 aes(y = ..density..,
                     fill = ..count..)) +
  scale_fill_gradient("Count", low = "blue", high = "red")+
  stat_function(fun = dnorm, 
                color = "black",
                size = 1.5,
                args = list(mean = mean(df_t$tmin), sd = sd(df_t$tmin)))+
  labs(title="Histogram minimalne vrednosti Temperature vazduha") +
  labs(x="Minimalna temperatura vazduha [°C]", y="Density") + 
  annotate("text", x = -10, y=0.05, label = paste("Mean:",round(mean(df_t$tmin),digits = 2), "°C"))+
  annotate("text", x = -10, y=0.045, label = paste("SD:",round(sd(df_t$tmin),digits = 2), "°C"))+
  scale_x_continuous(limits = c(-25, 50)) + 
  scale_y_continuous(limits = c(0, 0.065))+
  theme_bw()

grid.arrange(hist_tmean, hist_tmax, hist_tmin, ncol = 1, nrow = 3)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Hovmoller graficki prikaz
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ::::::::::::::::::::::::::::::::::::::::::::::::::::
# Hovmoller plot - latitude
# :::::::::::::::::::::::::::::::::::::s:::::::::::::::

lim_lat <- range(tempv$lat) # latitude range
lim_t <- range(tempv$julian) # time range
lat_axis <- seq(lim_lat[1], # latitude axis
                lim_lat[2],
                length=25)
t_axis <- seq(lim_t[1], # time axis
              lim_t[2],
              length=100)
lat_t_grid <- expand.grid(lat = lat_axis,
                          t = t_axis)
tempv_grid <- tempv
dists <- abs(outer(tempv$lat, lat_axis, "-"))
tempv_grid$lat <- lat_axis[apply(dists, 1, which.min)]

# ::::::::::::::::::::::::::::::::::::::::::::::::::::

tmean_lat_Hov <- group_by(tempv_grid, lat, julian) %>%
  summarise(tmean = mean(tmean))

tmean_Hovmoller_lat <- ggplot(tmean_lat_Hov) + # take data
  geom_tile(aes(x = lat, y = julian, fill = tmean)) + # plot
  fill_scale(name = "Prosecna vrednost Temperature [°C]") + # add color scale
  scale_y_reverse() + # rev y scale
  ylab("Dan") + # add y label
  xlab("Latitude (deg)") + # add x label
  ggtitle("Hovmoller plot - latitude")+
  theme_bw()+
  theme(legend.position="bottom")

tmax_lat_Hov <- group_by(tempv_grid, lat, julian) %>%
  summarise(tmax = mean(tmean))

tmax_Hovmoller_lat <- ggplot(tmax_lat_Hov) + # take data
  geom_tile(aes(x = lat, y = julian, fill = tmax)) + # plot
  fill_scale(name = "Maksimalna vrednost Temperature [°C]") + # add color scale
  scale_y_reverse() + # rev y scale
  ylab("Dan") + # add y label
  xlab("Latitude (deg)") + # add x label
  ggtitle("Hovmoller plot - latitude")+
  theme_bw()+
  theme(legend.position="bottom")

tmin_lat_Hov <- group_by(tempv_grid, lat, julian) %>%
  summarise(tmin = mean(tmin))

tmin_Hovmoller_lat <- ggplot(tmin_lat_Hov) + # take data
  geom_tile(aes(x = lat, y = julian, fill = tmin)) + # plot
  fill_scale(name = "Minimalna vrednost temperature [°C]") + # add color scale
  scale_y_reverse() + # rev y scale
  ylab("Dan") + # add y label
  xlab("Latitude (deg)") + # add x label
  ggtitle("Hovmoller plot - latitude")+
  theme_bw()+
  theme(legend.position="bottom")

grid.arrange(tmean_Hovmoller_lat, tmax_Hovmoller_lat, tmin_Hovmoller_lat, ncol = 2, nrow = 2)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::
# Hovmoller plot - longitude
# ::::::::::::::::::::::::::::::::::::::::::::::::::::

lim_lon <- range(tempv$lon) # latitude range
lim_t <- range(tempv$julian) # time range
lon_axis <- seq(lim_lon[1], # latitude axis
                lim_lon[2],
                length=25)

t_axis <- seq(lim_t[1], # time axis
              lim_t[2],
              length=100)
lon_t_grid <- expand.grid(lon = lon_axis,
                          t = t_axis)
tempvv_grid <- tempv
dists <- abs(outer(tempv$lon, lon_axis, "-"))
tempvv_grid$lon <- lon_axis[apply(dists, 1, which.min)]

# ::::::::::::::::::::::::::::::::::::::::::::::::::::

tmean_lon_Hov <- group_by(tempvv_grid, lon, julian) %>%
  summarise(tmean = mean(tmean))
tmean_Hovmoller_lon <- ggplot(tmean_lon_Hov) + # take data
  geom_tile(aes(x = lon, y = julian, fill = tmean)) + # plot
  fill_scale(name = "Prosecna vrednost Temperature [°C]") + # add color scale
  scale_y_reverse() + # rev y scale
  ylab("Dan") + # add y label
  xlab("Longitude (deg)") + # add x label
  ggtitle("Hovmoller plot - longitude")+
  theme_bw()+
  theme(legend.position="bottom")

tmax_lon_Hov <- group_by(tempvv_grid, lon, julian) %>%
  summarise(tmax = mean(tmax))
tmax_Hovmoller_lon <- ggplot(tmax_lon_Hov) + # take data
  geom_tile(aes(x = lon, y = julian, fill = tmax)) + # plot
  fill_scale(name = "Maksimalna vrednost Temperature [°C]") + # add color scale
  scale_y_reverse() + # rev y scale
  ylab("Dan") + # add y label
  xlab("Longitude (deg)") + # add x label
  ggtitle("Hovmoller plot - longitude")+
  theme_bw()+
  theme(legend.position="bottom")

tmin_lon_Hov <- group_by(tempvv_grid, lon, julian) %>%
  summarise(tmin = mean(tmin))
tmin_Hovmoller_lon <- ggplot(tmin_lon_Hov) + # take data
  geom_tile(aes(x = lon, y = julian, fill = tmin)) + # plot
  fill_scale(name = "Minimalna vrednost temperature [°C]") + # add color scale
  scale_y_reverse() + # rev y scale
  ylab("Dan") + # add y label
  xlab("Longitude (deg)") + # add x label
  ggtitle("Hovmoller plot - longitude")+
  theme_bw()+
  theme(legend.position="bottom")

grid.arrange(tmean_Hovmoller_lon, tmax_Hovmoller_lon, tmin_Hovmoller_lon, nrow = 2, ncol = 2)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::
# Data preparation and transformation
# ::::::::::::::::::::::::::::::::::::::::::::::::::::

# Preklop sa prediktorima
r <- raster("D:/Doktorske_akademske_studije_Geodezija_i_geoinformatika/Semestar_2/Prostorno_vremenska_statisitika/Temperature [Sale]/dem/dem.sdat") # DEM
e <- extract(r,stfdf_temp@sp) 
stfdf_temp@sp$dem=e

r <- raster("D:/Doktorske_akademske_studije_Geodezija_i_geoinformatika/Semestar_2/Prostorno_vremenska_statisitika/Temperature [Sale]/dem/twi.sdat") # Croatian TWI
e <- extract(r,stfdf_temp@sp)
stfdf_temp@sp$twi=e

# U slucaju da je potrebno ukoliniti stanice bez merenja 

# ggplot(data = tempv, aes(x = julian, y = tmean))+
#   geom_line(colour = "red", alpha = 0.5)+
#   facet_wrap(~ staid, ncol = 8)+
#   ggtitle("Vremenska linija Prosecnih temperatura po stanicama")+
#   theme_bw()
# 
# stsdf_utm@sp <- stsdf_utm@sp[!(stsdf_utm@sp$staid == 13619), ] # brisanje stanice bez merenja


# ::::::::::::::::::::::::::::::::::::::::::::::::::::
# Geometrical temperature trend and residuals
# ::::::::::::::::::::::::::::::::::::::::::::::::::::
#year=2009:2018
#for(i in year){
#  time=seq(as.Date(paste(i,"-01-01", sep="")), as.Date(paste(i,"-12-31", sep="")), by="day")
#}

time <- stfdf_temp@time
days=gsub("-","",time,fixed=TRUE)
daysNum = length(time)

var_mean = "mean"
var_max = "max"
var_min = "min"

# Koeficijenti
data(tregcoef)
coef_tmean = as.vector(tregcoef$tmeanGSODECAD_noMODIS) # coefficients for tmean
coef_tmax = as.vector(tregcoef$tmaxGSODECAD_noMODIS) # coefficients for tmean
coef_tmin = as.vector(tregcoef$tminGSODECAD_noMODIS) # coefficients for tmean

# Geometrijski temperaturni trend
gtt_mean <- tgeom2STFDF(stfdf_temp@sp, time = time, variable = var_mean)
gtt_max <- tgeom2STFDF(stfdf_temp@sp, time = time, variable = var_max)
gtt_min <- tgeom2STFDF(stfdf_temp@sp, time = time, variable = var_min)

# ::::: # Kontrola (jednak prostorni domen)
identical(gtt_mean@sp, stfdf_temp@sp)

# Racunanje reziduala za TMEAN
stfdf_temp@data$res_tmean = stfdf_temp@data$tmean - (coef_tmean[1] + coef_tmean[2]*as.numeric(gtt_mean@data$temp_geo)+coef_tmean[3]*rep(as.numeric(stfdf_temp@sp$dem),daysNum) + coef_tmean[4]*rep(as.numeric(stfdf_temp@sp$twi),daysNum))

# Racunanje reziduala za TMAX
stfdf_temp@data$res_tmax = stfdf_temp@data$tmax - (coef_tmax[1] + coef_tmax[2]*as.numeric(gtt_max@data$temp_geo)+coef_tmax[3]*rep(as.numeric(stfdf_temp@sp$dem),daysNum) + coef_tmax[4]*rep(as.numeric(stfdf_temp@sp$twi),daysNum))

# Racunanje reziduala za TMIN
stfdf_temp@data$res_tmin = stfdf_temp@data$tmin - (coef_tmin[1] + coef_tmin[2]*as.numeric(gtt_min@data$temp_geo)+coef_tmin[3]*rep(as.numeric(stfdf_temp@sp$dem),daysNum) + coef_tmin[4]*rep(as.numeric(stfdf_temp@sp$twi),daysNum))

# Transformacija u STSDF klasu
stsdf_temp <- as(stfdf_temp, "STSDF")
# Transformacija podataka u UTM projekciju
stsdf_utm <- spTransform(stfdf_temp, CRS("+init=epsg:32634"))

# ::::: # Kontrola
length(stsdf_utm@data$tmean)
length(gtt_mean@data$temp_geo)
length(rep(as.numeric(stsdf_utm@sp$dem),daysNum))

stsdf_utm@data$res_tmean

summary(stfdf_temp)
summary(stsdf_utm)

summary(stfdf_temp@data$res_tmean)
summary(stsdf_utm@data$res_tmean)






































