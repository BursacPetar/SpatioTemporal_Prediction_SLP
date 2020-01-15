# Prostorno - vremenska statistika
# Seminarski zadatak

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

wdir <- 'd:/Doktorske_akademske_studije_Geodezija_i_geoinformatika/Semestar_2/Prostorno_vremenska_statisitika/Seminarski_zadatak/'
setwd(wdir)
getwd()

# Ucitavanje podataka
load(file = "Ulazni_podaci/ogimet_st.rda") # stations - stanice
load(file = "Ulazni_podaci/ogimet_slp_stfdf.rda") # stfdf - STFDF sa pritiskom (slp) za 10 godina (2009-2018)
load(file = "Ulazni_podaci/ogimet_obs.rda") # observations - merenja

# ===================================================================================
# Data analysis, preparation and wrangling
# ===================================================================================

crs <- "+init=epsg:4326"
coords <- c("lon", "lat")
stations_sf <- st_as_sf(stations, crs = crs, coords = coords)

#st_write(stations_sf, "Ulazni_podaci/Stanice.shp")

# Prostorna dispozicija mernih stanica vazdusnog pritiska
dat=map_data("world")
srb_plot <- ggplot(stations_sf)+
  geom_sf(aes(colour=altitude, size=altitude))+
  guides(size=FALSE)+
  labs(x = "Longitude [deg]",
       y = "Latitude [deg]",
       title = "Prostorna dispozicija mernih stanica",
       caption = "Coordinate Reference System - WGS84")+
  #geom_sf_text(aes(label = staid))+
  #geom_sf_text_repel(aes(label = staid), nudge_x = -0.2)+
  geom_map(data=dat[dat$region=="Serbia",], map=dat[dat$region=="Serbia",],
           aes(x=long, y=lat, map_id=region),
           color="white", fill="#7f7f7f", size=0.05, alpha=1/4) +
  theme_bw()+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) # temu bolje napisi (kao kartu Debeljace)

# Vizuelizacija stanica
srb_plot_web <- mapview(stations_sf, z = "altitude", layer.name = "Visina mernih stanica [m]")
rm(srb_plot_web)

# Broj stanica po drzavama
stations %>% group_by(country) %>%
 summarise(n = n())

dim(stfdf)
# stfdf - full grid - kombinacija bilo kog sp objekta i bilo kog xts objekta kako bi se repreyentovale sve moguce lokacije
# u dostupnom prostornom i vremenskom domenu.
class(stfdf)
typeof(stfdf)

slotNames(stfdf)
summary(stfdf@data)[7] # ukupno 23527 vrmenskih trenutaka bez merenja, od
count(stfdf@data)

str(stfdf)

#stfdf %>% group_by(timeIndex) %>%
#  summarise(mean_slp = mean(slp))
stfdf@sp@proj4string

# =============
# Visualization
# =============
observations # kolona slp [jedinica: hPa] (kilo-Paskala)

# Hovmoller diagram
stplot(stfdf[], mode="xt", scaleX=0, col.regions=terrain.colors(100))

# Time-series plots
stplot(stfdf[], mode="ts", auto.key=FALSE) # kumulativno na svim stanicama na jednom plot-u
# sa auto.key se uklanja legenda

stplot(stfdf[1:9,"2009::2018",'slp'], mode = "ts") # primer selekcije podataka u prostornom i vremenskom smislu

stplot(stfdf[], mode = "tp") # kumulativno po stanicama, svaka poseban plot
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Plot za dati vremenski interval
world <- getMap(resolution="low")
stplot(stfdf[,"2011-08-23::2011-08-28"],col.regions= brewer.pal(6,"Spectral"),cuts=6,sp.layout=list(world,first=TRUE))
stplot(stfdf[,"2011-08-23::2011-08-28"], mode = "ts", auto.key=FALSE)

# ======================================
# Exploring spatio-temporal dependencies
# ======================================
var <- variogramST(slp ~ 1, stfdf[])
a1<-plot(var) # sample variograms at each time lag
b1<-plot(var, map=FALSE) # ST-variogram map
c1<-plot(var, wireframe=TRUE) # ST-variogram wireframe
d1 <- plot(var, diff = TRUE)
grid.arrange(a1,b1,c1, nrow = 1, ncol = 3)

# =============================================================================
head(observations)
head(stations)

df_slp <- as.data.frame(stfdf)
df_slp <- df_slp[!is.na(df_slp$slp), ]
df_slp$time <- str_sub(df_slp$time, end=-1)

obs <- df_slp
sta <- stations

obs$date_1 <- obs$time
obs %<>% separate(date_1, "-", into = c("y", "m", "d"))

# Join sa prostornim podacima - lokacijama stanica

#slpr <- left_join(obs, sta, by = "staid")
slpr <- obs

slpr %<>% mutate(y = as.numeric(y),
                 m = as.numeric(m),
                 d = as.numeric(d))

# =============================================================================
# Converting to Julian date
# ==========================
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

toJulian(2009, 12, 1)

slpr$julian <- toJulian(slpr$y, slpr$m, slpr$d)
# ==============================================================================

# Sumarni prikaz po godinama
summ <- group_by(slpr, y) %>%
  summarise(mean_slp = mean(slp))

summ$mean_slp

# Vrednost normalnog vazdusnog pritiska je 1013mb na morskom nivou
# Primer: sve stanice u junu mesecu po godinama gde je pritisak bio iznad normale
obs_6 <- filter(slpr , m == 6)
summ_6 <- group_by(obs_6, y, staid) %>%
  summarise(dani_iznad = sum(slp > 1013))

# Medijana ukupnog broja dana za sve godine
median(summ_6$dani_iznad)

ggplot(summ_6, aes(x = staid, y = dani_iznad, group = y, colour = y))+
  geom_line() +
  labs(x = "Stanica", y = "Broj dana", title = "Pritisak iznad normale u junu mesecu od 2009-2018")+
  theme(axis.text.x=element_blank())

#The first task faced by the spatio-temporal modeler is data visualization. This is an
#important preliminary task that needs to be carried out prior to the exploratory-data-analysis
#stage and the modeling stages.

#slpr_1 <- subset(slpr, d %in% c('01', '15', '30'))

# Time-series plots
UIDs <- unique(slpr$staid) # extract IDs
UIDs_sub <- sample(UIDs, 10) # sample 10 IDs
slpr_sub <- filter(slpr, staid %in% UIDs_sub) # subset data

slpr_TS <- ggplot(slpr_sub) +
  geom_line(aes(x = julian, y = slp, color = 'blue', alpha = 0.3)) + # line plot of z against t
  facet_wrap(~staid, ncol = 5) + # facet by station
  xlab("Day number (days)") + # x label
  ylab("SLP (mb)") + # y label
  theme_bw() + # BW theme
  theme(panel.spacing = unit(1, "lines"),legend.position = "none") # facet spacing

slpr_TS
# =============================================================================
# Hovmoller plots

lim_lat <- range(slpr$lat) # latitude range
lim_t <- range(slpr$julian) # time range
lat_axis <- seq(lim_lat[1], # latitude axis
                lim_lat[2],
                length=25)
t_axis <- seq(lim_t[1], # time axis
              lim_t[2],
              length=100)
lat_t_grid <- expand.grid(lat = lat_axis,
                          t = t_axis)


slpr_grid <- slpr
dists <- abs(outer(slpr$lat, lat_axis, "-"))
slpr_grid$lat <- lat_axis[apply(dists, 1, which.min)]

slpr_lat_Hov <- group_by(slpr_grid, lat, julian) %>%
  summarise(slp = mean(slp))
#Every latitude and time band contains at least one data point, so that the Hovmoller
#plot contains no missing points on the established grid.

Hovmoller_lat <- ggplot(slpr_lat_Hov) + # take data
  geom_tile(aes(x = lat, y = julian, fill = slp)) + # plot
  fill_scale(name = "Normalni vazdusni pritisak [mb]") + # add color scale
  scale_y_reverse() + # rev y scale
  ylab("Dan") + # add y label
  xlab("Latitude (deg)") + # add x label
  ggtitle("Hovmoller plot - latitude")+
  theme_bw()


lim_lon <- range(slpr$lon) # latitude range
lim_t <- range(slpr$julian) # time range
lon_axis <- seq(lim_lon[1], # latitude axis
                lim_lon[2],
                length=25)
t_axis <- seq(lim_t[1], # time axis
              lim_t[2],
              length=100)
lon_t_grid <- expand.grid(lon = lon_axis,
                          t = t_axis)

slpr_grid <- slpr
dists <- abs(outer(slpr$lon, lon_axis, "-"))
slpr_grid$lon <- lon_axis[apply(dists, 1, which.min)]

slpr_lon_Hov <- group_by(slpr_grid, lon, julian) %>%
  summarise(slp = mean(slp))

Hovmoller_lon <- ggplot(slpr_lon_Hov) + # take data
  geom_tile(aes(x = lon, y = julian, fill = slp)) + # plot
  fill_scale(name = "Normalni vazdusni pritisak [mb]") + # add color scale
  scale_y_reverse() + # rev y scale
  ylab("Dan") + # add y label
  xlab("Longitude (deg)") + # add x label
  ggtitle("Hovmoller plot - longitude")+
  theme_bw()

grid.arrange(Hovmoller_lon, Hovmoller_lat, nrow = 1, ncol = 2)
# =============================================================================
# Exploratory Data Analysis
# =============================================================================
set.seed(1)

# Nova varijabla koja broji dane od 1-3655 (ukupan broj dana)
head(slpr)
slpr$t <- slpr$julian - 2454832
range(slpr$t)

# =============================================================================
# Empirijska prostorna sredina (Empirical Spatial Mean)

# The empirical spatial mean is a
# spatial quantity that can be stored in a new data frame that contains the spatial locations
# and the respective average pressure value at each location.
spat_av <- group_by(slpr, lat, lon) %>% # group by lon-lat
  summarise(spm_emp = mean(slp)) # mean for each lon-lat

lat_means <- ggplot(spat_av) +
  geom_point(aes(lat, spm_emp, colour = spm_emp)) +
  xlab("Latitude (deg)") +
  ylab("Average sea level pressure (mb)") +
  ggtitle("Empirical Spatial Mean - latitude")+
  theme_bw()
lon_means <- ggplot(spat_av) +
  geom_point(aes(lon, spm_emp,colour = spm_emp)) +
  xlab("Longitude (deg)") +
  ylab("Average sea level pressure (mb)") +
  ggtitle("Empirical Spatial Mean - longitude")+
  theme_bw()

grid.arrange(lon_means, lat_means, nrow = 1, ncol = 2)

# Scatter plot with correlation coefficient
#:::::::::::::::::::::::::::::::::::::::::::::::::
gg1 <- ggscatter(spat_av, x = "lon", y = "spm_emp",
                 color = "red", alpha = 0.5,
                 add = "loess",  # Add regressin line ("reg.line") or local regression fitting ("loess")
                 add.params = list(color = "blue", fill = "orange"), # Customize reg. line
                 conf.int = TRUE)+ # Add confidence interval
  stat_cor(method = "pearson")+# Add correlation coefficient
  xlab("Longitude (deg)") +
  ylab("Average sea level pressure (mb)") +
  ggtitle("Empirical Spatial Mean - longitude")

#:::::::::::::::::::::::::::::::::::::::::::::::::
gg2 <- ggscatter(spat_av, x = "lat", y = "spm_emp",
                 color = "red", alpha = 0.5,
                 add = "loess",
                 add.params = list(color = "blue", fill = "orange"),
                 conf.int = TRUE)+
  stat_cor(method = "pearson")+
  xlab("Latitude (deg)") +
  ylab("Average sea level pressure (mb)") +
  ggtitle("Empirical Spatial Mean - latitude")

grid.arrange(gg1, gg2, ncol = 2, nrow = 1)

# =============================================================================
# Emprijska vremenska sredina (Empirical Temporal Means)

slpr_dan <- slpr %>% 
  dplyr::group_by(y, m, d) %>%
  summarise(slp_dan = mean(slp), 
            julian = mean(julian))
# Po danu
slp_av_d <- group_by(slpr,d) %>%
  summarise(meanSlp = mean(slp))

em_mean_d <- ggplot() +
  geom_line(data = slpr, aes(x = d, y = slp, group = staid),
            colour = "red", alpha = 0.04) +
  geom_line(data = slp_av_d, aes(x = d, y = meanSlp)) +
  geom_boxplot(data = slpr_dan, aes(x=d, y=slp_dan, group = d, colour =d), show.legend = FALSE)+
  xlab("Dan") + ylab("Normalni vazdusni pritisak [mb]") +
  ggtitle("Emprijska vremenska sredina - po danu")+
  theme_bw()

# Po mesecu
slp_av_m <- group_by(slpr,m) %>%
  summarise(meanSlp = mean(slp))

em_mean_m <- ggplot() +
  geom_line(data = slpr, aes(x = m, y = slp, group = staid),
            colour = "red", alpha = 0.04) +
  geom_line(data = slp_av_m, aes(x = m, y = meanSlp)) +
  geom_boxplot(data = slpr_dan, aes(x=m, y=slp_dan, group = m, colour =m), show.legend = FALSE)+
  xlab("Mesec") + ylab("Normalni vazdusni pritisak [mb]") +
  ggtitle("Emprijska vremenska sredina - po mesecu")+
  scale_x_discrete(limits=c(0:12))+
  theme_bw()

# Po godini
slp_av_y <- group_by(slpr,y) %>%
  summarise(meanSlp = mean(slp))

em_mean_y <- ggplot() +
  geom_line(data = slpr, aes(x = y, y = slp, group = staid),
            colour = "red", alpha = 0.04) +
  geom_line(data = slp_av_y, aes(x = y, y = meanSlp)) +
  geom_boxplot(data = slpr_dan, aes(x=y, y=slp_dan, group = y, colour =y), show.legend = FALSE) +
  xlab("Godina") + ylab("Normalni vazdusni pritisak [mb]") +
  ggtitle("Emprijska vremenska sredina - po godini") +
  theme_bw()

grid.arrange(em_mean_d, em_mean_m, em_mean_y, ncol = 1, nrow = 3)

# Plot po uzoru u knjzi
slp_av <- group_by(slpr, julian) %>%
  summarise(meanSlp = mean(slp))
gSLPav <- ggplot() +
  geom_line(data = slpr, aes(x = julian, y = slp, group = staid),
            colour = "red", alpha = 0.04) +
  geom_line(data = slp_av, aes(x = julian, y = meanSlp)) +
  xlab("Julian date (2454000 = 2009)") + ylab("Average sea level pressure (mb)") +
  ggtitle("Empirical Temporal Means")+
  theme_bw()
gSLPav

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ggplot(data = slpr, aes(x = julian, y = slp))+
  geom_line(colour = "red", alpha = 0.5)+
  facet_wrap(~ staid, ncol = 8)+
  ggtitle("Vremenska linija merenja Normalnog vazdusnog pritiska po stanicama")
  theme_bw()
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  #gSLPav_1 <- ggplot() +
#  geom_line(data = slp_av, aes(x = julian, y = meanSlp),
#            colour = "DeepSkyBlue", alpha = 0.5) +
#  xlab("Year") + ylab("Average sea level pressure (mb)") +
#  theme_bw()

#grid.arrange(gSLPav, gSLPav_1, nrow = 2, ncol = 1)

#na_count <-sapply(slpr, function(y) sum(length(which(is.na(y)))))
#na_count <- data.frame(na_count)
#na_count_obs <-sapply(stfdf@data, function(y) sum(length(which(is.na(y))))) %>% as.data.frame()
#length(slpr$slp)
#length(stfdf@data$slp) - length(obs_1$slp)
#length(observations$slp)
#23527
#obs_1 <- observations
#obs_1 <- obs_1[obs_1$slp > 870 , ]
#length(obs_1$slp)
#summary(observations$slp)
#summary(stfdf@data$slp)

# Empirijska kovarijansa (Empirical Covariances)

# Linearni model (Linear model) - regresioni model
lm <- lm(slp ~ lat + lon + t, data = slpr)
summary(lm)

par(mfrow = c(2, 2))
plot(lm)
par(mfrow = c(1, 1))

anova(lm) 	         # ANOVA tabela
coefficients(lm)    # Koeficijenti modela
confint(lm) 	       # Intervali poverenja za koef.
deviance(lm) 	     # Suma kvadrata reziduala
effects(lm)         # Vektor ortogonalnog efekta
fitted(lm) 	       # Vektor prediktovane Y vrednosti
residuals(lm)       # Reziduali
summary(lm)	       # R2; F statistika; standardna greška reziduala
vcov(lm) 	         # Varijans-kovarijans matrica parametara

# Scatter plot with correlation coefficient
#:::::::::::::::::::::::::::::::::::::::::::::::::
sp1 <- ggscatter(slpr, x = "lon", y = "slp",
                color = "red", alpha = 0.05,
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE)+ # Add confidence interval
  stat_cor(method = "pearson")# Add correlation coefficient

#:::::::::::::::::::::::::::::::::::::::::::::::::
sp2 <- ggscatter(slpr, x = "lat", y = "slp",
                color = "red", alpha = 0.05,
                add = "reg.line",
                add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE)+
  stat_cor(method = "pearson")

#:::::::::::::::::::::::::::::::::::::::::::::::::
sp3 <- ggscatter(slpr, x = "julian", y = "slp",
                color = "red", alpha = 0.05,
                add = "reg.line",
                add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE)+
  stat_cor(method = "pearson")

#:::::::::::::::::::::::::::::::::::::::::::::::::
slpr$julian2 <- slpr$julian^2
sp4 <- ggscatter(slpr, x = ("julian2"), y = "slp",
                color = "red", alpha = 0.05,
                add = "reg.line",
                add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE)+
  stat_cor(method = "pearson")

#:::::::::::::::::::::::::::::::::::::::::::::::::
grid.arrange(sp1, sp2, sp3, sp4, nrow = 2, ncol = 2)

# Linearni model (Linear model) - regresioni model
lm <- lm(slp ~ lat + lon + julian + julian2, data = slpr)
summary(lm)

# Values per station
sta_mean <- group_by(slpr, staid) %>%
  summarise(slp_mean = mean(slp))

sta_max <- group_by(slpr, staid) %>%
  summarise(slp_max = max(slp))

sta_min <- group_by(slpr, staid) %>%
  summarise(slp_min = min(slp))

slpr$st_mean <- sta_mean$slp_mean[match(slpr$staid, sta_mean$staid) ]
slpr$st_max <- sta_max$slp_max[match(slpr$staid, sta_max$staid) ]
slpr$st_min <- sta_min$slp_min[match(slpr$staid, sta_min$staid) ]

lm <- lm(slp ~ lat + lon + julian + julian2 + st_mean + st_max + st_min, data = slpr)
summary(lm)

# ====================================================================================
# Preparing covariates
# ====================================================================================
library(raster)
library(stars)

dem <- raster::stack("Ulazni_podaci/SRTM_dem_90m/SRTM_DEM_90m_clipped.tif")

dem_stars <- st_as_stars(dem)
ggplot()+
  geom_stars(data = dem_stars)

slpr_sp <- slpr
coordinates(slpr_sp) <- ~lon+lat
proj4string(slpr_sp) <- CRS("+init=epsg:4326")

slpr_sp@data <- cbind(slpr_sp@data, elevation = raster::extract(dem, slpr_sp, df = TRUE))

slpr_sp@data %<>% rename(elevation = elevation.SRTM_DEM_90m_clipped) 

# ::::::::::::::::
# Plot - statioins
slpr_sf <- st_as_sf(slpr_sp)
slpr_sf_stat <- group_by(slpr_sf, staid) %>%
  summarise(elevation = max(elevation))
ggplot()+
  geom_sf(data = slpr_sf_stat,aes(colour = elevation))
# ::::::::::::::::

lm <- lm(slp ~ lat + lon + julian + julian2 + elevation, data = slpr_sp)
summary(lm)

# :::::::::::::::::::::::::::::::::::::::::::::
# Za plot
dem_coarse_res <- raster::stack("Ulazni_podaci/SRTM_dem_90m/SRTM_DEM.tif")
dem_c_stars <- st_as_stars(dem_coarse_res)

ggplot()+
  geom_stars(data = dem_c_stars)+
  geom_sf(data = slpr_sf, aes(colour = "cyan"))+
  ggtitle("SRTM_DEM as covariate")+
  xlab("Longitutde [deg]")+
  ylab("Latitude [deg]")+
  guides(colour = FALSE)+
  theme_bw()

# :::::::::::::::::::::::::::::::::::::::::::::
# Empirijska kovarijansa (Empirical Covariances)
# :::::::::::::::::::::::::::::::::::::::::::::

# Linearni model (Linear model) - regresioni model
slpr$elevation <- slpr_sp@data$elevation
lm <- lm(slp ~ lat + lon + julian + julian2 + elevation, data = slpr_long)
summary(lm)

slpr_long$residuals <- lm$residuals
names(lm)

# REMOVE any records that contain missing observations
slpr_long <- slpr
slpr_long <- slpr_long %>%
  filter(!is.na(slp))
dim(slpr_long)
dim(slpr)

na_count <-sapply(stfdf@data, function(y) sum(length(which(is.na(y))))) %>% as.data.frame()

stfdf

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Variogram
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
var <- variogramST(slp ~ 1, stfdf[])
a1<-plot(var) # sample variograms at each time lag
b1<-plot(var, map=FALSE) # ST-variogram map
c1<-plot(var, wireframe=TRUE) # ST-variogram wireframe
d1 <- plot(var, diff = TRUE)
grid.arrange(a1,b1,c1, nrow = 1, ncol = 3)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Variogram Modelling

# As in a normal 2D kriging experiment, at this point we need to fit a model to our variogram. 
# For doing so we will use the function vgmST and fit.StVariogram, 
# which are the spatio-temporal matches for vgm and fit.variogram.

# Four basic space-time model types supported by gstat
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
?vgmST
?fit.StVariogram

# Primer
# demo(stkrige) 

# Estimation of the spatio-temporal anisotropy without an underlying spatio-temporal model.
estiStAni(var, c(0, 10), "metric",
          vgm(psill=4,"Sph", range=3, nugget=0),
          vgm(psill=8,"Sph", range=3, nugget=0) )

# Lower (and upper) bounds
pars.l <- c(sill.s = 0, range.s = 1, nugget.s = 0,
            sill.t = 0, range.t = 1, nugget.t = 0,
            sill.st = 0, range.st = 1, nugget.st = 0,
            anis = 0)

sumMetric <- vgmST("sumMetric",
                   space = vgm(psill=3,"Sph", range=3, nugget=0.1),
                   time = vgm(psill=7,"Sph", range=3, nugget=0.1),
                   joint = vgm(psill=1,"Sph", range=2, nugget=0.1),
                   stAni=5)

sumMetric_vgm <- fit.StVariogram(var, sumMetric, method="L-BFGS-B",lower=pars.l)

attr(sumMetric_vgm, "MSE")
sumMetric_vgm

sumMetric_vgm$space$range 
sumMetric_vgm$joint$range 
sumMetric_vgm$stAni 

pp1 <- plot(var, wireframe = T,
            zlab=NULL,
            xlab=list("distance (km)", rot=30, cex=0.7),
            ylab=list("time lag (days)", rot=-35, cex=0.7),
            scales=list(arrows=F, z = list(distance = 5), cex=0.5), 
            main="3d sample variogram")
pp2 <- plot(var, list(sumMetric_vgm), wireframe = T, 
            xlab=list("distance (km)", rot=30, cex=0.7),
            ylab=list("time lag (days)", rot=-35, cex=0.7),
            scales=list(arrows=F, z = list(distance = 5), cex=0.5), main="3d fitted sum-metric variogram")
grid.arrange(pp1, pp2, ncol = 2, nrow = 1)

pp3 <- plot(var, map = F)
pp4 <- plot(var, list(sumMetric_vgm), map = F)
grid.arrange(pp3, pp4, ncol = 2, nrow = 1)

grid.arrange(pp1, pp2, pp3, pp4, ncol = 2, nrow = 2)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Grid - new data
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

grid <- as(dem_coarse_res, "SpatialPointsDataFrame")
# plot(grid)
proj4string(grid)
proj4string(grid) <- proj4string(stfdf@sp) 

#st_make_grid()
#tm.grid <- seq(as.POSIXct('2009-01-01 UTC'),as.POSIXct('2018-12-31 UTC'),length.out=3650)
#summary(tm.grid)
#summary(stfdf@time)
#range(stfdf@time)
grid@data <- grid@data[-1]
grid_for_ok <- STF(grid, 
                   stfdf@time, 
                   stfdf@endTime)
proj4string(stfdf@sp)
proj4string(grid_for_ok@sp) <- proj4string(stfdf@sp) 

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Ordinary global Spatio-Temporal Kriging
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ok_stfdf <- krigeST(as.formula(slp ~ 1),
                    data = stfdf,
                    newdata = grid_for_ok,
                    modelList = sumMetric_vgm,
                    computeVar = TRUE)

na_count_obs <-sapply(stfdf@data, function(y) sum(length(which(is.na(y))))) %>% as.data.frame()

# :::::::::::::::::::::::::::::::::::::::::
# Remove stations without observations
# count NAs per stations
temp = stfdf
summary(temp)
nrowsp <- length(temp@sp)
# 65 ukupan broj stanica
numNA <- apply(matrix(temp@data[,'slp'],
                      nrow=nrowsp,
                      byrow=F), 
               MARGIN=1,
               FUN=function(x) sum(is.na(x)))
sum(numNA)
# 23527
# Remove stations 
rem <- numNA != length(time)
temp <-  temp[rem, drop=F]
summary(temp)
nrowsp_temp <- length(temp@sp)
# 62 - ukupan broj stanica
na_count_temp <-sapply(temp@data, function(y) sum(length(which(is.na(y))))) %>% as.data.frame()
# 23524

# :::::::::::::::::::::::::::::::::::::::::
# Stations removed - remove NA observations
#nrowobs <- length(temp@data$slp)
## 65 ukupan broj stanica
#numNA_slp <- apply(matrix(temp@data[,'slp'],
#                      nrow=nrowobs,
#                      byrow=T), 
#               MARGIN=1,
#               FUN=function(x) is.na(x))
#sum(numNA_slp)
#rem_obs <- numNA_slp != length(time)
#temp <-  temp[rem_obs, drop=F]
temp@data %<>% 
  drop_na(slp)

temp = stfdf

#temp_stsdf <- as(temp, "STSDF") 
#summary(temp_stsdf)

# Replace NA attribute values; disaggregation time series
temp <- na.locf(temp, na.rm = TRUE)
summary(temp)


stplot(temp[], mode="xt", scaleX=0, col.regions=terrain.colors(100))
stplot(temp[], mode="ts", auto.key=FALSE) # kumulativno na svim stanicama na jednom plot-u
stplot(temp[1:9,"2009::2018",'slp'], mode = "ts") # primer selekcije podataka u prostornom i vremenskom smislu
stplot(temp[], mode = "tp") # kumulativno po stanicama, svaka poseban plot

temp_1 <- temp
ok_stfdf <- krigeST(as.formula(slp ~ 1),
                    data = temp,
                    newdata = STF(grid, 
                                  temp_1@time,
                                  temp_1@endTime),
                    modelList = sumMetric_vgm,
                    computeVar = TRUE)

# :::::::::::::::
# NOVI GRID - 1KM
# :::::::::::::::

dem_1km <- raster::stack("Ulazni_podaci/SRTM_dem_90m/SRTM_DEM_90m_clipped_4326_1km.tif")
dem_1km_stars <- st_as_stars(dem_1km)

ggplot()+
  geom_stars(data = dem_1km_stars)+
  ggtitle("SRTM_DEM as covariate")+
  xlab("Longitutde [deg]")+
  ylab("Latitude [deg]")+
  theme_bw()

# Make grid from dem 1km, as new locations for prediction
grid_1km <- as(dem_1km, "SpatialPointsDataFrame")
# plot(grid)
proj4string(grid_1km)
proj4string(grid_1km) <- proj4string(stfdf@sp) 

grid_1km@data <- grid_1km@data[-1]
grid_1km_ok <- STF(grid_1km, 
                   temp@time, 
                   temp@endTime)

ok_stfdf <- krigeST(as.formula(slp ~ 1),
                    data = temp,
                    newdata = STF(grid, 
                                  temp_1@time,
                                  temp_1@endTime),
                    modelList = sumMetric_vgm,
                    computeVar = F)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

show.vgms()
vgm()

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Subset data
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
slp_df <- as.data.frame(stfdf)
slpr_dan <- slpr %>% 
  dplyr::group_by(y, m, d) %>%
  summarise(slp_dan = mean(slp), 
            julian = mean(julian))

slpr_dan$vreme <- paste(slpr_dan$y, "-", slpr_dan$m, "-", slpr_dan$d)

ggp1 <- ggplot()+
  geom_line(data = slpr_dan, aes(x = julian, y = slp_dan, group = y, colour = y))+
  scale_colour_gradient(low="blue", high="red")+
  xlab("Julian date [2455000 = 2010 year]") + ylab("Sea level pressure [mb]") +
  ggtitle("Timeline of sea level pressure")+
  theme_bw()

ggp3 <- ggplot(slpr_dan, aes(x=y, y=slp_dan, group = y, colour =y)) +
  geom_boxplot()+
  scale_colour_gradient(low="blue", high="red")+
  xlab("Year") + ylab("Sea level pressure [mb]") +
  ggtitle("Sea level pressure per year")+
  theme_bw()
  
grid.arrange(ggp1, ggp3, nrow = 2, ncol = 1)

stfdf_sub <- stfdf[, "2017-01-01::2017-12-31"]
dim(stfdf_sub)
summary(stfdf_sub)
#grid_1km_ok <- STF(grid_1km, 
#                   stfdf_sub@time, 
#                   stfdf_sub@endTime)

#summary(grid_1km_ok)

temp <- na.locf(stfdf, na.rm = TRUE)
summary(temp)
temp_sub <- temp[ ,"2017-01-01::2017-12-31"]
summary(temp_sub)
grid_1km_ok_sub <- STF(grid_1km, 
                   temp_sub@time, 
                   temp_sub@endTime)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Variogram Modelling
var_sub <- variogramST(slp ~ 1, temp_sub)
a2<-plot(var_sub) # sample variograms at each time lag
b2<-plot(var_sub, map=FALSE) # ST-variogram map
c2<-plot(var_sub, wireframe=TRUE) # ST-variogram wireframe
d2 <- plot(var_sub, diff = TRUE)
grid.arrange(a2,b2,c2, nrow = 1, ncol = 3)
# Estimation of the spatio-temporal anisotropy without an underlying spatio-temporal model.
estiStAni(var_sub, c(0, 10), "metric",
          vgm(psill=4,"Sph", range=3, nugget=0),
          vgm(psill=8,"Sph", range=3, nugget=0) )
# Lower (and upper) bounds
pars.l <- c(sill.s = 0, range.s = 1, nugget.s = 0,
            sill.t = 0, range.t = 1, nugget.t = 0,
            sill.st = 0, range.st = 1, nugget.st = 0,
            anis = 0)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Variogram Modelling
sumMetric <- vgmST("sumMetric",
                   space = vgm(psill=3,"Sph", range=3, nugget=0.1),
                   time = vgm(psill=7,"Sph", range=3, nugget=0.1),
                   joint = vgm(psill=1,"Sph", range=2, nugget=0.1),
                   stAni=5)
sumMetric_vgm_sub <- fit.StVariogram(var_sub, 
                                 sumMetric, 
                                 method="L-BFGS-B",
                                 lower=pars.l)

ppp1 <- plot(var_sub, wireframe = T,
            zlab=NULL,
            xlab=list("distance (km)", rot=30, cex=0.7),
            ylab=list("time lag (days)", rot=-35, cex=0.7),
            scales=list(arrows=F, z = list(distance = 5), cex=0.5), 
            main="3d sample variogram")
ppp2 <- plot(var_sub, list(sumMetric_vgm_sub), wireframe = T, 
            xlab=list("distance (km)", rot=30, cex=0.7),
            ylab=list("time lag (days)", rot=-35, cex=0.7),
            scales=list(arrows=F, z = list(distance = 5), cex=0.5), main="3d fitted sum-metric variogram")
grid.arrange(ppp1, ppp2, ncol = 2, nrow = 1)

ppp3 <- plot(var_sub, map = F)
ppp4 <- plot(var_sub, list(sumMetric_vgm_sub), map = F)
grid.arrange(ppp3, ppp4, ncol = 2, nrow = 1)

grid.arrange(ppp1, ppp2, ppp3, ppp4, ncol = 2, nrow = 2)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Ordinary Kriging
ok_stfdf_2009 <- krigeST(as.formula(slp ~ 1),
                    data = temp_sub,
                    newdata = grid_1km_ok_sub,
                    modelList = sumMetric_vgm_sub,
                    computeVar = F)


dff <- na.omit(stfdf@data)

histogram <- ggplot(data=dff, aes(slp)) + 
  geom_histogram(binwidth = 1,
                 col="red", 
                 aes(y = ..density..,
                     fill = ..count..)) +
  scale_fill_gradient("Count", low = "blue", high = "red")+
  stat_function(fun = dnorm, 
                color = "black",
                size = 1.5,
                args = list(mean = mean(dff$slp), sd = sd(dff$slp)))+
  labs(title="Histogram of Sea level pressure") +
  labs(x="Sea level pressure [mb]", y="Count") + 
  annotate("text", x = 980, y=0.06, label = paste("Mean:",round(mean(dff$slp),digits = 2), "mb"))+
  annotate("text", x = 980, y=0.055, label = paste("SD:",round(sd(dff$slp),digits = 2), "mb"))+
  theme_bw()

histogram








