# Separable Variogram ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
separable <- vgmST("separable", 
                   space = vgm(psill = var_sp[[1]], "Exp", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget = 0.1),
                   time = vgm(psill = var_t[[1]], "Exp", range = 9, nugget = 0.1), 
                   sill=var_j[[1]])

plot(variogram_OK, separable, map=F) 


separable_Vgm <- fit.StVariogram(variogram_OK, separable)
attr(separable_Vgm,"MSE")

summary(stsdf_utm)

# There are several possibilities that lead to a singular covariance matrix. Two common ones:

#  duplicate observations (identical location & time stamp),
#  a variogram model that does not discriminate observations sufficiently, leading to near-perfectly correlated observations.

# Trazenje duplikata
dupl <- zerodist(stsdf_utm@sp) # find point pairs with equal spatial coordinates



k_var_ok_exp <- estiStAni(variogram_OK , c(0, 10), "metric", # Estimation of the spatio-temporal anisotropy without an underlying spatio-temporal model.
                          vgm(psill = var_sp[[1]], "Lin", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget = 0.1),
                          vgm(psill = var_t[[1]], "Exp", range = 9, nugget = 0.1))
pars.l_exp <- c(sill.s = 0, range.s = 1, nugget.s = 0, # Lower (and upper) bounds
                sill.t = 0, range.t = 1, nugget.t = 0,
                sill.st = 0, range.st = 1, nugget.st = 0,
                anis = 0)
sumMetric_var_ok_exp <- vgmST("sumMetric", # Variogram Modelling
                              space = vgm(psill = var_sp[[1]], "Lin", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget = 0.1),
                              time = vgm(psill = var_t[[1]], "Exp", range = 9, nugget = 0.1),
                              joint = vgm(psill = var_j[[1]], "Exp", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget=0.1),
                              stAni = k_var_ok_exp)
sumMetric_vgm_ok_exp <- fit.StVariogram(variogram_OK, # Fitting
                                        sumMetric_var_ok_exp, 
                                        method = "L-BFGS-B",
                                        lower = pars.l_exp)
attr(sumMetric_vgm_ok_exp,"MSE")

#year <- 2009
#time <- seq(as.Date(paste(year,"-01-01", sep = "")), as.Date(paste(year,"-12-31", sep="")), by="day")
#days=gsub("-", "", time, fixed = TRUE)

temp <- na.locf(stsdf_utm, na.rm = TRUE)

df_utm <- as.data.frame(stsdf_utm)



ggplot(data = df_utm, aes(x = timeIndex, y = slp))+
  geom_line(colour = "red", alpha = 0.5)+
  facet_wrap(~ staid, ncol = 8)+
  ggtitle("Vremenska linija merenja Normalnog vazdusnog pritiska po stanicama")+
  theme_bw()


df_utm <- df_utm[!(df_utm$staid == 13619), ]

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

k_var_ok_temp <- estiStAni(variogram_temp , c(0, 10), "metric", # Estimation of the spatio-temporal anisotropy without an underlying spatio-temporal model.
                           vgm(psill = var_sp[[1]], "Sph", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget = 0.1),
                           vgm(psill = var_t[[1]], "Sph", range = 9, nugget = 0.1))
pars.l_exp <- c(sill.s = 0, range.s = 1, nugget.s = 0, # Lower (and upper) bounds
                sill.t = 0, range.t = 1, nugget.t = 0,
                sill.st = 0, range.st = 1, nugget.st = 0,
                anis = 0)
sumMetric_var_ok_temp <- vgmST("sumMetric", # Variogram Modelling
                               space = vgm(psill = var_sp[[1]], "Sph", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget = 1),
                               time = vgm(psill = var_t[[1]], "Sph", range = 9, nugget = 0.1),
                               joint = vgm(psill = var_j[[1]], "Sph", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget=1),
                               stAni = 2)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Sa 64 stanica, bez stanice 13619 (koja nije imala opazanja)
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
temp <- stsdf_utm
temp@sp <- temp@sp[!(temp@sp$staid == 13619), ]
summary(temp)


summary(temp@data$slp)
hist(temp@data$slp)
all(is.na(temp@data$slp))

variogram_temp <- variogramST(slp ~ 1, temp)
a1111 <- plot(variogram_temp) # sample variograms at each time lag
b1111 <- plot(variogram_temp, map=FALSE) # ST-variogram map
c1111 <- plot(variogram_temp, wireframe=TRUE) # ST-variogram wireframe
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(a1111, b1111 ,c1111 , nrow = 2, ncol = 2)




estiStAni(variogram_temp , c(0, 10), "metric",
          vgm(psill=7,"Sph", range=180000, nugget=0),
          vgm(psill=7,"Sph", range=8, nugget=0) )

sumMetric_var_ok_temp <- vgmST("sumMetric",
                               space = vgm(psill=7,"Sph", range=180000, nugget=0.1),
                               time = vgm(psill=7,"Sph", range=8, nugget=0.1),
                               joint = vgm(psill=7,"Sph", range=180000, nugget=0.1),
                               stAni=5)


sumMetric_vgm_ok_temp <- fit.StVariogram(variogram_temp, # Fitting
                                         sumMetric_var_ok_temp, 
                                         method = "L-BFGS-B",
                                         lower = pars.l_exp)
attr(sumMetric_vgm_ok_temp,"MSE")


plot(variogram_temp, sumMetric_vgm_ok_temp, map=F) 
plot(variogram_temp, sumMetric_vgm_ok_temp, map=T) 

plot(variogram_temp, map=F)


obs = as(temp[,i_1[3]:ip1[3],'slp', drop=F],"STSDF")

x3 <- krigeST(as.formula("slp~1"),
              data=obs,
              newdata=STF(grid_1km,
                          temp@time[3],
                          temp@endTime[3]),  
              modelList = sumMetric_vgm_ok_temp,
              computeVar=FALSE)
hist(x3@data$var1.pred)
stplot(x3)


obs200 = as(temp[,i_1[200]:ip1[200],'slp', drop=F],"STSDF")

x200 <- krigeST(as.formula("slp~1"),
                data=obs200,
                newdata=STF(grid_1km,
                            temp@time[200],
                            temp@endTime[200]),  
                modelList = sumMetric_vgm_ok_temp,
                computeVar=FALSE)
hist(x3@data$var1.pred)
stplot(x200)





# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Dobar kod pre funkcije st_krige_prediction()
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


obs1 = as(stsdf_utm[,i_1[1]:ip1[1],'slp', drop=F],"STSDF")

x31 <- krigeST(as.formula("slp~1"),
               data=obs1,
               newdata=STF(grid_1km,
                           stsdf_utm@time[1],
                           stsdf_utm@endTime[1]),  
               modelList = sumMetric_vgm_ok,
               computeVar=FALSE)

#px31 <- stplot(x31)

obs3 = as(stsdf_utm[,i_1[3]:ip1[3],'slp', drop=F],"STSDF")

x33 <- krigeST(as.formula("slp~1"),
               data=obs3,
               newdata=STF(grid_1km,
                           stsdf_utm@time[3],
                           stsdf_utm@endTime[3]),  
               modelList = sumMetric_vgm_ok,
               computeVar=TRUE)

#grid.arrange(px31, px33, ncol = 2)

var1temp <- as.data.frame(x33) 
crs1 <- "+init=epsg:32634"
coords1 <- c("x", "y")
var1 <- st_as_sf(var1temp, crs = crs1, coords = coords1)

rast1 <- rasterize(var1, rast, var1$var1.pred, fun=mean) 
star_pred <- st_as_stars(rast1)

z1 <- ggplot()+
  geom_stars(data = star_pred)+
  scale_fill_viridis(option="B", limits = c(1010, 1030)) + # D = viridis color ramp, B = inferno, A = magma, C, E
  ggtitle(paste("Prediction plot", var1$time[1]))+
  xlab("Easting_UTM [m]")+
  ylab("Northing_UTM [m]")+
  labs(fill = "SLP [mb]")+
  theme_bw()+
  theme(legend.position="bottom")

rast2 <- rasterize(var1, rast, var1$var1.var, fun=mean) 
star_var <- st_as_stars(rast2)

z2 <- ggplot()+
  geom_stars(data = star_var)+
  scale_fill_viridis(option="D") + # D = viridis color ramp, B = inferno, A = magma, C, E
  ggtitle(paste("Variance plot", var1$time[1]))+
  xlab("Easting_UTM [m]")+
  ylab("Northing_UTM [m]")+
  labs(fill = "SLP [mb]")+
  theme_bw()+
  theme(legend.position="bottom")

grid.arrange(z1, z2, ncol = 2)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

obs200 = as(stsdf_utm[,i_1[200]:ip1[200],'slp', drop=F],"STSDF")

x200 <- krigeST(as.formula("slp~1"),
                data=obs200,
                newdata=STF(grid_1km,
                            stsdf_utm@time[200],
                            stsdf_utm@endTime[200]),  
                modelList = sumMetric_vgm_ok,
                computeVar=TRUE)


ok200temp <- as.data.frame(x200) 
crs1 <- "+init=epsg:32634"
coords1 <- c("x", "y")
ok200 <- st_as_sf(ok200temp, crs = crs1, coords = coords1)

rast200 <- rasterize(ok200, rast, ok200$var1.pred, fun = mean) 
star200_pred <- st_as_stars(rast200)

z3 <- ggplot()+
  geom_stars(data = star200_pred)+
  scale_fill_viridis(option = "B", limits = c(1010, 1030)) + # D = viridis color ramp, B = inferno, A = magma, C, E
  ggtitle(paste("Prediction plot", ok200$time[1]))+
  xlab("Easting_UTM [m]")+
  ylab("Northing_UTM [m]")+
  labs(fill = "SLP [mb]")+
  theme_bw()+
  theme(legend.position="bottom")

rast2200 <- rasterize(ok200, rast, ok200$var1.var, fun = mean) 
star200_var <- st_as_stars(rast2200)

z4 <- ggplot()+
  geom_stars(data = star200_var)+
  scale_fill_viridis(option="D") + # D = viridis color ramp, B = inferno, A = magma, C, E
  ggtitle(paste("Variance plot", ok200$time[1]))+
  xlab("Easting_UTM [m]")+
  ylab("Northing_UTM [m]")+
  labs(fill = "SLP [mb]")+
  theme_bw()+
  theme(legend.position="bottom")

grid.arrange(z3, z4, ncol = 2)


grid.arrange(z1, z2, z3, z4, ncol = 2, nrow = 2)


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Mapview
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ras1 <- (as(day_1[[3]], 'Raster'))
crs(ras1) <- crs("+init=epsg:32634")

ras2 <- (as(day_2[[3]], 'Raster'))
crs(ras2) <- crs("+init=epsg:32634")


pal <- colorRampPalette(brewer.pal(7, "PuBu"))

m1 <- mapview(day_1[[4]], layer.name = "day_1", col.regions = pal, query.type = c("mousemove"), query.digits = 2)
m2 <- mapview(day_2[[4]], layer.name = "day_2", col.regions = pal, query.type = c("mousemove"), query.digits = 2)

#sync(m1,m2)

#?mapview::sync
#install.packages('leafsync')
library(leafsync)
leafsync::sync(m1, m2, m3, m4, ncol = 2)




























