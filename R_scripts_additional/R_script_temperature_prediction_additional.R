
stsdf_temp <-  as(stfdf_temp, "STSDF")

d_temp <- na.omit(tempv)

d_temp <-  as.data.frame(stsdf_temp)
d_temp$time <- str_sub(d_temp$time, end=-1)
obs_temp <- d_temp
sta_temp <- as.data.frame(stsdf_temp@sp)

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

d_temp <- d_temp[!is.na(d_temp$tmax), ]
d_temp <- d_temp[!is.na(d_temp$tmin), ]
d_temp <- d_temp[!is.na(d_temp$tmean), ]

lim_lon <- range(d_temp$lon) # latitude range
lim_t <- range(d_temp$julian) # time range
lon_axis <- seq(lim_lon[1], # latitude axis
                lim_lon[2],
                length=25)

t_axis <- seq(lim_t[1], # time axis
              lim_t[2],
              length=100)
lon_t_grid <- expand.grid(lon = lon_axis,
                          t = t_axis)
tempvv_grid <- d_temp
dists <- abs(outer(d_temp$lon, lon_axis, "-"))
tempvv_grid$lon <- lon_axis[apply(dists, 1, which.min)]

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
