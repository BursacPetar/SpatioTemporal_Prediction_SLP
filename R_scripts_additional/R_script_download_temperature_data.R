library(spacetime)
library(rgdal)
library(zoo)
library(raster)
library(meteo)
library(doParallel)
library(R.utils)
library(plyr)
library(XML)
library(rvest)
library(xml2)

setwd("D:/Doktorske_akademske_studije_Geodezija_i_geoinformatika/Semestar_2/Prostorno_vremenska_statisitika/Temperature [Sale]/folder/")
getwd()

borders <- getData('GADM', country='SRB', level=0)

### OGIMET ############
### OGIMET download ###
### OGIMET stations download ###
countries = c("serbia", "hungary", "romania", "bulgaria", "kosovo", "macedonia", "albania",
              "montenegro", "croatia", "bosnia", "slovenia", "austria", "slovakia")

registerDoParallel(cores=detectCores()-1)
stations <- foreach (country_id = 1:length(countries), .packages = c("XML")) %dopar% {
  # for (country in countries){
  country = countries[country_id]
  url = paste("http://ftp.ogimet.com/display_stations.php?lang=en&tipo=AND&isyn=&oaci=&nombre=&estado=", country, "&Send=Send", sep="")
  html_table <- readHTMLTable(url, which=2)
  html_table$Latitude <- as.character(html_table$Latitude)
  html_table$Latitude <- ifelse(nchar(html_table$Latitude)==6,
                                as.integer(substr(html_table$Latitude, 1, 2)) +
                                  as.integer(substr(html_table$Latitude, 4, 5)) / 60,
                                as.integer(substr(html_table$Latitude, 1, 2)) +
                                  as.integer(substr(html_table$Latitude, 4, 5)) / 60 + 
                                  as.integer(substr(html_table$Latitude, 7, 8)) / 3600)
  html_table$Longitude <- as.character(html_table$Longitude)
  html_table$Longitude <- ifelse(nchar(html_table$Longitude)==7,
                                 as.integer(substr(html_table$Longitude, 1, 3)) +
                                   as.integer(substr(html_table$Longitude, 5, 6)) / 60,
                                 as.integer(substr(html_table$Longitude, 1, 3)) +
                                   as.integer(substr(html_table$Longitude, 5, 6)) / 60 + 
                                   as.integer(substr(html_table$Longitude, 8, 9)) / 3600)
  return(html_table)
  
}
stopImplicitCluster()

stations <- do.call("rbind", stations)
bbox(borders)
serbia = point.in.polygon(stations[,"Longitude"], stations[,"Latitude"], c(18,23.8,23.8,18),
                          c(41.4,41.4,47,47))
stations = stations[ serbia!=0, ] # 164
stations$`WMO INDEX` <- as.character(stations$`WMO INDEX`)
stations$ICAO <- as.character(stations$ICAO)
stations$Name <- as.character(stations$Name)
stations$Country <- as.character(stations$Country)
stations$Altitude <- as.numeric(as.character(stations$Altitude))
stations$Notes <- as.character(stations$Notes)
stations <- stations[stations$`WMO INDEX` != "----" & stations$`WMO INDEX` != "-----", ]

stations = stations[!is.na(stations$Longitude),]
stations <- stations[, c(1,6,5,3,4,7)]
names(stations) <- c("staid", "lon", "lat", "name","country", "altitude")
stations <- stations[!duplicated(stations$staid), ]
stations <- stations[stations$staid %in% unique(observations$staid), ]

save(stations, file = "ogimet_st.rda")
load(file = "ogimet_st.rda")


### OGIMET observations download ###
years = 2009:2018

time<-seq(as.Date(paste(as.numeric(min(years)),"-01-02", sep="")), as.Date(paste(max(years)+1,"-01-01", sep="")), by="day")
days<-gsub("-","",time,fixed=TRUE)
daysNum <- length(time)
stNum <- nrow(stations)

registerDoParallel(cores=detectCores()-1)
observations <- foreach (d = 1:daysNum, .packages = c("XML", "xml2", "rvest", "spacetime", "rgdal", "zoo", "meteo", "R.utils", "plyr")) %dopar% {
  day <- time[d]
  obs <- matrix(NA, nrow = stNum, ncol = 5)
  for (st in 1:stNum) {
    wmo <- stations[st, "staid"]
    if (format(day,"%m")=="01" & format(day,"%d")=="01") {
      url = paste("http://ftp.ogimet.com/cgi-bin/gsynres?lang=en&ind=", wmo, "&ndays=1&ano=", substr(day-1, 1, 4), "&mes=", substr(day-1, 6, 7), "&day=", substr(day-1, 9, 10), "&hora=23&ord=REV&Send=Send&Send=Send", sep="")
    } else {
      url = paste("http://ftp.ogimet.com/cgi-bin/gsynres?lang=en&ind=", wmo, "&ndays=1&ano=", substr(day, 1, 4), "&mes=", substr(day, 6, 7), "&day=", substr(day, 9, 10), "&hora=00&ord=REV&Send=Send&Send=Send", sep="")
    }
    file<-read_html(url)
    tables<-html_nodes(file, "table")
    # html_table <- readHTMLTable(url, which=3)
    if (length(tables) == 3) {
      html_table <- html_table(tables[3], fill = TRUE)[[1]]
      if (length(html_table) != 1 & (sum(names(html_table)=="Temperature(C)", na.rm = T)>0)) { # || sum(names(html_table)=="Prec.(mm)", na.rm = T)>0)) {
        # if (grepl("Gust", as.character(html_table[2, 6]))) {
        obs_st <- c(wmo, as.character(day))
        if (sum(names(html_table)=="Temperature(C)", na.rm = T)==3) {
          obs_st <- c(obs_st, c(as.numeric(html_table[2, 2]),
                                as.numeric(html_table[2, 3]),
                                as.numeric(html_table[2, 4])))
        } else {
          obs_st <- c(obs_st, c(NA, NA, NA))
        }
        # if (any(names(html_table)=="Prec.(mm)")) {
        #   obs_st[6] <- (as.numeric(html_table[, "Prec.(mm)"][2]))
        # } else {
        #   obs_st[6] <- NA
        # }
        obs[st, ] <- obs_st
      } else {
        next
      }
    } else {
      next
    }
  }
  return(obs)
}
stopImplicitCluster()
observations <- as.data.frame(do.call("rbind", observations))
names(observations) <- c("wmo", "date", "tmax", "tmin", "tmean") # , "prcp")

#### ------------> NAPRAVI KOPIJU !!!


observations <- observations[!(is.na(observations$tmax) & is.na(observations$tmin) & is.na(observations$tmean)), ]# & is.na(observations$prcp)), ]
observations$wmo <- as.character(observations$wmo)
observations$date <- as.Date(as.character(observations$date)) - 1
observations$tmax <- as.numeric(as.character(observations$tmax))
observations$tmin <- as.numeric(as.character(observations$tmin))
observations$tmean <- as.numeric(as.character(observations$tmean))
# observations$prcp <- as.numeric(as.character(observations$prcp))

# observations[is.na(observations$tmax), ]
# observations[observations$tmax < observations$tmin, ]
# observations[observations$tmax < observations$tmean, ]
# observations[observations$tmean < observations$tmin, ]
# summary(observations)

save(observations, file = "ogimet_obs.rda")
load("ogimet_obs.rda")
load("ogimet_st.rda")
### stfdf creation ###
stfdf <- meteo2STFDF ( obs      = observations,
                       stations = stations,
                       crs      = CRS("+proj=longlat +datum=WGS84"),
                       obs.staid.time=c(1,2),
                       stations.staid.lon.lat=c(1,2,3)
)

# remove duplicates
stfdf = rm.dupl(stfdf, zcol = 1, zero.tol = 5) # 5 km
stfdf_temp <- stfdf

plot(stfdf@sp)
plot(borders, add=T)
save(stfdf_temp, file="ogimet_temp_stfdf.rda")
summary(stfdf)


load("ogimet_temp_stfdf.rda")

