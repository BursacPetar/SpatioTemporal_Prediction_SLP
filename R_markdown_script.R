#' ---
#' title: "Spatio-Temporal prediction of Sea level pressure using Geostatistics and Machine learning methods"
#' author: "Student: Petar Bursac; Professor: Milan Kilibarda"
#' date: "`r format(Sys.time(), '%d %B %Y')`"
#' output:   
#'    html_document:
#'      theme: "simplex"
#'      highlight: tango
#'      toc: true
#'      toc_depth: 5
#'      toc_float: true
#'      fig_caption: yes
#' ---
#'<style>
#'body {
#'  text-align: justify}
#'.html-widget {
#'margin: auto;
#'}
#'</style>
#'
#'  <img src="Grb_Gradjevinski.png" align="center" alt="logo" width="150" height = "150" style = "border: none; fixed: right;">
#'  <img src="R_logo.png" align="center" alt="logo" width="150" height = "150" style = "border: none; fixed: right;">
#' 
#'
#+ include = TRUE, echo = FALSE, results = 'hide', warning = FALSE, message = FALSE
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
#' 
#' 
#' 
#' 
#+ include = FALSE, echo = FALSE
wdir <- 'd:/Doktorske_akademske_studije_Geodezija_i_geoinformatika/Semestar_2/Prostorno_vremenska_statisitika/Seminarski_zadatak/'
setwd(wdir)
getwd()

mycolors=c("#f32440","#2185ef","#d421ef")
"%!in%" = Negate("%in%")

# Plot settings
my_theme <- function(base_size = 10, base_family = "sans"){
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      axis.title = element_text(size = 13),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#fffcfc"),
      strip.background = element_rect(fill = "#820000", color = "#820000", size =0.5),
      strip.text = element_text(face = "bold", size = 10, color = "white"),
      legend.position = "bottom",
      legend.justification = "center",
      legend.background = element_blank(),
      panel.border = element_rect(color = "grey30", fill = NA, size = 0.5)
    )
}
theme_set(my_theme())

#'
#'
#' # **Uvod**
#'      
#' Cilj ovog rada je analiza i predikcija dnevnih karata Normalnog vazdusnog pritiska (Sea level pressure) na podrucju teritorije R. Srbije, na osnovu merenja dostupnih sa 65 mernih stanica.
#' Prostorni domen predikcije definisan je pravilnim gridom cija prostorna dispozicija odgovara prostornom rasporedu mernih stanica, dok temporalni domen merenja obuhvata period od 10 godina (2009-2018),
#' sa dostupnim merenjima na dnevnom nivou izrazenih kao srednje vrednosti dnevnog prikupljanja podataka sa senzora. U radu je koriscen programski jezik R sa dostupnim biblioteka kojima su implementirane tehnike Prostorno vremenske predikcije (gstat, GSIF, randomForest, caret, ranger).
#'   
#' 
#' # **Metodologija**
#'      
#' Metodologija podrazumeva jasno definisan niz postupaka u cilju izrade karata Normalnog vazdusnog pritiska sa prostornom rezolucijom od 1km.
#' Metodologija obuhvata koriscenje sledecih pristupa u oceni i predikciji posmatrane promenljive:      
#' 
#' - Oridnary kriging - Obicni kriging       
#' - Regression kriging - Regresioni kriging (Ocena trenda linearnom regresijom, kvanitifikacija i modelovanje prostorne korelacije na osnovu reziduala i predikcija vrednosti reziduala)       
#' - Machine learning techique Random Forest - Koricenje tehnika Masinskog ucenja, tj. Random Forest klasifikatora u cilju ocene trenda, a zatim koriscenje Kriging tehnika u cilju izrade karata predikcije.       
#'
#' Kao ocena kvaliteta predikcije izvrsen je postupak krosvalidacije na dva nacina:
#' 
#' - Leave-one station out - sukcesivno iskljucivanje pojedinacnih stanica iz modela i ponavljanje procesa predikcije
#' - 5-fold neasted cross-validation - ugnježdena krosvalidacija sa 5 foldova, od kojih u svakom narednom koraku jedan fold se koristi za izgradnju modela (trening), dok se preostalih 4 koriste za ocenu i predikciju [Pejovic et al., 2018]
#'
#' # **Analiza i istraživanje seta podataka**
#' 
#' Analizom dostupnog seta podataka na raspolaganju su sledeći podaci:
#' 
#' - 65 mernih stanica sa vrednostima koordinata tacaka u WGS84 referentnom koordinatnom sistemu
#' - Vremenska serija podataka koja odgovara periodu od 01.01.2009 - 31.12.1018 godine (10 godina - 3652 dana). Kao jedinica mere definisana je [dan]
#' - Vrednosti posmatrane promenljive SLP - Sea Level Pressure u jedinicama mere [mb], gde je ukupan broj vremenskih trenutaka na svim stanicama bez vrednosti opazanja 23527.
#' 

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
load(file = "Ulazni_podaci/ogimet_st.rda") # stations - stanice
load(file = "Ulazni_podaci/ogimet_slp_stfdf.rda") # stfdf - STFDF sa pritiskom (slp) za 10 godina (2009-2018)
load(file = "Ulazni_podaci/ogimet_obs.rda") # observations - merenja

# ===================================================================================
# Data analysis, preparation and wrangling
# ===================================================================================
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
po_drzavama <- stations %>% group_by(country) %>% # Broj stanica po drzavama
  summarise(n = n())
po_drzavama %>%
  kable(caption = "Broj mernih stanica po drzavama", digits = 4) %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
#' ### Prostorna dispozicija mernih stanica vazdušnog pritiska
#' Za potrebe analize seta podataka i istraživanja (sagledavanja lokacija mernih stanica, njihove udaljenosti i prostornog rasporeda, visina na kojima se nalaze), 
#' veoma je bitno u početnoj fazi rada sagledati prostornu dispoziciju seta podataka, kao i postaviti početne hipoteze u vezi njihovih lokacija. Na sledećem sledećoj karti moguće je sagledati prostornu dispoziciju 65 mernih stanica. 
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
dat=map_data("world")
crs <- "+init=epsg:4326"
coords <- c("lon", "lat")
stations_sf <- st_as_sf(stations, crs = crs, coords = coords)
stations_sf %<>% dplyr::rename("Visina [m]" = altitude)
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 8, fig.height = 8, fig.align='center'
ggplot(stations_sf)+
  geom_sf(aes(colour=`Visina [m]`, size=`Visina [m]`))+
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
#' Prostorna dispozicija mernih stanica Normalnog vazdusnog pritiska data je i putem web kartografskog prikaza, koji je kreiran koristeći R programski paket "mapview" [Appelhans et al., 2018]
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 6, fig.height = 6, fig.align='center'
mapview(stations_sf, z = "Visina [m]", layer.name = "Visina mernih stanica [m]")

#' ### Vremenska serija podataka
#' 
#' Analizom seta podataka, kreiran je grafički prikaz vremenske serije podataka. Na gornjem delu grafičkog prikaza data je vremenska linija vrednosti Normalnog vazdušnog pritiska, po godinama, 
#' sa tim što su po $\mathbf{x}$ osi date vrednosti vremenskih trenutaka izraženih po Julianskom kalendaru (godina 2455000 po Julianskom kalendaru odgovara 2009 godini po Gregorijanskom kalendaru),
#' Na donjem delu grafika dat je boxplot grafički prikaz Normalnog vazdušnog pritiska izraženog po godinama.
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
df_slp <- as.data.frame(stfdf)
df_slp <- df_slp[!is.na(df_slp$slp), ]
df_slp$time <- str_sub(df_slp$time, end=-1)
obs <- df_slp
sta <- stations
obs$date_1 <- obs$time
obs %<>% separate(date_1, "-", into = c("y", "m", "d"))
slpr <- obs
slpr %<>% mutate(y = as.numeric(y),
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
slpr$julian <- toJulian(slpr$y, slpr$m, slpr$d)
slpr$t <- slpr$julian - 2454832
slpr_dan <- slpr %>% 
  dplyr::group_by(y, m, d) %>%
  summarise(slp_dan = mean(slp), 
            julian = mean(julian))
slpr_dan$vreme <- paste(slpr_dan$y, "-", slpr_dan$m, "-", slpr_dan$d)
ggp1 <- ggplot()+
  geom_line(data = slpr_dan, aes(x = julian, y = slp_dan, group = y, colour = y))+
  scale_colour_gradient(low="blue", high="red")+
  xlab("Julian [2455000 = 2010 godina]") + ylab("Normalni vazdusni pritisak [mb]") +
  ggtitle("Vremenska linija merenja Normalnog vazdusnog pritiska")+
  theme_bw()
ggp3 <- ggplot(slpr_dan, aes(x=y, y=slp_dan, group = y, colour =y)) +
  geom_boxplot()+
  scale_colour_gradient(low="blue", high="red")+
  xlab("Godina") + ylab("Normalni vazdusni pritisak [mb]") +
  ggtitle("Normalni vazdusni pritisak po godinama")+
  theme_bw()
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(ggp1, ggp3, nrow = 2, ncol = 1)

#' Analizom grafičkih prikaza ne uočavaju se velike promenu u intenzitetu pritiska, ako se vrši poređenje po godinama. 
#' Zaključuje se da su vrednosti u potpunosti slučajnog karaktera, iz statističkog ugla posmatranja.
#' Najmanja izmerena vrednost zabeležena je 2015. godine.     
#' 
#' Na sledećem grafiku je moguće ostvariti uvid u vremensku liniju podataka po stanicama:
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
ggplot(data = slpr, aes(x = julian, y = slp))+
  geom_line(colour = "red", alpha = 0.5)+
  facet_wrap(~ staid, ncol = 8)+
  ggtitle("Vremenska linija merenja Normalnog vazdusnog pritiska po stanicama")+
  theme_bw()

#' Analizom podataka i uvidom u grafički prikaz vremesnkih serija podataka po stanicama, jasno se uočava da na pojedinim stanicama postoje vremenski trenuci bez registrovanih vrednosti vazdušnog pritiska. 
#' Stanica 13619 je u potpunosti bez merenja za dati vremenski period. U kasnijem nastavku rada dati vremenski trenuci bez merenja, kao i stanica, su uklonjeni iz seta podataka.        
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
#' ### Histogram
#' 
#' Na sledećem grafičkom prikazu dat je histogram frekvencija merenih vrednosti vazdušnog pritiska, kao i kriva normalne rapodele koja odgovara datom uzorku.
#' Uočava se da se podaci u potpunosti pokoravaju normalnoj raspodeli, sa srednjom vrednošću i merama odstupanja vrednosti slučajne promenjive od njenog matematičkog očekivanja - standardnom devijacijom.
#' Kriva raspodele je blago pomerena, i postoje uzorci u setu podataka čiji intenzitet je izvan matematčkog očekivanja. Kasnijom obradom podataka potrebno je ostvariti detaljniji uvid u set podataka kao i izvršiti statisitčko testiranje hipoteza da u skupu podataka nema statistički značajnih grubih grešaka.
#' Glavni preduslov za korišćenje tehnika predikcije - kriginga je da podaci, tj. uzorak merenja pripada normalnoj raspodeli, što je i ispunjeno.  
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
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
  labs(title="Histogram merenja Normalnog vazdusnog pritiska") +
  labs(x="Normalni vazdusni pritisak [mb]", y="Count") + 
  annotate("text", x = 980, y=0.06, label = paste("Mean:",round(mean(dff$slp),digits = 2), "mb"))+
  annotate("text", x = 980, y=0.055, label = paste("SD:",round(sd(dff$slp),digits = 2), "mb"))+
  theme_bw()
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 6, fig.height = 6, fig.align='center'
histogram
#' ### Hovmoller grafički prikaz
#' 
#' Hovmöller plot - grafički prikaz je dvodimenzionalni prostorno-vremenski način vizuelizacije. Prostorni domen je predstavljen na jednoj osi (projektovan ili osrednjen), dok je po drugoj osi definisano vreme.
#' Ovakav način vizuelizacije je tradicionalno primenjivan u naukama i oblastima koje se bave praćenjem stanja atmosfere, mora i okeana, kako bi se na veoma efikasan način predstavila - vizuelizovala propagacije entiteta. [Wikle et al., 2019]
#' Pomoću ovakvog načina vuzuelizacije moguće je sagledati korelacije u prostorno-vremenskom domenu, tj. npr. geografske širine ili dužine i vremena, vremenskih trenutaka opažanja posmatrane promenljive. 
#' Na ovaj način moguće je ostvariti i preliminarni uvid u postajenje trenda koji se opisuje i koordinatama tačaka.
#' 
#' 
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
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
Hovmoller_lat <- ggplot(slpr_lat_Hov) + # take data
  geom_tile(aes(x = lat, y = julian, fill = slp)) + # plot
  fill_scale(name = "Normalni vazdusni pritisak [mb]") + # add color scale
  scale_y_reverse() + # rev y scale
  ylab("Dan") + # add y label
  xlab("Latitude (deg)") + # add x label
  ggtitle("Hovmoller plot - latitude")+
  theme_bw()+
  theme(legend.position="bottom")
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
  theme_bw()+
  theme(legend.position="bottom")
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(Hovmoller_lat, Hovmoller_lon, nrow = 1, ncol = 2)
#'
#' Analzom Hovmöller plot-ova po koordinatnim osama definisanog referentnog koordinatanog sistema (WGS84), kao i vremena - vremenskih trenutaka opažanja,
#' uočava se da ne postoji znatna prostorno-vremenska korelacija, koja bi opisivala korelisanost Normalnog vazdušnog pritiska opažanog tokom vremena, kao i lokacija stanica (njihove geografske dužine i širine).
#'
#' ### Empirijska prostorna sredina (*Empirical Spatial Mean*)
#' 
#' Veoma koristan način analize prostorno-vrmenskih podatka je računanje empirijske sredine. Ako se uzme u obzir da su dostupna opažanja u definisanom prostornom i vremenskom domenu, empirijska prostorna sredina za svaku pojedinu lokaciju - mernu stanicu je sračunata agregacijom vrednosti tokom vremena.
#' Empirijsku prostornu sredinu je moguće vizuelizovati, kao što je prikazano na sledećem scatter plot-u (i Howmoller plot je jedan od načina vizuelizacije). Cilj ove analize je pre svega kvantifikacija prostorno-vrmenske zavisnosti, korelisanosti, kao i ocene i vrednovanja potencijalnog trenda u podacima.
#' 
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
spat_av <- group_by(slpr, lat, lon) %>% # group by lon-lat
  summarise(spm_emp = mean(slp)) # mean for each lon-lat
gg1 <- ggscatter(spat_av, x = "lon", y = "spm_emp",
                 color = "red", alpha = 0.5,
                 add = "loess",  # Add regressin line ("reg.line") or local regression fitting ("loess")
                 add.params = list(color = "blue", fill = "orange"), # Customize reg. line
                 conf.int = TRUE)+ # Add confidence interval
  stat_cor(method = "pearson")+# Add correlation coefficient
  xlab("Longitude (deg)") +
  ylab("Average sea level pressure (mb)") +
  ggtitle("Empirical Spatial Mean - longitude")
gg2 <- ggscatter(spat_av, x = "lat", y = "spm_emp",
                 color = "red", alpha = 0.5,
                 add = "loess",
                 add.params = list(color = "blue", fill = "orange"),
                 conf.int = TRUE)+
  stat_cor(method = "pearson")+
  xlab("Latitude (deg)") +
  ylab("Average sea level pressure (mb)") +
  ggtitle("Empirical Spatial Mean - latitude")
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 6, fig.align='center'
grid.arrange(gg1, gg2, ncol = 2, nrow = 1)
#'
#' Na datom grafičkom prikazu vrednosti pritiska u odnosu na geografsku širinu i dužinu, koje su dobijene agregacijom tokom vremena, 
#' jasno se uočava da postoji zavisnost između geografske širine i vazdušnog pritiska. 
#' Veća vrednost koeficijenta determinacije ovo potvđuje, ali još uvek nedovoljno, čime se postavlja hipoteza za dalja istraživanja i statisitičko testiranje.
#'
#' ### Emprijska vremenska sredina (*Empirical Temporal Means*)
#'
#' Empirijska vremenska sredina predstavlja agregaciju vrednosti promenljive na lokacijama od interesa tokom vremena, sa tim što sad agregacija se vrši u vremenskom domenu koji je posmatran, za razliku od empirijske prostorne sredine gde su lokacije mernih stanica bile od interesa.
#' Agregacija podataka je sračunata tokom godina, meseca i dana, a dostupna je na sledećem grafičkom prikazu.
#'
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
slpr_dan <- slpr %>% 
  dplyr::group_by(y, m, d) %>%
  summarise(slp_dan = mean(slp), 
            julian = mean(julian))
slp_av_d <- group_by(slpr,d) %>% # Po danu
  summarise(meanSlp = mean(slp))
em_mean_d <- ggplot() +
  geom_line(data = slpr, aes(x = d, y = slp, group = staid),
            colour = "red", alpha = 0.04) +
  geom_line(data = slp_av_d, aes(x = d, y = meanSlp)) +
  geom_boxplot(data = slpr_dan, aes(x=d, y=slp_dan, group = d, colour =d), show.legend = FALSE)+
  xlab("Dan") + ylab("Normalni vazdusni pritisak [mb]") +
  ggtitle("Emprijska vremenska sredina - po danu")+
  theme_bw()
slp_av_m <- group_by(slpr,m) %>% # Po mesecu
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
slp_av_y <- group_by(slpr,y) %>% # Po godini
  summarise(meanSlp = mean(slp))
em_mean_y <- ggplot() +
  geom_line(data = slpr, aes(x = y, y = slp, group = staid),
            colour = "red", alpha = 0.04) +
  geom_line(data = slp_av_y, aes(x = y, y = meanSlp)) +
  geom_boxplot(data = slpr_dan, aes(x=y, y=slp_dan, group = y, colour =y), show.legend = FALSE) +
  xlab("Godina") + ylab("Normalni vazdusni pritisak [mb]") +
  ggtitle("Emprijska vremenska sredina - po godini") +
  theme_bw()
#' Graficki prikaz koji sadrzi srednje vrednosti normalnog vazdusnog pritiska iskazanom po danu, mesecu i godini
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(em_mean_d, em_mean_m, em_mean_y, ncol = 1, nrow = 3)
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
#'
#' Uvidom u vremesnke serije podataka, osrednjene tokom definisanih vremenskih intervala, kao i prostorno sumiranih, jasno se uočava pojava paterna i uticaja doba godine (seasonal pattern) i drugih vremenskih prilika i atmosferskih promenljivih na vazdušni pritisak.
#' Uočava se periodična promena raspona vrednosti na par godina, ali godišnje srednje vrednosti su približno jednakih vrednosti.
#' 
#' ### Prosta linearna regresija (*Simple Linear Regression*)
#' 
#' Cilj proste linearne regresije je opisivanje i kvatifikovanje veze izmedju dve neprekidne slucajne velicine. 
#' Ako je:
#' 
#' * $\mathbf{y_i}$ - vrednost zavisne promenljive i-te opservacije u uzorku
#' * $\mathbf{x_i}$ - vrednost prediktora (nezavisne promenljive) i-te opservacije u uzorku
#' 
#' Dati uzorak neće pripadati tačno lineranoj pravoj, ali dovoljno je da bude blizu sa nekim manjim odstupanjem, pa da linearna veza bude dovoljno objašnjenje. Dakle, model je:     
#' 
#' $$ \begin{aligned} \mathbf{y_i=b_0+b_1x_i+e_i} \end{aligned} $$      
#' 
#' Parameter $\mathbf{b_0}$ se naziva intercept - presek, a $\mathbf{b_1}$ slope - nagib. Ovaj model se naziva prosta linearna regresija (Simple Linear Regression). 
#' Ako je u pitanju linearna zavisnost jedne zavisne sa više nezavisnih promenljivih radi se o višestrukoj linearnoj regresiji (Multiple Linear Regression). 
#' Iako se koriste nazivi zavisna i nezavisna promenljiva, vrednosti $\mathbf{x_i}$ se ne smatraju slučajnim u modelu. Slučajan je niz grešaka $\mathbf{e_i}$, pa odatle i niz $\mathbf{y_i}$.
#' 
#' Ako je dat dvodimenzioni uzorak (x,y) sa grafika raspršenja (scatter plot) dobija se informacija o postojanju i obliku zavisnosti.
#' 
#' Scatter plot sa koeficijentima korelacije:
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
sp1 <- ggscatter(slpr, x = "lon", y = "slp",
                 color = "red", alpha = 0.05,
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE)+ # Add confidence interval
  stat_cor(method = "pearson")# Add correlation coefficient
sp2 <- ggscatter(slpr, x = "lat", y = "slp",
                 color = "red", alpha = 0.05,
                 add = "reg.line",
                 add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE)+
  stat_cor(method = "pearson")
sp3 <- ggscatter(slpr, x = "julian", y = "slp",
                 color = "red", alpha = 0.05,
                 add = "reg.line",
                 add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE)+
  stat_cor(method = "pearson")
slpr$julian2 <- slpr$julian^2
sp4 <- ggscatter(slpr, x = ("julian2"), y = "slp",
                 color = "red", alpha = 0.05,
                 add = "reg.line",
                 add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE)+
  stat_cor(method = "pearson")
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(sp1, sp2, sp3, sp4, nrow = 2, ncol = 2)
#' 
#' 
#' Kod linearne regresije parametri modela su odabrani tako da je suma kvadrata reziduala najmanja. Postavlja se pitanje i koliko je ta minimalna suma mala, jer ona govori o tome u 
#' kojoj meri je linerani model odgovarajući, koliko je disperzije početnog skupa objašnjeno uzimanjem u obzir linerane veze sa prediktorom. Jedna mera koja govori o tome koliko je smanjeno
#' srednje kvadratno odstupanje je koeficijent determinacije $\mathbf{R^2}$. Koficijent determinacije je izražen na sledeći način:
#' 
#' $$ \begin{aligned} \mathbf{R^2 = RSS/SSTO = 1 - SSE/SSTO} \end{aligned} $$
#' 
#' gde je:
#' 
#' * SSR - Regression Sum of Squares - suma kvadrata reziduala regresije
#' * SSTO - Total Sum of Squares - ukupna suma kvadrata reziduala
#' * SSE - Error Sum of Squares - suma grešaka kvadrata reziduala    
#' 
#' Ako je $\mathbf{R^2 = 1}$ to znači da sve tačke $\mathbf{y_i}$ pripadaju regresionoj pravoj. Ako je $\mathbf{R^2 = 0}$ onda je ocenjena linija horizontalna, pa prediktor ne objašnjava varijabilnost u zavisnoj promenljivoj ništa bolje od srednje vrednosti.
#' Nije moguće na osnovu $\mathbf{R^2}$ porediti modele koji imaju i one koji nemaju intercept. Naziv $\mathbf{R^2}$ potiče od toga što je ovako definisana vrednost jednaka +/− korenu koeficijenta korelacije između prediktora i zavisne promenljive.
#' Vrednost $\mathbf{R^2}$ se može često zloupotrebiti i pogrešno shvatiti, što je opisano sledećim primerima: 
#' 
#' * $\mathbf{R^2}$ meri linearnu vezu između datih promenljivih. Ako se dobije da je $\mathbf{R^2 = 0}$ to znači da nema linearne zavisnosti, ali možda postoji neka nelinearna zavisnost. 
#' * Velika vrednost $\mathbf{R^2}$ ne mora da znači da je linearna funkcija najbolja veza. Ovo se može videti sa grafika. 
#' * Jedna neobična vrednost u uzorku može puno da utiče na vrednost $\mathbf{R^2}$. 
#' * Velika vrednost $\mathbf{R^2}$ ne znači prediktivnu moć.
#' 
#' ### Grid - domen predikcije 
#' 
#' Za potrebe kreiranja predikcije vrednosti Normalnog vazdušnog pritiska nad područjem od interesa i u vremenskom domenu, koji je definisan kao dan, 
#' kreiran je grid - skup tačaka na pravilnom i jednakom rastojanju. Prostorna rezolucija grid-a je 1km, dok je vremenska rezolucija 1 dan, posmatrano za vremsnki interval 01.01.2009-31.12.2018.
#' Ciljni grid je kreiran kao objekat STF klase podataka, paketa "spacetime" (Spatio-Temporal Feature). 
#' 
#'  
#' # **Prostorno-vremenski obični kriging (*Spatio-temporal Ordinary Kriging*)**    
#' 
#' U okviru ovog poglavlja dat je opis, metodologija i rezultati predikcije Normalnog vazdušnog pritiska kao dnevnih karata na području od interesa, obuhvaćenog mernim stanicama. 
#' Metodologija se zasniva na korišćenju Kriging tehnika predikcije u prostorno-vremenskom domenu.
#' 
#' Sva merenja – opažanja na površi Zemlje, iako je ponekad ignorisano, imaju prostorno –
#' vremensku komponentu (referencu). Prostorno vremenska komponenta opažnja je određena sa
#' minimum četri parametra:
#' 
#' - Geografska – prostorna lokacija (geografskim koordinatama ili u definisanom koordinatnom sistemu u projekciji).
#' - Visina na samoj površi Zemlje u odnosu na definisanu referentnu nivosku površ (visina).
#' - Vreme merenja – opažnja (godina, mesec, dan, sat, itd.).
#' - Prostorno – vremensku podršku – određenost (veličinu bloka pridruženih merenjima, vremenski interval merenja;).
#' 
#' Na osnovu prethodno definisanih komponenti moguće je definisati čitavi niz metodologije obrade i analize podataka u mnogim
#' oblastima, a zajedničko za tu grupu nauka je da se analize definšu kao
#' prostorno – vremenske analize podataka (*Spatio – temporal Data Analysis*). Grubo rečeno prostorno
#' vremenske analize podataka se mogu posmatrati kao kombinacija dve
#' osnovne nauke: geoinformatika i prostorno vremenska statistika. Druga reč koja se može upotrebiti za termin
#' prostorne statistike, a isto tako dodeliti i vremensku dimenziju, je nauka koja se zove Geostatistika (*Geostatistics*).  
#' 
#' 
#' Geostatistika je grana matematičke statistike koja se bavi analizom i interpretacijom
#' geografski referenciranih podataka, tj. prostorno određenih i kojima je moguće dodeliti vremensku dimenziju.
#' Geostatistika pre svega ima fokus na prostonim podacima koji su kontinualni (prostorna
#' kontinualnost) i pritom se postavljaju osnovna pitanja na koja je potrebno da odgovori: *kako se promenljiva ponaša u prostorno – vremenskom koordinatnom sistemu, šta je uzrok
#' i kontrola varijacije, gde je potrebno utorkovati određenu pojavu da bi se opisala prostorna promenljivost, koliko uzoraka je potrebno, kolika je vrednost promenljive na nekoj novoj lokaciji, kolika je nesigurnost ocenjenih vrednosti.*       
#' Postavljena pitanja opisiju glavni zadatak na koji Geostatistika kao nauka je potrebno da odgovori. Sama metodologija je jasno definisana i osnovne
#' hipoteze koje se postavljaju pred svaku promenljivu koja egzistira u prostoru su sledeće:
#' 
#' - Određivanje modela (*Model estimation*) – zaključci o parametrima modela.
#' - Predikcija – interpolacija vrednosti promenljive (*Prediction*) – zaključivanje o prediktovanoj vrednosti posmatrane promenljive.
#' - Testiranje hipoteza (*Hypothesis testing*).      
#' 
#' Kao sinonim za Geostatističku interpolaciju može se reći da je to **Kriging**, kao statističku
#' prostornu predikciju. Geostatističke metode su prvi razradili francuski geomatematičar Georges
#' Matheron i južnoafricki rudarski inženjer D.G. Krige. Zato se ove metode u literaturi nazivaju
#' kriging. U geodeziju je geostatistiku uveo austrijski geodeta Moritz (pod nazivom kolokacija) za
#' potrebe odredivanja anomalija sile zemljine teže i otklona vertikala.       
#' U oviru ovog rada date su osnovne teroijske postavke na kojima se zaniva Kriging kao
#' metoda predikcije prostornih promenljivih. Dat je osvrt na kvantifikaciju prostorne strukture
#' promenljivih, definisanje variograma i komponenti variograma, kao i koncepti Kriging
#' predikcije.
#' Eksperimentalni deo rada je urađen dostupnim programskim paketima, kojima je moguće izvršiti geostatističku analizu i
#' interpretaciju u okviru programskog jezika R. Osnovne funkcionalnosti su dostupne putem paketa "gstat", koji podržava i rad sa prostorno vremenskim promenljivima, čija implementacija je moguća kroz definisane klase paketa "spacetime":
#' 
#' - STF - *Spatio-Temporal Full grid layout* - Prostorno - vremenski potpuni grid 
#' - STS - *Spatio-Temporal Sparse grid layout* - Prostorno - vremenski redak grid 
#' - STI - *Spatio-Temporal Irregular layout* - Prostorno - vremenski nepravilan grid 
#' - STT - *Spatio-Temporal Trajectory* - Prostorno - vremenska trajektorija
#' 
#' Svaka od pomenutih klasa omogućuje i proširenje na -DF - Data Frame u cilju čuvanja osobina entiteta, atributa. Forma u kojoj je podaci mogu čuvati je:
#'
#' - *long formats* - potpuna kombinacija prostornog i vremenskog domena
#' - *space-wide* - gde različite kolone predstavljaju različite lokacije u prostoru
#' - *time-wide* - gde različite kolone različite trenutke vremena
#' 
#' Evidentno da definisani pristup čuvanju i obradi podataka putem dostupnih paketa, u poptunosti predstavlja kompatibilno i efikasno rešenje za potrebe prostorno-vremenske predikcije koristeći Geostatističke metode.
#' 
#' ### Variogram - modelovanje i kvantifikacija prostorne korelisanosti
#' 
#' Interpolacija primenom geostatističkih metoda se izvodi u dva koraka:
#' - Kvantifikacija prostorne strukture površi koja se modelira (na osnovu ulaznih podataka), i
#' ▪ Predikcije, tj. ocene vrednosti funkcije površi u zadatim tačkama.
#' 
#' Kvantifikacija prostorne strukture površi koja se modelira izvodi se empirijskim
#' određivanjem kovarijacione funkcije, odnosno poluvariograma kojima se opisuje prostorna
#' zavisnost vrednosti funkcije u tačkama površi. Osnovni problem kod praktične primene ovih
#' metoda je upravo izbor odgovarajuće kovarijacione funkcije (funkcije poluvariograma), tj. izbor
#' odgovarajućeg modela funkcije i empirijsko određivanje parametara te funkcije.
#' Varijacije vrednosti promenljive koja se posmatra, najbolje se mogu opisati stohastičkim
#' površima, gde se posmatrana pojava posmatra kao regionalizovana promenljiva. Vrednost
#' promenljive $\mathbf{Z}$ na lokaciji $\mathbf{x}$ se tada posmatra kao realizacija jedne slučajne promenljive, $\mathbf{Z(x)}$ je
#' slučajna funkcija, a cela površ slučajni proces. Teorija regionalizovane promenljive pretpostavlja da se prostorne varijacije mogu
#' razložiti na tri nkomponente:
#' 
#' $$ \begin{aligned} \mathbf{Z(x) = m(x) + ε’(x) + ε’’(x)} \end{aligned} $$
#' 
#' - determinističke varijacije (trend) – $\mathbf{m(x)}$
#' - slučajne, ali prostorno autokorelisane varijacije – $\mathbf{ε’(x)}$
#' - prostorno nekorelisani šum (očekivana vrednost šuma je $\mathbf{0}$, a disperzija je $\mathbf{σ^2)}$ - $\mathbf{ε’’(x)}$
#'
#' Za slučajne promenljive u prostoru se kaže da su regionalizovane. Za ekstrakciju prostornih
#' osobina regionalizovanih promenljivih možemo koristiti geostatističke mere. Jednom
#' kvantifikovane, osobine regionalizovanih promenljivih mogu se koristiti u mnogo aplikacija.
#' Geostatistička analiza koristi prostorne autokorelacione informacije u procesu kriging
#' interpolacije. *Fenomen da su geografski bliže pojave međusobno više korelisane nego one koje
#' su dalje.* Autokorelacija je statistička veza između merenih tačaka, čija korelacija zavisi od
#' udaljenja i pravca pojedinačnih lokacija. Mi znamo iz posmatranja realnog sveta da prostorna
#' autokorelacija postoji zato što generalno zapažamo da pojave koje su međusobno bliže su mnogo
#' više slične nego one koje su udaljenije. Kako se udaljenje povećava, prostorna autokorelacija
#' opada.      
#' Variografija (*Variography*) je proces kvantifikacije modela prostorne zavisnosti koja
#' proizilazi iz podataka i prostorne strukture. Za određivanje veličine interpolacije na specifičnoj
#' lokaciji, kriging koristi postavljen model u variografiji, prostornu konfiguraciju podataka i
#' vrednosti merenog uzorka tačaka oko lokacije za predikciju.
#' Jedna od najvažnijih mera koja se koristi za razumevanje prostorne strukture
#' regionalizovanih promenljivih je **variogram**, koji se može korsititi za odnos varijanse unutar dela
#' prostora (i autokorelacije) između uzorka. Poluvarijansa daje nepristrasan opis razmere i
#' obrasca/šablona/uzorka prostorne promenljive unutar regiona.
#' Geostatistički model zahteva kvantitativno merenje prostorne korelacije (*spatial corelation*) radi
#' ocene i simulacija i naviše korišćen alat kojim se meri, kvantifikuje prostorna korelacija je
#' **poluvariogram**.
#' Variogram je grafikon kojim se rastojanja prikazuju u formi korelacije. Na osonovu seta ulaznih
#' podataka, merenja iz realnog sveta moguće je dobiti eksperimentalni variogram. Osnovne
#' karakteristike su da na kraćim rastojanjima korelacija je veća, dok povećanjem rastojanja
#' korelacija opada. Na određenim rastojanjima među podacima ne postoji korelacija.
#'
#' Dve bitne pretpostvke je potrebno uzeti u obzir:
#' 
#' - γ – poluvarijansa je funkcija rastojanja – pretpostvka **stacionarnosti** (*Stationarity assumption*).
#' - Funkcija samo rastojanja, ali ne i pravca – pretpostvka **izotropije** (*Isotropy assumption*).      
#' 
#' U slučaju da pravac ulazi u obzir, da prostorne razlike u određenim pravcima naglašene u
#' odnosu na druge pravce, dolazi do pojave anizotropije (*anistropy assumption*). Ako postoji način
#' moguće je odkloniti trendom. U slučaju pojave anizotropije definišu se azimuti i zatim definišu
#' direkcioni variogram ili kompletni rasterski variogram (zonalni).
#'
#' Definisana terminologija i pretpostavke su u potpunosti proširive u domenu vremena, tako da sada u slučaju Prosotorno-vremenskih podataka, potrebno je i vreme uzeti u obzir. 
#' Definicija poluvarijanse se dopunjuje u pogledu vrmemna, tako da sada poluvariogram dobija treću dimenziju kojom se opisuje i međusobna zavisnost podataka u vremenu, a ne samo u prostoru. 
#' U narednom delu biće detaljnije objašnjen koncept prostorno-vremenske zavisnosti, kao i način kvantifikacije i izrade empirijskog prostorno-vremenskog poluvariograma na osnovu merenja vazdušnog pritiska.
#'
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
stsdf <- as(stfdf, "STSDF")
#'
#' #### Transformacija podataka u koordinatni sistem u projekciji
#' 
#' Koordinate lokacija mernih stanica su izražene u elipsoidnim koordinatama. U cilju izrade variograma i kvantifikacije prostorne korelisanosti, potrebno je da rastojanja budu izražena u metrima, a ne u stepenima. U cilju toga izvršena je
#' transformacija koordinata mernih stanica iz geografskog-elipsoidnog koordinatnog sistema WGS84 u koordinatni sistem u projekciji (UTM - Universal Mercator Projection - 34N zona), koji je dat sledećom notacijom u obliku EPSG koda (EPSG:32634):
#' 
#' $$ \begin{aligned} \mathbf{+init=epsg:32634 +proj=utm +zone=34 +datum=WGS84} \end{aligned} $$
#' $$ \begin{aligned} \mathbf{+units=m +nodefs +ellps=WGS84 +towgs84=0,0,0} \end{aligned} $$
#' 
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
stsdf_utm <- spTransform(stsdf, CRS("+init=epsg:32634"))
stsdf_utm@sp <- stsdf_utm@sp[!(stsdf_utm@sp$staid == 13619), ] # brisanje stanice bez merenja
#' #### Variogram 
#' Na sledećem grafičkom prikazu dat je prostorno-vremenski variogram opažanja vazdušnog pritiska:
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE

variogram_OK <- variogramST(slp ~ 1, stsdf_utm)
a1 <- plot(variogram_OK) # sample variograms at each time lag
b1 <- plot(variogram_OK, map=FALSE) # ST-variogram map
c1 <- plot(variogram_OK, wireframe=TRUE) # ST-variogram wireframe
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(a1, b1 ,c1 , nrow = 2, ncol = 2)
#'                
#' ### Empirijski prostorno-vremenski variogram       
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
#' Nakon kreiranja prostorno-vremenskog variograma, koji opisuje set podataka i zavisnost - korelaciju između merenja na različitim lokacijama i različitim vremenskim trenucima, potrebno je izvršiti fitovanje modela - određivanje matematičke funcije i njenih parametara koji na najbolji mogući način opisuju set podataka.
#' Fitovanje variograma se postiže određivanjem inicijalnih parametara funkcije na osnovu vizuelnog pregleda podataka. Konačna empirijska funkcija, kao i njeni parametri se određuju automatski kroz više iteracija.
#' R programski paket "gstat" nudi više modela za fitovanje variograma prostorno-vremenskih podataka: **separable, product sum, metric, sum metric, and simple sum metric**. Više informacija o svakom modelu ponaosob se može pronaću u radu [Graler et al., 2015].
#' Za modelovanje prostorno-vremenskog variograma vrednosti merenja vazdušnog pritiska, korišćn je *sum metric* model, koji podrazumeva uključivanje prostorne i vremenske komponente putem modela kovarijansi - poluvariograma, kao i zajedničke komponente sa vrednovanjam parametra anizotropije.
#' Prostorno vremenska funkcija kovarijanse data je na sledeći način: 
#' $$ \begin{aligned} \mathbf{Csm (h, u) = Cs (h) + Ct (u) + Cjoint (sqrt(h^2 + (k*u)^2))} \end{aligned} $$
#'
#'gde je:    
#'
#' - $\mathbf{h}$ - rastojanje u prostornom domenu [m]
#' - $\mathbf{u}$ - rastojanje u vremenskom domenu [dan]
#' - $\mathbf{C}$ - matrica kovarijansi (s - prostorne, t - vremenske i joint - prostorno-vremenske komponente)
#' - $\mathbf{k}$ - parametar anizotropije
#'
#'
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
#' U cilju određivanja inicijalnih vrednosti parametara variograma potrebno je zadati početne vrednosti parametra psill - parcijalni prag i range - domet.
#' Vrednost parcijalnog praga je približno jednaka vrednosti varijanse seta podataka.
#' Za početnu vrednost dometa može uzeti 1/4 vrednosti dijagonale minimalnog obuhvatnog pravougaonika područja od interesa obuhvaćenog stanicama, sa tim da se ta vrednost uzima samo za prostornu komponentu.      
#'         
#' Srednja vrednost varijanse vrednosti Normalnog vazdušnog pritiska na osnovu vrednosti po stanicama [mb]
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
df_utm <- as.data.frame(stsdf_utm)
var_loc <- df_utm %>% 
  group_by(staid) %>%
  summarise(var = var(slp))
var_loc %<>% mutate(var = coalesce(var, 0)) # Replace NA with 0
var_sp <- mean(var_loc$var)
var_sp %<>% as.data.frame() %>% 
  mutate(units = "mb") %>% 
  rename(Variance_s = ".")
var_sp %>%
  kable(caption = "Varijansa merenja u prostornom domenu", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)
#' Srednja vrednost varijanse vrednosti Normalnog vazdušnog pritiska na svim stanicama za posmatrani vremenski period [mb]
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
var_t <- var(stsdf_utm@data$slp)
var_t %<>% as.data.frame() %>% 
  mutate(units = "mb") %>% 
  rename(Variance_t = ".")
var_t %>%
  kable(caption = "Varijansa merenja u vremenskom domenu", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)
#' Srednja vrednost varijanse vrednosti Normalnog vazdušnog pritiska u prostorno-vremenskom domenu [mb]
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
var_j <- (var_sp[[1]] + var_t[[1]])/2
var_j %<>% as.data.frame() %>% 
  mutate(units = "mb") %>% 
  rename(Variance_j = ".")
var_j %>%
  kable(caption = "Varijansa merenja u prostorno-vremenskom domenu", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
k_var_OK <- estiStAni(variogram_OK , c(0, 10), "metric", # Estimation of the spatio-temporal anisotropy without an underlying spatio-temporal model.
                      vgm(psill = var_sp[[1]]/2, "Sph", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget = 0),
                      vgm(psill = var_t[[1]]/2, "Sph", range = 9, nugget = 0))
pars.l <- c(sill.s = 0, range.s = 1, nugget.s = 0, # Lower (and upper) bounds
            sill.t = 0, range.t = 1, nugget.t = 0,
            sill.st = 0, range.st = 1, nugget.st = 0,
            anis = 0)
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
sumMetric_var_OK <- vgmST("sumMetric", # Variogram Modelling
                   space = vgm(psill = var_sp[[1]]/2, "Sph", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget = 0.1),
                   time = vgm(psill = var_t[[1]]/2, "Sph", range = 9, nugget = 0.1),
                   joint = vgm(psill = var_j[[1]]/2, "Sph", range = sqrt(diff(stsdf_utm@sp@bbox[1,])^2 + diff(stsdf_utm@sp@bbox[2,])^2)/4, nugget=0.1),
                   stAni = k_var_OK)
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
sumMetric_vgm_OK <- fit.StVariogram(variogram_OK, # Fitting
                                     sumMetric_var_OK, 
                                     method = "L-BFGS-B",
                                     lower = pars.l)
ppp1 <- plot(variogram_OK, wireframe = T, # Plot
             zlab=NULL,
             xlab=list("Distance [m]", rot=30, cex=0.7),
             ylab=list("Time lag [day]", rot=-35, cex=0.7),
             scales=list(arrows=F, z = list(distance = 5), cex=0.5), 
             main="3D variogram")
ppp2 <- plot(variogram_OK, list(sumMetric_vgm_OK), wireframe = T, 
             xlab=list("Distance [m]", rot=30, cex=0.7),
             ylab=list("Time lag [day]", rot=-35, cex=0.7),
             scales=list(arrows=F, z = list(distance = 5), cex=0.5), main="3D fitovan sum-metric variogram")
#grid.arrange(ppp1, ppp2, ncol = 2, nrow = 1)
ppp3 <- plot(variogram_OK, 
             map = F)
ppp4 <- plot(variogram_OK, 
             list(sumMetric_vgm_OK), 
             map = F)
#grid.arrange(ppp3, ppp4, ncol = 2, nrow = 1)
#' Prikaz empirijskog prostorno-vremenskog variograma, fitovanog korišćenjem *sum-metric* modela. Model je dobijen korišćenjem funkcija vgmST() i fit.StVariogram(), R programskog paketa "gstat". 
#' Približna vrednost koeficijenta anizotropije je sračunata empirijski korišćenjem funkcije estiStAni().
#' Za fitovanje modela variograma korišćena je sferna funkcija kod sve tri komponente: prostorne, vremenske i zajedničke - prostorno-vremenske.
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(ppp1, ppp2, ppp3, ppp4, ncol = 2, nrow = 2)
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
#' Mean Squared Error (MSE) - srednja kvadratna greška optimizacije modela variograma po metodi najmanjih kvadrata
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
mse_ok <- attr(sumMetric_vgm_OK, "MSE")
mse_ok %<>% as.data.frame() %>% 
  mutate(units = "mb^2") %>% 
  rename("Mean Squared Error" = ".")
mse_ok %>%
  kable(caption = "Srednja kvadratna greska fitovanja variograma", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)
# sumMetric_vgm_OK$stAni/1000 # /24/3600
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
spaceC <- as.data.frame(sumMetric_vgm_OK$space[1:3])
timeC <- as.data.frame(sumMetric_vgm_OK$time[1:3])
jointC <- as.data.frame(sumMetric_vgm_OK$joint[1:3])
anis <- as.data.frame(sumMetric_vgm_OK$stAni)
anis <- anis %>% 
  rename(anisotropy = `sumMetric_vgm_OK$stAni`) %>%
  mutate(units = "km/day")

spaceC %>%
  kable(caption = "Prostorna komponenta", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)
timeC %>%
  kable(caption = "Vremenska komponenta", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)
jointC %>%
  kable(caption = "Prostorno-vremenska komponenta", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)
anis %>%
  kable(caption = "Koeficijent anizotropije", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)
#'
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
#' ### Predikcija
#' 
#' Predikcija nad setom podataka u cilju dobijanja dnevnih karata Normalnog vazdušnog pritiska je urađena korišćenjem Običnog prostorno-vremenskog kriginga. Kriging tehnike su implementirane u okviru R programskog paketa "gstat", uključijući i rad sa
#' podacima u formatu podataka definisanim paketom "spacetime". Prostorni domen predikcije je definisan gridom sa rastojanjem tačaka od 1km, dok je vremenski domen predikcije definisan vremenskim periodom dostupnog seta podataka, za period od 2009-2018 godine.
#'    
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
daysNum <- length(stsdf_utm@time)
i_1=c(rep(1,1),1:(daysNum -1)) ### 2 days ###
ip1=c(1:daysNum)
dem_1km_utm34 <- raster::stack("d:/Doktorske_akademske_studije_Geodezija_i_geoinformatika/Semestar_2/Prostorno_vremenska_statisitika/Seminarski_zadatak/R [scripts]/Grid/Grid_1km_UTM34.tif")
grid_1km <- as(dem_1km_utm34, "SpatialPointsDataFrame") # Make grid from dem 1km, as new locations for prediction
proj4string(grid_1km) <- proj4string(stsdf_utm@sp)
grid_1km@data <- grid_1km@data[-1]
grid_1km <- as(grid_1km,"SpatialPoints")
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
#xxx <- lapply(1:daysNum, function(i) {
#  
#  obs = as(stsdf_utm[,i_1[i]:ip1[i],'slp', drop=F],"STSDF")
#  
#  krigeST(as.formula("slp~1"),
#          data=obs,
#          newdata=STF(grid_1km,
#                      stsdf_utm@time[i],
#                      stsdf_utm@endTime[i]),  
#          modelList = sumMetric_vgm_OK,
#          computeVar=FALSE)
#} )
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
v1 <- stsdf_utm@time$timeIndex['2009-01-01']
v2 <- stsdf_utm@time$timeIndex['2010-01-01']
v3 <- stsdf_utm@time$timeIndex['2011-01-01']
v4 <- stsdf_utm@time$timeIndex['2012-01-01']
v5 <- stsdf_utm@time$timeIndex['2013-01-01']
v6 <- stsdf_utm@time$timeIndex['2014-01-01']
v7 <- stsdf_utm@time$timeIndex['2015-01-01']
v8 <- stsdf_utm@time$timeIndex['2016-01-01']
v9 <- stsdf_utm@time$timeIndex['2017-01-01']
v10 <- stsdf_utm@time$timeIndex['2018-01-01']
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
rast <- raster()
extent(rast) <- extent(dem_1km_utm34) # this might be unnecessary
ncol(rast) <- ncol(dem_1km_utm34) # this is one way of assigning cell size / resolution
nrow(rast) <- nrow(dem_1km_utm34)
crs(rast) <- crs("+init=epsg:32634")
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
st_krige_prediction <- function(day_number = day_number, data = data, variogram_ok = variogram_ok, grid1km = grid1km, compute_variance = TRUE, raster_temp = raster_temp){
  daysNum <- length(data@time)
  i_1=c(rep(1,1),1:(daysNum -1))
  ip1=c(1:daysNum)
  
  obs = as(data[,i_1[day_number]:ip1[day_number],'slp', drop=F],"STSDF")
  
  krige_result <- krigeST(as.formula("slp~1"),
                 data = obs,
                 newdata = STF(grid1km,
                             data@time[day_number],
                             data@endTime[day_number]),  
                 modelList = variogram_ok,
                 computeVar = compute_variance)
  
  crs <- "+init=epsg:32634"
  coords <- c("x", "y")
  
  oktemp <- as.data.frame(krige_result)
  ok_sf <- st_as_sf(oktemp, 
                    crs = crs, 
                    coords = coords)
  
  rast_pred <- rasterize(ok_sf, raster_temp, ok_sf$var1.pred, fun = mean) 
  star_pred <- st_as_stars(rast_pred)
  
  rast_vari <- rasterize(ok_sf, raster_temp, ok_sf$var1.var, fun = mean) 
  star_vari <- st_as_stars(rast_vari)
  
  pred_plot <- ggplot()+
    geom_stars(data = star_pred)+
    scale_fill_distiller(palette = "PuBu") +
    ggtitle(paste("Prediction plot", ok_sf$time[1]))+
    xlab("Easting_UTM [m]")+
    ylab("Northing_UTM [m]")+
    labs(fill = "SLP [mb]")+
    theme_bw()+
    theme(legend.position="bottom")
  
  vari_plot <- ggplot()+
    geom_stars(data = star_vari)+
    scale_fill_viridis(option="D") + # D = viridis color ramp, B = inferno, A = magma, C, E
    ggtitle(paste("Variance plot", ok_sf$time[1]))+
    xlab("Easting_UTM [m]")+
    ylab("Northing_UTM [m]")+
    labs(fill = "SLP [mb]")+
    theme_bw()+
    theme(legend.position="bottom")
  
  raster_day <- (as(star_pred, 'Raster'))
  
  pal <- colorRampPalette(brewer.pal(7, "PuBu"))
  mapview_day <- mapview(raster_day, layer.name = paste(ok_sf$time[1], " SLP [mb]"), col.regions = pal, query.type = c("mousemove"), query.digits = 2)
  
  list_results <- list(pred_plot, vari_plot, star_pred, mapview_day)
  return(list_results)
}
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
day_1 <- st_krige_prediction(day_number = v1[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

day_2 <- st_krige_prediction(day_number = v2[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

day_3 <- st_krige_prediction(day_number = v3[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

day_4 <- st_krige_prediction(day_number = v4[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

day_5 <- st_krige_prediction(day_number = v5[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

day_6 <- st_krige_prediction(day_number = v6[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

day_7 <- st_krige_prediction(day_number = v7[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

day_8 <- st_krige_prediction(day_number = v8[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

day_9 <- st_krige_prediction(day_number = v9[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

day_10 <- st_krige_prediction(day_number = v10[[1]], data = stsdf_utm, variogram_ok = sumMetric_vgm_OK, grid1km = grid_1km, compute_variance = TRUE, raster_temp = rast)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_1[[1]], day_1[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_2[[1]], day_2[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_3[[1]], day_3[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_4[[1]], day_4[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_5[[1]], day_5[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_6[[1]], day_6[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_7[[1]], day_7[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_8[[1]], day_8[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_9[[1]], day_9[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 8, fig.align='center'
grid.arrange(day_10[[1]], day_10[[2]], ncol = 2, nrow = 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.align='center'
leafsync::sync(day_1[[4]], 
               day_2[[4]], 
               day_3[[4]], 
               day_4[[4]], 
               day_5[[4]], 
               day_6[[4]], 
               day_7[[4]], 
               day_8[[4]], 
               day_9[[4]], 
               day_10[[4]], 
               ncol = 2)
#'      
#'      
#' ### Krosvalidacija 
#' Ocenu kvaliteta modela i predikcije kontinualne promenljive je moguće oceniti unutrašnjim i spoljašnjim merama kvaliteta. Kao unutrašnje mere kvaliteta, koriste se koeficijent determinacije lineranog modela i 
#' varijansa kriging predikcije. Kod spoljašnjih mera ocena kvaliteta modela predikcije, može se koristiti nezavisan set podataka kojim se model ocenjuje tokom ili nakon izgradnje modela. Kada nije dostupan set podataka za validaciju modela, 
#' moguće je iskoristiti isti set podataka, ali isključivanjem pojedinih prostornih i/ili vremesnkih entiteta - merenja iz izgradnje modela. Takva metoda ocene kvaliteta i validacije modela predikcije naziva se $\mathbf{krosvalidacija}$.      
#' Krosvalidacijom se vrši ocena modela predikcije na taj način što se za izgradnju modela koristi potpun set podataka, a zatim prilikom predikcije u zadatom prostorno-vremenskom domenu jedan i/ili više merenja na jednoj ili više lokacija uzoraka se ne koriste.
#' Koriste se sledeće metode krosvalidacije:      
#'
#' - $\mathbf{LOOCV}$ - Leave-one-out cross validation - isključivanje jedne stanice ili entiteta u vremenskom domenu tokom predikcije.
#' - $\mathbf{5-FOLD}$ - Krosvalidacija sa isključivanjem 5 merenja u prostornom i/ili vremenskom domenu tokom predikcije ili deljenjem seta podataka na 5 grupa - foldova od kojih se jedan koristi za validaciju seta podataka, dok se preostali koriste za izgradnju modela.
#' - $\mathbf{5-FOLD-neasted-CV}$ - Ugnježdena krosvalidacija sa 5 foldova, gde se podrazumeva kreiranje 5 grupa podataka, gde prilikom svake iteracije se jedna od grupa koristi za validaciju a preostale za izgradnju modela, i tako sukcesivno u što više nivoa i kreiranja podgrupa (sub-folds).
#' 
#' Date metode krosvalidacije nisu primenljive za modele koji ne podrazumevaju teoriju geostatistike, tj. teoriju prostorne korelacije. U ovom radu je korišćena prva metoda krosvalidacije - LOOCV u prostornom domenu sa pojedinačnim isključivanjem stanica iz predikcije. Korišćen je identičan model - prostorno-vrmenski variogram dobijen na osnovu potpunog skupa merenja. 
#' Sa obzirom da su dostupne 64 merne stanice sa opažanjima od 10 godina - ukupno 3652 dana, dobijeno je ukupno 233728 (64x3652) karata predikcije. Na osnovu dobijenih podataka sračunate su mere kvaliteta krosvalidacije: koeficijent determinacije i ukupna srednje kvadratna greška predikcije.
#' 
#'
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
daysNum <- length(stsdf_utm@time) # ----- >>>   daysNum = 6
sp.nmax = 20 # umesto da se koriste sve stanice prilikom predikcije jedne tacke, koristi se samo lokalno okruzenje od npr 20 tacaka
i_1 = c(1, 1:(daysNum -1))
ip1 = 1:daysNum
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
temp = stsdf_utm
nrowsp <- length(temp@sp)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
numNA <- apply(matrix(temp@data[,'slp'], # count NAs per stations
                      nrow=nrowsp,byrow=F), MARGIN=1,
               FUN=function(x) sum(is.na(x)))
time <- stsdf_utm@time # Remove stations out of covariates # ----- >>>   time <- stsdf_utm@time[1:6]
rem <- numNA != length(time)
temp <-  temp[rem,drop=F]
N_POINTS <- length(temp@sp@coords[,1])

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
registerDoParallel(cores=detectCores()-1)
cv <- foreach(i = 1:N_POINTS, .packages = c("raster","spacetime","gstat","rgdal","meteo")) %dopar% {
  
  st= temp@sp
  st$dist=spDists(temp@sp,temp@sp[i,])
  tmp_st<-st[order(st$'dist') , ]
  ### remove target station. Instead of sp.nmax all of the stations can be used, but accuracy stays the same ###
  local_t= row.names(tmp_st[2:sp.nmax,] )
  
  xxx = as.list ( rep(NA, length(time) ) )
  for( ii in 1:length(time) ) {
    data=temp[local_t, i_1[ii]:ip1[ii],'slp',drop=F]
    nrowsp <- length(data@sp)
    # count NAs per stations
    numNA <- apply(matrix(data@data[,'slp'],
                          nrow=nrowsp,byrow=F), MARGIN=1,
                   FUN=function(x) sum(is.na(x)))
    # Remove stations out of covariates
    rem <- !numNA > 0
    data <- data[rem,drop=F]
    
    xxx[[ii]] = krigeST(as.formula("slp~1"),
                      data = data, 
                      newdata = STF(as(temp@sp[i,],"SpatialPoints"),
                                  temp@time[ii],  
                                  temp@endTime[ii]),     
                      modelList = sumMetric_vgm_OK,
                      computeVar = FALSE)@data[,1]
  } # end of  for
  ret=unlist(xxx) 
}
stopImplicitCluster()
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
cv <- do.call(rbind,cv)
cv <- as.vector(cv)
cv.slp_ok <- temp
cv.slp_ok <- as(cv.slp_ok, 'STFDF')
cv.slp_ok$pred.cv <- cv 
cv.slp_ok$resid.cv <- cv.slp_ok@data$slp - cv.slp_ok@data$pred.cv

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
rmse_total = sqrt(sum((cv.slp_ok$resid.cv)^2, na.rm = T)/(length(cv.slp_ok$resid.cv[!is.na(cv.slp_ok$resid.cv)]))) # 1.4916
#' Srednja kvadratna greška - $\mathbf{RMSE}$ (Root Mean Square Error) predikcije, dobijena krosvalidacijom nad setom podataka, metodom LOOCV:
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
rmse_total %<>% as.data.frame() %>% 
  mutate(units = "mb") %>% 
  dplyr::rename(RMSE = ".")
rmse_total %>%
  kable(caption = "Srednje kvadratna greška - RMSE", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
data = cv.slp_ok@data[complete.cases(cv.slp_ok@data[, c(1,3)]),]
tss = t(data$slp - mean(data$slp)) %*% (data$slp - mean(data$slp))
ess = t(data$resid.cv) %*% (data$resid.cv)
r2 = (tss-ess)/tss #0.9575941
#' Koeficijent determinacije $\mathbf{R^2}$ predikcije, dobijen krosvalidacijom nad setom podataka, metodom LOOCV:
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
r2 %<>% as.data.frame() %>% 
  mutate(units = "/") %>% 
  dplyr::rename(R2 = "V1")
r2 %>%
  kable(caption = "Koeficijent determinacije - R2", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
month_results <- matrix(nrow = 12, ncol = 5) ### per month ###
for (m in 1:12){
  all <- cv.slp_ok[, month(cv.slp_ok@time)==m ]
  rmse = sqrt(sum((all$resid.cv)^2, na.rm = T)/(length(all$resid.cv[!is.na(all$resid.cv)])))
  max = max(all$slp, na.rm = T)
  min = min(all$slp, na.rm = T)
  range = max - min
  mean <- mean(all$slp, na.rm = T)
  #print(round(c(rmse, max, min, mean, range),2))
  month_results[m, ] <- round(c(rmse, max, min, mean, range),2)
}
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
month_results <- as.data.frame(month_results)
month_results <- dplyr::rename(month_results, `RMSE [mb]` = "V1", `max` = "V2", `min` = "V3", `mean` = "V4", `range` = "V5")
month_results$month <- 1:12
month_results$Mesec <- c("Jan", "Feb", "Mar", "Apr", "Maj", "Jun", "Jul", "Avg", "Sep", "Okt", "Nov", "Dec")
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
mo_res_plot <- ggplot(month_results, aes(x = Mesec, y= mean, group = Mesec, fill = `RMSE [mb]`)) +
  geom_boxplot(aes(ymax = max+rmse, ymin = min-rmse, middle = mean, lower = min, upper = max), stat = "identity")+
  scale_fill_gradient(low="blue", high="red")+
  xlab("Mesec") + ylab("Normalni vazdusni pritisak [mb]") +
  ggtitle("Krosvalidacija [Leave-one-station out]")+
  scale_x_discrete(limits=month_results$Mesec)+
  theme_bw()+
  theme(legend.position="bottom")
#' Analizom rezultata krosvalidacije u vremenskom domenu, nakon agregacije u prostornom domenu, moguće je ostvariti uvid u ocenu kvaliteta kao i opseg vrednosti po mesecima. 
#' Očigledno je da postoji zakonitost u podacima, kojom se opisuje manji opseg vrednosti, kao i greška predikcije u letnjim mesecima. 
#' U zimskim mesecima postoji veća varijabilnost seta podataka, kao i da je greška predikcije veća.    
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 10, fig.height = 4, fig.align='center'
mo_res_plot
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
res_month_cv <- month_results %>%
  dplyr::rename(MESEC = Mesec, `MAX [mb]` = max, `MIN [mb]` = min, `MEAN [mb]` = mean, `RANGE [mb]` = range)
res_month_cv %<>% dplyr::select(-month) 
res_month_cv %<>% dplyr::select(MESEC, `RMSE [mb]`, `MAX [mb]`, `MIN [mb]`, `MEAN [mb]`, `RANGE [mb]`) 
res_month_cv %>%
  kable(caption = "Rezultati krosvalidacije u vremenskom domenu po mesecima:", digits = 4, align = "c") %>%
  kable_styling(bootstrap_options = "striped", full_width = TRUE)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
stNum = length(cv.slp_ok@sp)
results = matrix(nrow = stNum, ncol = 3)
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
for (st in 1:stNum){
  rmse = sqrt(sum((cv.slp_ok[st, ]$resid.cv)^2, na.rm = T)/(length(cv.slp_ok[st, ]$resid.cv[!is.na(cv.slp_ok[st, ]$resid.cv)])))
  max = max(cv.slp_ok[st, ]$resid.cv, na.rm = T)
  min = min(cv.slp_ok[st, ]$resid.cv, na.rm = T)
  #h = cv.slp_ok[st, , drop=F]@sp$h
  #dem = cv.temp_gl[st, , drop=F]@sp$dem
  results[st,] = c(rmse, max, min)#, h, dem)
}
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
colnames(results) = c("rmse", "max", "min")
results_order = results[order(results[,1]),]
cv.slp_ok@sp$rmse = results[,1]
cv.slp_ok@sp$max = results[,2]
cv.slp_ok@sp$min = results[,3]
global_loo_plot = cv.slp_ok@sp[!is.na(cv.slp_ok@sp$rmse),]
global_loo_plot$rmse = round(global_loo_plot$rmse, 1)
global_loo_plot$max = round(global_loo_plot$max, 1)
global_loo_plot$min = round(global_loo_plot$min, 1)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
global_loo_plot = global_loo_plot[order(global_loo_plot@data$rmse, decreasing = T), ]
global_loo_plot = global_loo_plot[order(global_loo_plot@data$max, decreasing = T), ]
global_loo_plot = global_loo_plot[order(global_loo_plot@data$min, decreasing = T), ]

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
bubbles_global_loo = bubbleSP(global_loo_plot, zcol="rmse", scale_e = 100)
bubbles_global_loo_max = bubbleSP(global_loo_plot, zcol="max", scale_e = 100)
bubbles_global_loo_min = bubbleSP(global_loo_plot, zcol="min", scale_e = 100)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
bubble_sf_rmse <- st_as_sf(bubbles_global_loo)
bubble_sf_max <- st_as_sf(bubbles_global_loo_max)
bubble_sf_min <- st_as_sf(bubbles_global_loo_min)

#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
bubble_sf_rmse %<>% dplyr::rename("RMSE [mb]" = rmse) 
bubble_sf_rmse <- st_transform(bubble_sf_rmse, crs="+init=epsg:4326")
bubble_sf_max %<>% dplyr::rename("MAX [mb]" = max) 
bubble_sf_max <- st_transform(bubble_sf_max, crs="+init=epsg:4326")
bubble_sf_min %<>% dplyr::rename("MIN [mb]" = min) 
bubble_sf_min <- st_transform(bubble_sf_min, crs="+init=epsg:4326")

#' Agregacijom seta podataka u prostornom domenu, moguće je ostvariti uvid u grešku predikcije po stanicama, na kojima je opažan vazdušni pritisak za dati vremenski period.
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
ok_plot_rmse <- ggplot(bubble_sf_rmse)+
  geom_sf(aes( fill = `RMSE [mb]`))+ #colour = "red",
  guides(colour=FALSE, size=FALSE)+
  scale_fill_distiller(palette = "OrRd", direction = 1)+
  labs(x = "Longitude [deg]",
       y = "Latitude [deg]",
       title = "RMSE vrednosti na stanicama nakon krosvalidacije",
       caption = "Coordinate Reference System - WGS84")+
  geom_map(data=dat[dat$region=="Serbia",], map=dat[dat$region=="Serbia",],
           aes(x=long, y=lat, map_id=region),
           color="white", fill="#7f7f7f", size=0.05, alpha=1/4) +
  theme_bw()+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) 
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
ok_plot_max <- ggplot(bubble_sf_max)+
  geom_sf(aes( fill = `MAX [mb]`))+ #colour = "red",
  guides(colour=FALSE, size=FALSE)+
  scale_fill_distiller(palette = "OrRd", direction = 1)+
  labs(x = "Longitude [deg]",
       y = "Latitude [deg]",
       title = "Maksimalne vrednosti razlika na stanicama nakon krosvalidacije",
       caption = "Coordinate Reference System - WGS84")+
  geom_map(data=dat[dat$region=="Serbia",], map=dat[dat$region=="Serbia",],
           aes(x=long, y=lat, map_id=region),
           color="white", fill="#7f7f7f", size=0.05, alpha=1/4) +
  theme_bw()+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) 
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE
ok_plot_min <- ggplot(bubble_sf_min)+
  geom_sf(aes( fill = `MIN [mb]`))+ #colour = "red",
  guides(colour=FALSE, size=FALSE)+
  scale_fill_distiller(palette = "OrRd", direction = 1)+
  labs(x = "Longitude [deg]",
       y = "Latitude [deg]",
       title = "Minimalne vrednosti razlika na stanicama nakon krosvalidacije",
       caption = "Coordinate Reference System - WGS84")+
  geom_map(data=dat[dat$region=="Serbia",], map=dat[dat$region=="Serbia",],
           aes(x=long, y=lat, map_id=region),
           color="white", fill="#7f7f7f", size=0.05, alpha=1/4) +
  theme_bw()+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) 
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 8, fig.height = 8, fig.align='center'
ok_plot_rmse
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 8, fig.height = 8, fig.align='center'
ok_plot_max
#+ echo = FALSE, warning = FALSE, message = FALSE, cache = TRUE, fig.width = 8, fig.height = 8, fig.align='center'
ok_plot_min

#' # **Prediktori (*Covariates, Explanatory variables*)**
#' ### Digitalni model terena (*DEM - Digital Elevation Model*)
#' ### Temperatura (vremenska serija dnevnih karata temperature vazduha 2009-2018)

load("d:/Doktorske_akademske_studije_Geodezija_i_geoinformatika/Semestar_2/Prostorno_vremenska_statisitika/Temperature [Sale]/folder/ogimet_temp_stfdf.rda")



#' # **Prostorno-vremenski Regresioni kriging (*Spatio-temporal Regression Kriging*)**
#' ### Ocena trenda - linearni regresioni model
#' ### Variogram - modelovanje i kvantifikacija prostorne korelisanosti     
#' ### Empirijski prostorno-vremenski variogram       
#' ### Predikcija 
#' ### Krosvalidacija 

#' # **Prostorno-vremenska predikcija koristeći *Random Forest* kao tehniku Mašinskog učenja (*Spatio-temporal Random Forest*)**
#' ### Krosvalidacija - *5 fold crossvalidation*
#' ### Ocena trenda - *variable importance*
#' ### Predikcija

