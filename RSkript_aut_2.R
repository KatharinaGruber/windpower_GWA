library(lubridate)
library(tibble)
library(feather)
library(dplyr)
library(tidyverse)
library(Metrics)
library(reshape2)
library(ggplot2)
library(BBmisc)
library(readxl)
library(hash)
library(gtools)
library(plotly)
library(raster)
library(rgdal)
library(ncdf4)
library(httr)
library(parallel)
library(forecast)
library(tseries)
library(fitdistrplus)
library(zoo)

library(rjson)
library(stringr)
library(tidyr)
library(reticulate)

# base directory for simulation
dirbase <- "C:/Users/KatharinaG/Documents/DOK/IEWT"
# base ERA5 data directory
direrabase <- "C:/Users/KatharinaG/Documents/era5/data_aut"
# directory where era5 data per point are stored
direra_aut <- "C:/Users/KatharinaG/Documents/era5/data_aut/erapoints"
# dirwindatlas
dirwindatlas <- "C:/Users/KatharinaG/Documents/DOK/IEWT/windatlas/aut"
# directory where observed wind power generation data from ÖMAG are stored
diromag <- paste0(dirbase,"/wind_power_oemag")


source("C:/Users/KatharinaG/Documents/DOK/IEWT/ERA5_data.R")
source("C:/Users/KatharinaG/Documents/DOK/IEWT/functions_aut2.R")





########################################################################################################################
##### DOWNLOAD OF ERA5 DATA WITH PYTHON ################################################################################
########################################################################################################################

####the boundary of the box to download
####
lon1<- 8
lat1<- 45.5
lon2<- 18
lat2<- 50




# split date sequence to several parts so reduce memory requirements
date_seq<-list(seq(as.POSIXct("2000-01-01",tz="UTC"),as.POSIXct("2005-12-31",tz="UTC"),by="d"),seq(as.POSIXct("2006-01-01",tz="UTC"),as.POSIXct("2011-12-31",tz="UTC"),by="d"),seq(as.POSIXct("2012-01-01",tz="UTC"),as.POSIXct("2017-12-31",tz="UTC"),by="d"))
setwd(direrabase)

lapply(date_seq,convertERAFeather,"u10","10m_u")
lapply(date_seq,convertERAFeather,"u100","100m_u")
lapply(date_seq,convertERAFeather,"v10","10m_v")
lapply(date_seq,convertERAFeather,"v100","100m_v")


lonlat<-read_feather(paste0(paste("./feather/u10",format(date_seq[[1]][1],"%Y%m%d"),format(date_seq[[1]][length(date_seq[[1]])],"%Y%m%d"),sep="_"),"/lonlat.feather"))
names(lonlat) <- c("long","lat")
write_feather(lonlat,paste(direrabase,"/lonlat.feather",sep=""))

ERADate <- seq(date_seq[[1]][1],date_seq[[length(date_seq)]][length(date_seq[[length(date_seq)]])],by="h")
hours <- as.POSIXct(rep(ERADate[length(ERADate)],23),tz="UTC")
hours <- hours + (1:23)*3600
ERADate <- c(ERADate,hours)
write_feather(as.data.frame(ERADate),paste(direrabase,"/ERADate.feather",sep=""))




ERADate <- read_feather(paste(direrabase,"/ERADate.feather",sep=""))
LonLat <- read_feather(paste(direrabase,"/lonlat.feather",sep=""))

pnames <- c("u10","u100","v10","v100")


invisible(apply(LonLat,1,saveERAPoint,pnames,date_seq))




###############################################################################################################################
##### DOWNLOAD IGWINDKRAFT DATA ###############################################################################################
###############################################################################################################################






# import IG-Wind data
json_file <- "https://www.igwindkraft.at/src_project/external/maps/generated/gmaps_daten.js"
lines <- readLines(json_file)
lines[1] <- sub(".* = (.*)", "\\1", lines[1])
lines[length(lines)] <- sub(";", "", lines[length(lines)])
json_data <- fromJSON(paste(lines, collapse="\n"))


# extract data into data frame
col_names <- c("Name","Betreiber1","Betreiber2","n_Anlagen","KW","Type","Jahr","x","Lat","Long","url","Hersteller","Nabenhoehe","Rotordurchmesser")
wind_turbines_data        <- data.frame(matrix(ncol = 14, nrow = 0))
names(wind_turbines_data) <- col_names


for (i in seq(json_data[[2]]$places)){
  for(j in seq(col_names)){
    if (!is.null(json_data[[2]]$places[[i]]$data[[j]])){
      wind_turbines_data[i,j] <- json_data[[2]]$places[[i]]$data[[j]]
    } else {
      wind_turbines_data[i,j] <- NA 
    }
  }
}

wind_turbines_data$x <- NULL
wind_turbines_data <- wind_turbines_data %>% mutate(Name_Save=Name)
name<-wind_turbines_data$Name
name<-str_replace(name,"ä","ae")
name<-str_replace(name,"ö","oe")
name<-str_replace(name,"ü","ue")
name<-str_replace(name,"ß","sz")
name<-str_replace(name," ","_")
name<-str_replace(name,"\\(","-")
name<-str_replace(name,"\\)","-")

wind_turbines<-wind_turbines_data
wind_turbines$Name<-name

wind_turbines<-wind_turbines %>% mutate(Park=matrix(unlist(strsplit(wind_turbines$Name,",")), 
                                                    ncol=2, 
                                                    byrow=TRUE)[,1]) %>% mutate(Park=str_replace(Park,"/",""))

wind_turbines<-remove_erroneous_data(wind_turbines)

wind_turbines<-wind_turbines %>% mutate(Nabenhoehe=as.numeric(Nabenhoehe),Rotordurchmesser=as.numeric(Rotordurchmesser))

setwd(dirbase)
save(wind_turbines,file="wind_turbines_AUT.RData")



#################################################################################################################################
##### CALCULATE CAPACITY CORRECTION FACTOR ######################################################################################
#################################################################################################################################
date.start <- as.POSIXct("2003-01-01",tz="UTC")


load(paste0(dirbase,"/wind_turbines_AUT.RData"))
# prepare wind turbine data
# add missing commissioning dates
# Oberwaltersdorf is missing turbine type, add "V112"
wind_turbines$Type[grep("Oberwaltersdorf",wind_turbines$Name)] <- "V112"
# load information of rest missing data
mis_comdate <- read.csv(file=paste0(dirbase,"/missing_comdates_aut.csv"))
# add information
for(i in c(1:length(mis_comdate$location))){
  ind <- grep(mis_comdate$location[i],wind_turbines$Name)
  ind <- ind[which(wind_turbines$Jahr[ind]=="")]
  types <- wind_turbines$Type[ind]
  wind_turbines$Jahr[ind[which(types==mis_comdate$type[i])]] <- mis_comdate$year[i]
  wind_turbines$Nabenhoehe[ind[which(types==mis_comdate$type[i])]] <- mis_comdate$height[i]
  wind_turbines$Rotordurchmesser[ind[which(types==mis_comdate$type[i])]] <- mis_comdate$diameter[i]
}
# exact commissioning date is not known, therefore set to middle of year (1st of july)
wind_turbines <- data.frame(wind_turbines,comdate=as.POSIXct(paste(wind_turbines$Jahr,"-07-01 00:00:00",sep=""),tz="UTC",format="%Y-%m-%d %H:%M:%OS"))
wind_turbines <- wind_turbines[which(wind_turbines$comdate < as.POSIXct("2018-01-01 00:00:00",tz="UTC")),]
# make all dates that are before start date to start date
wind_turbines$comdate[which(wind_turbines$comdate<date.start)] <- date.start

# save completed wind turbine data
save(wind_turbines,file=paste0(dirbase,"/wind_turbines_AUT_complete.RData"))

capdate_wp <- aggregate(wind_turbines$KW,by=list(wind_turbines$comdate),sum)
capdate_wp[,2] <- cumsum(capdate_wp[,2])/1000
names(capdate_wp) <- c("time","cap_MW")

# prepare time series
caps1 <- data.frame(time=seq(as.POSIXct("2003-01-01",tz="UTC"),as.POSIXct("2017-12-31",tz="UTC"),by="d"))
caps1$cap <- rep(NA,length(caps1$time))
caps1$cap[match(capdate_wp$time,caps1$time)] <- capdate_wp$cap_MW
caps1$cap <- na.locf(caps1$cap)
# cut before 2008 because ömag data start at 2008
caps1 <- caps1[which(caps1$time>=as.POSIXct("2008-01-01",tz="UTC")),]

# read oemag capacity data
caps_omag <- read.csv2(paste0(dirbase,"/inst_cap.csv"))
caps_omag$cap_MW <- as.numeric(as.vector(caps_omag$cap_MW))
caps_omag$comdate <- as.POSIXct(paste0(caps_omag$ï..year,"-",(caps_omag$quartal-1)*3+1,"-01 00:00:00"),tz="UTC")

capdate_omag <- data.frame(time=caps_omag$comdate,cap_MW=caps_omag$cap_MW)

# prepare time series
caps2 <- data.frame(time=caps1$time,cap=rep(NA,length(caps1$time)))
caps2$cap[match(capdate_omag$time,caps2$time)] <- capdate_omag$cap_MW
caps2$cap <- na.locf(caps2$cap)

cf_caps <- sum(caps2$cap)/sum(caps1$cap)

save(cf_caps,file=paste0(dirbase,"/cf_cap.RData"))

#################################################################################################################################
##### wind speed correction with WindAtlas data #################################################################################
#################################################################################################################################

# general data
LonLat <- read_feather(paste(direrabase,"/lonlat.feather",sep=""))
# start date 2003 because startdate of ÖMAG data
date.start <- as.POSIXct("2003-01-01",tz="UTC")
rad <- pi/180
# Enercon E-82
ratedpower <- 2000
height <- 108
windspeed <- c(0:25)
powercurve <- c(0,0,3,25,82,174,312,532,815,1180,1580,1810,1980,rep(2050,13))

# prepare wind atlas data (for faster loading)
tif = raster(paste(dirwindatlas,"/windatlas.tif",sep=""))
windatlas <- rasterToPoints(tif)
save(windatlas,file=paste(dirwindatlas,"/windatlas.RData",sep=""))

# calculate wind atlas corrected wind power generation
statpower_meanAPT <- calcstatpower(ratedpower,windspeed,powercurve,useWA=1)
save(statpower_meanAPT,file=paste0(dirbase,"/results/statpowerWA.RData"))
rm(statpower_meanAPT)

# calculate wind power generation without wind power correction
statpower <- calcstatpower(ratedpower,windspeed,powercurve,useWA=0)
save(statpower,file=paste0(dirbase,"/results/statpower.RData"))
rm(statpower)



# sum over austria
load(paste0(dirbase,"/results/statpowerWA.RData"))
wpaut_meanAPT <- statpower_meanAPT[[1]]
for(i in c(2:length(statpower_meanAPT))){
  wpaut_meanAPT[,2] <- wpaut_meanAPT[,2] + statpower_meanAPT[[i]][,2]
}
names(wpaut_meanAPT) <- c("time","wp_kwh")


load(paste0(dirbase,"/results/statpower.RData"))
wpaut <- statpower[[1]]
for(i in c(2:length(statpower))){
  wpaut[,2] <- wpaut[,2] + statpower[[i]][,2]
}
names(wpaut) <- c("time","wp_kwh")

save(wpaut,wpaut_meanAPT,file=paste0(dirbase,"/results/wp_austria.RData"))



#################################################################################################################################
##### analyse results ###########################################################################################################
#################################################################################################################################

# load results
load(paste0(dirbase,"/results/wp_austria.RData"))
# correct capacities
load(paste0(dirbase,"/cf_cap.RData"))
wpaut$wp_kwh <- wpaut$wp_kwh*cf_caps/1000
names(wpaut) <- c("time","wp_MWh")
wpaut_meanAPT$wp_kwh <- wpaut_meanAPT$wp_kwh*cf_caps/1000
names(wpaut_meanAPT) <- c("time","wp_MWh")

# load oemag data
files <- list.files(diromag)

years <- c(2003:2017)
for(i in c(1:length(files))){
  if(years[i]<2012){
    file1 <- read_xls(paste0(diromag,"/",files[i]))
  }else{
    file1 <- read_xlsx(paste0(diromag,"/",files[i]))
  }
  if(years[i]==2005){
    prody <- as.numeric(unlist(file1[-c(1:4),4]))
  }else{
    prody <- as.numeric(unlist(file1[-c(1:3),4]))
  }
 
  time <- rep(seq(as.POSIXct(paste0(years[i],"-01-01 00:00:00",tz="UTC")),as.POSIXct(paste0(years[i],"-12-31 23:00:00",tz="UTC")),by="h"),each=4)
  prod_h <- aggregate(prody,by=list(time),sum)
  if(i == 1){
    prod_omag <- prod_h
  } else {
    prod_omag <- rbind(prod_omag,prod_h)
  }
  
}

prod_omag[,2] <- prod_omag[,2]/1000
names(prod_omag) <- c("time","wp_MWh")



# aggregate daily for also a daily comparison
prod_omag_daily <- aggregate(prod_omag$wp_MWh,by=list(format(prod_omag$time,"%Y%m%d")),sum)
  
wpaut_meanAPT_daily <- aggregate(wpaut_meanAPT$wp_MWh,by=list(format(wpaut_meanAPT$time,"%Y%m%d")),sum)
wpaut_daily <- aggregate(wpaut$wp_MWh,by=list(format(wpaut$time,"%Y%m%d")),sum)

comp <- data.frame(time=seq(as.POSIXct("2003-01-01",tz="UTC"),as.POSIXct("2017-12-31",tz="UTC"),by="d"),
                   obs=prod_omag_daily[,2],
                   wp=wpaut_daily[,2],
                   wpc=wpaut_meanAPT_daily[,2])

comp[,2:4] <- comp[,2:4]/1000


# calculate statistics
stats <- data.frame(type=c("RMSE","MBE","mean"),wp=NA,wpc=NA,obs=NA)
# RMSE
stats[1,2:3] <- c(round(rmse(comp$obs,comp$wp),2),
                  round(rmse(comp$obs,comp$wpc),2))
#MBE
stats[2,2:3] <- c(round(sum(comp$wp-comp$obs),2),
                  round(sum(comp$wpc-comp$obs),2))
# mean
stats[3,2:4] <- c(round(mean(comp$wp),2),
                  round(mean(comp$wpc),2),
                  round(mean(comp$obs),2))

write.table(stats,file=paste0(dirbase,"/stats.csv"),sep=";")


comp_t <- gather(comp,"wp","wpc","obs",key="type",value="wp_GWh")
#comp_t$wp_GWh <- comp_t$wp_GWh/1000

# boxplots
ggplot(data=comp_t,aes(x=type,y=wp_GWh)) +
  geom_boxplot() +
  xlab("data (observed/simulated corrected/uncorrected") +
  scale_y_continuous(name="daily wind power generation [GWh]") +
  labs(title="Comparison of simulated and observed daily wind power generation")
ggsave(paste(dirbase,"/daily_comp.png",sep=""), width = 8, height = 5.5)



