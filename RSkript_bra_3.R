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
dirbase <- "C:/..."
# base ERA5 data directory
direrabase <- "C:/..."
# directory where era5 data per point are stored
direra_bra <- "C:/..."
# dirwindatlas
dirwindatlas <- "C:/..."
# directory where observed wind power generation data from ONS are stored
dirwindprod <- "C:/..."
# directory where capacities provided by ONS are stored
dircaps <- "C:/..."
# directory where windparkdata from thewindpower.net are stored
dirwindparks <- "C:/..."
# directory where windparkdata from thewindpower.net for selected wind parks are stored
dirwindparks_sel <- "C:/..."


source("C:/.../ERA5_data_bra.R")
source("C:/.../functions_bra2.R")





########################################################################################################################
##### DOWNLOAD OF ERA5 DATA WITH PYTHON ################################################################################
########################################################################################################################


repl_python()

import cdsapi
import numpy as np
import calendar
import os

# define directory in which data shall be stored
os.chdir("C:/Users/KatharinaG/Documents/DOK/era5/data_aut")

c = cdsapi.Client()



# define the years you want to download
yearstart = 2000
yearend = 2017
# define the start and end month you want to download
monthstart = 1
monthend = 12
# define the start and end day you want to download
daystart = 1
dayend = 31

# define spatial limits of download (eg. around Austria)
lon1 = -74.1
lon2 = -36
lat1 = -33
lat2 = 5.5

# create lists
years = np.array(range(yearstart,yearend+1),dtype="str")
lonlat = [lat2, lon1, lat1, lon2]


# define variable to download (eg. u component 100m wind)
vars = np.array(['100m_u_component_of_wind','10m_u_component_of_wind','100m_v_component_of_wind','10m_v_component_of_wind'])

for var in vars:
  for year in years:
  if (int(year)==yearstart) and (int(year)==yearend):
  months = np.array(range(monthstart,monthend+1),dtype="str")
elif (year == yearstart) :
  months = np.array(range(monthstart,13),dtype="str")
elif (year == yearend):
  months = np.array(range(1,monthend + 1),dtype="str")
else:
  months = np.array(range(1,13),dtype="str")

for month in months:
  
  if int(month) < 10:
  m = '0' + month
else:
  m = month

if(int(year) == yearstart) and (int(year) == yearend) and (int(month) == monthstart) and (int(month) == monthend):
  days = list(np.array(range(daystart,dayend+1),dtype="str"))
elif (int(year) == yearstart) and (int(month) == monthstart):
  days = list(np.array(range(daystart,calendar.monthrange(int(year),int(month))[1]+1),dtype="str"))
elif (int(year) == yearend) and (int(month) == monthend):
  days = list(np.array(range(1,dayend+1),dtype="str"))
else:
  days = list(np.array(range(1,calendar.monthrange(int(year),int(month))[1]+1),dtype="str"))

for day in days:
  if int(day) < 10:
  d = '0' + day
else:
  d = day
if not (('era5_' + var + '_' + year + m + d + '.nc') in os.listdir()):                
  c.retrieve(
    'reanalysis-era5-single-levels',
    {
      'variable': var,
      'product_type':'reanalysis',
      'year': year,
      'month': month,
      'day': day,
      'time':[
        '00:00','01:00','02:00',
        '03:00','04:00','05:00',
        '06:00','07:00','08:00',
        '09:00','10:00','11:00',
        '12:00','13:00','14:00',
        '15:00','16:00','17:00',
        '18:00','19:00','20:00',
        '21:00','22:00','23:00'
        ],
      'area': lonlat,
      'format':'netcdf',
    },
    'era5_' + var + '_' + year + m + d + '.nc')

exit




# split date sequence to several parts so reduce memory requirements
# as data for Brazil are relatively big, it is necessary to splice the time steps into months for loading, but they are joined to years within
# the function where they are transformed from file per timestep to file per point
# therefore each year is handled seperately automatically by a loop
startyear <- 2006
endyear <- 2017
date_seq <- list()
for(year in c(startyear:endyear)){
  ds <- list()
  for(month in c(1:12)){
    # determine number of days in month
    daymon <- days_in_month(as.POSIXct(paste0(year,"-",month,"-01"),tz="UTC"))
    ds[[month]] <- seq(as.POSIXct(paste0(year,"-",month,"-01"),tz="UTC"),
                       as.POSIXct(paste0(year,"-",month,"-",daymon),tz="UTC"), by="d")
  }
  date_seq[[year-startyear+1]] <- ds
}

setwd(direrabase)
# join date sequences  
lapply(date_seq,convertERAFeather,"u10","10m_u")
lapply(date_seq,convertERAFeather,"v10","10m_v")
lapply(date_seq,convertERAFeather,"u100","100m_u")
lapply(date_seq,convertERAFeather,"v100","100m_v")




# create new yearly date_seq, as data were stored
date_seq2 <- list()
for(i in c(1:length(date_seq))){
  date_seq2[[i]] <- do.call("c",date_seq[[i]])
}

lonlat<-read_feather(paste0(paste("./feather/u10",format(date_seq2[[1]][1],"%Y%m%d"),format(date_seq2[[1]][length(date_seq2[[1]])],"%Y%m%d"),sep="_"),"/lonlat.feather"))
names(lonlat) <- c("long","lat")
write_feather(lonlat,paste(direrabase,"/lonlat.feather",sep=""))

ERADate <- seq(as.POSIXct(paste0(startyear,"-01-01 00:00:00"),tz="UTC"),as.POSIXct(paste0(endyear,"-12-31 23:00:00"),tz="UTC"),by="h")
write_feather(as.data.frame(ERADate),paste(direrabase,"/ERADate.feather",sep=""))




ERADate <- read_feather(paste(direrabase,"/ERADate.feather",sep=""))
LonLat <- read_feather(paste(direrabase,"/lonlat.feather",sep=""))

pnames <- c("u10","u100","v10","v100")


# join different variables to common file
# Points need to be handled in batches, as otherwise memory requirements are too high
for(i in c(1:10)){
  r <- c(((i-1)*1000+1):(i*1000))
  LL <- LonLat2[r,]
  invisible(apply(LL,1,saveERAPoint,pnames,date_seq2))
}
r <- c(10000:length(LonLat2[,1]))
LL <- LonLat2[r,]
invisible(apply(LL,1,saveERAPoint,pnames,date_seq2))




#################################################################################################################################
##### CALCULATE CAPACITY CORRECTION FACTOR ######################################################################################
#################################################################################################################################
date.start <- as.POSIXct("2006-01-01",tz="UTC")
date.end <- as.POSIXct("2017-12-31",tz="UTC")


# read wind park data from The Wind Power
load(paste(dirwindparks,"/windparks_complete.RData",sep=""))
# extract only commissiongs and capacities of NE and S
NEwp <- windparks[windparks$state %in% c("Bahia","Ceará","Paraíba","Pernambuco","Piaui","Rio Grande do Norte","Sergipe"),c(3,9,10,11)]
Swp <- windparks[windparks$state %in% c("Paraná","Santa Catarina","Rio Grande do Sul"),c(3,9,10,11)]
# create datetime commissioning dates from year month and day
Bwp_cd <- data.frame(cap=windparks$cap,comdate=as.POSIXct(paste(windparks$year,"-",windparks$month,"-",windparks$day," 00:00:00",sep=""),tz="UTC"))
NEwp_cd <- data.frame(cap=NEwp$cap,comdate=as.POSIXct(paste(NEwp$year,"-",NEwp$month,"-",NEwp$day," 00:00:00",sep=""),tz="UTC"))
Swp_cd <- data.frame(cap=Swp$cap,comdate=as.POSIXct(paste(Swp$year,"-",Swp$month,"-",Swp$day," 00:00:00",sep=""),tz="UTC"))
# agregate by date to avoid multiple same date stamps
Bwp_ag <- aggregate(Bwp_cd$cap,by=list(Bwp_cd$comdate),sum)
NEwp_ag <- aggregate(NEwp_cd$cap,by=list(NEwp_cd$comdate),sum)
Swp_ag <- aggregate(Swp_cd$cap,by=list(Swp_cd$comdate),sum)
names(Bwp_ag) <- c("comdate","cap")
names(NEwp_ag) <- c("comdate","cap")
names(Swp_ag) <- c("comdate","cap")
# calculate cumulative installed capacities and divide by 1000 to get from kW to MW
Bcap_WP <- data.frame(commissioning=Bwp_ag$comdate,cap=cumsum(Bwp_ag$cap)/1000)
NEcap_WP <- data.frame(commissioning=NEwp_ag$comdate,cap=cumsum(NEwp_ag$cap)/1000)
Scap_WP <- data.frame(commissioning=Swp_ag$comdate,cap=cumsum(Swp_ag$cap)/1000)
# cut at 2006 and 2017
bef06_B <- Bcap_WP[which(Bcap_WP$commissioning<date.start),]
bef06_NE <- NEcap_WP[which(NEcap_WP$commissioning<date.start),]
bef06_S <- Scap_WP[which(Scap_WP$commissioning<date.start),]
Bcap_WP <- Bcap_WP[which((Bcap_WP$commissioning>=date.start)&(Bcap_WP$commissioning<=date.end)),]
NEcap_WP <- NEcap_WP[which((NEcap_WP$commissioning>=date.start)&(NEcap_WP$commissioning<=date.end)),]
Scap_WP <- Scap_WP[which((Scap_WP$commissioning>=date.start)&(Scap_WP$commissioning<=date.end)),]
# add capacity on 1.1.2006
Bcap_WP <- rbind(data.frame(commissioning=date.start,cap=tail(bef06_B$cap,1)),Bcap_WP)
NEcap_WP <- rbind(data.frame(commissioning=date.start,cap=tail(bef06_NE$cap,1)),NEcap_WP)
Scap_WP <- rbind(data.frame(commissioning=date.start,cap=tail(bef06_S$cap,1)),Scap_WP)

# read capacities from ONS
Bcap_ONS <- read.table(paste(dircaps,"/Brasil.csv",sep=""),sep=";",header=T,stringsAsFactors=F)
NEcap_ONS <- read.table(paste(dircaps,"/Nordeste.csv",sep=""),sep=";",header=T,stringsAsFactors=F)
Scap_ONS <- read.table(paste(dircaps,"/Sul.csv",sep=""),sep=";",header=T,stringsAsFactors=F)
# extract yearmonth and generation in MW
Bcap_ONS <- data.frame(comdate=as.POSIXct(as.vector(paste(substr(Bcap_ONS[,2],7,10),"-",substr(Bcap_ONS[,2],4,5),"-01",sep="")),tz="UTC"),cap_MW=as.numeric(gsub(",",".",Bcap_ONS[,7],fixed=T)))
NEcap_ONS <- data.frame(comdate=as.POSIXct(as.vector(paste(substr(NEcap_ONS[,2],7,10),"-",substr(NEcap_ONS[,2],4,5),"-01",sep="")),tz="UTC"),cap_MW=as.numeric(gsub(",",".",NEcap_ONS[,7],fixed=T)))
Scap_ONS <- data.frame(comdate=as.POSIXct(as.vector(paste(substr(Scap_ONS[,2],7,10),"-",substr(Scap_ONS[,2],4,5),"-01",sep="")),tz="UTC"),cap_MW=as.numeric(gsub(",",".",Scap_ONS[,7],fixed=T)))
# reverse to have earliest at top and most recent at bottom
Bcap_ONS <- Bcap_ONS[c(length(Bcap_ONS[,1]):1),]
NEcap_ONS <- NEcap_ONS[c(length(NEcap_ONS[,1]):1),]
Scap_ONS <- Scap_ONS[c(length(Scap_ONS[,1]):1),]
# add starting capacity
Bcap_ONS <- rbind(data.frame(comdate=date.start,cap_MW=0),Bcap_ONS)
NEcap_ONS <- rbind(data.frame(comdate=date.start,cap_MW=0),NEcap_ONS)
Scap_ONS <- rbind(data.frame(comdate=date.start,cap_MW=0),Scap_ONS)

# cut at 2017 because end of ERA5 data
Bcap_ONS <- Bcap_ONS[which(Bcap_ONS$comdate<=date.end),]
NEcap_ONS <- NEcap_ONS[which(NEcap_ONS$comdate<=date.end),]
Scap_ONS <- Scap_ONS[which(Scap_ONS$comdate<=date.end),]


# create data frame for all capacities to compare them
caps_df <- data.frame(time=seq(date.start,date.end,by="day"),Bons=NA,Bwp=NA,NEons=NA,NEwp=NA,Sons=NA,Swp=NA)
caps_df[match(Bcap_ONS[,1],caps_df[,1]),2] <- Bcap_ONS[,2]
caps_df[match(Bcap_WP[,1],caps_df[,1]),3] <- Bcap_WP[,2]
caps_df[match(NEcap_ONS[,1],caps_df[,1]),4] <- NEcap_ONS[,2]
caps_df[match(NEcap_WP[,1],caps_df[,1]),5] <- NEcap_WP[,2]
caps_df[match(Scap_ONS[,1],caps_df[,1]),6] <- Scap_ONS[,2]
caps_df[match(Scap_WP[,1],caps_df[,1]),7] <- Scap_WP[,2]
# fill NAs
caps_df[,2] <- na.locf(caps_df[,2])
caps_df[,3] <- na.locf(caps_df[,3])
caps_df[,4] <- na.locf(caps_df[,4])
caps_df[,5] <- na.locf(caps_df[,5])
caps_df[,6] <- na.locf(caps_df[,6])
caps_df[,7] <- na.locf(caps_df[,7])
# proportions
cfB <- sum(caps_df$Bons)/sum(caps_df$Bwp)
cfNE <- sum(caps_df$NEons)/sum(caps_df$NEwp)
cfS <- sum(caps_df$Sons)/sum(caps_df$Swp)

save(cfB,cfNE,cfS,file=paste(dircaps,"/cap_cfs.RData",sep=""))


#################################################################################################################################
##### wind speed correction with WindAtlas data #################################################################################
#################################################################################################################################

# general data
LonLat <- read_feather(paste(direrabase,"/lonlat.feather",sep=""))
# start date 2006 because startdate of ONS data
date.start <- as.POSIXct("2006-01-01",tz="UTC")
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
statpower_meanAPT <- calcstatpower(ratedpower,height,windspeed,powercurve,useWA=1)
save(statpower_meanAPT,file=paste0(dirbase,"/results/statpowerWA_bra.RData"))
rm(statpower_meanAPT)

# calculate wind power generation without wind power correction
statpower <- calcstatpower(ratedpower,height,windspeed,powercurve,useWA=0)
save(statpower,file=paste0(dirbase,"/results/statpower_bra.RData"))
rm(statpower)



# sum over brazil
load(paste0(dirbase,"/results/statpowerWA_bra.RData"))
wpbra_meanAPT <- statpower_meanAPT[[1]]
for(i in c(2:length(statpower_meanAPT))){
  wpbra_meanAPT[,2] <- wpbra_meanAPT[,2] + statpower_meanAPT[[i]][,2]
}
names(wpbra_meanAPT) <- c("time","wp_kwh")


load(paste0(dirbase,"/results/statpower_bra.RData"))
wpbra <- statpower[[1]]
for(i in c(2:length(statpower))){
  wpbra[,2] <- wpbra[,2] + statpower[[i]][,2]
}
names(wpbra) <- c("time","wp_kwh")

save(wpbra,wpbra_meanAPT,file=paste0(dirbase,"/results/wp_brasil.RData"))




#################################################################################################################################
##### analyse results ###########################################################################################################
#################################################################################################################################

# load results
load(paste0(dirbase,"/results/wp_brasil.RData"))
# correct capacities
load(paste0(dircaps,"/cap_cfs.RData"))
wpbra$wp_kwh <- wpbra$wp_kwh*cfB/10^6
names(wpbra) <- c("time","wp_GWh")
wpbra_meanAPT$wp_kwh <- wpbra_meanAPT$wp_kwh*cfB/10^6
names(wpbra_meanAPT) <- c("time","wp_GWh")

# load ONS data
prod <- read.csv2(paste0(dirwindprod,"/brasil_dia.csv"))



prod_ons <- prod[-1,c(1,8)]
names(prod_ons) <- c("time","wp_GWh")
prod_ons$wp_GWh <- as.numeric(paste(prod_ons$wp_GWh))

# aggregate simulated time series daily
wpbra_meanAPT_daily <- aggregate(wpbra_meanAPT$wp_GWh,by=list(format(wpbra_meanAPT$time,"%Y%m%d")),sum)
wpbra_daily <- aggregate(wpbra$wp_GWh,by=list(format(wpbra$time,"%Y%m%d")),sum)
names(wpbra_meanAPT_daily) <- c("time","wp")
names(wpbra_daily) <- c("time","wp")

# bring dates to datetime format
prod_ons$time <- as.POSIXct(paste0(substr(prod_ons$time,7,10),"-",substr(prod_ons$time,4,5),"-",substr(prod_ons$time,1,2)),tz="UTC")
wpbra_meanAPT_daily$time <- as.POSIXct(paste0(substr(wpbra_meanAPT_daily$time,1,4),"-",substr(wpbra_meanAPT_daily$time,5,6),"-",substr(wpbra_meanAPT_daily$time,7,8)),tz="UTC")
wpbra_daily$time <- as.POSIXct(paste0(substr(wpbra_daily$time,1,4),"-",substr(wpbra_daily$time,5,6),"-",substr(wpbra_daily$time,7,8)),tz="UTC")


start <- max(prod_ons$time[1],wpbra_meanAPT$time[1])
end <- min(tail(prod_ons$time,1),tail(wpbra_meanAPT$time,1))

comp <- data.frame(time=seq(start,end,by="d"),
                   obs=prod_ons[which((prod_ons$time>=start)&(prod_ons$time<=end)),2],
                   wp=wpbra_daily[which((wpbra_daily$time>=start)&(wpbra_daily$time<=end)),2],
                   wpc=wpbra_meanAPT_daily[which((wpbra_meanAPT_daily$time>=start)&(wpbra_meanAPT_daily$time<=end)),2])


# calculate statistics
stats <- data.frame(type=c("RMSE","MBE","mean"),wp=NA,wpc=NA,obs=NA)
# RMSE
stats[1,2:3] <- c(round(rmse(comp$obs,comp$wp),2),
                  round(rmse(comp$obs,comp$wpc),2))
#MBE
stats[2,2:3] <- c(round(sum(comp$wp-comp$obs)/length(comp$wp),2),
                  round(sum(comp$wpc-comp$obs)/length(comp$wp),2))
# mean
stats[3,2:4] <- c(round(mean(comp$wp),2),
                  round(mean(comp$wpc),2),
                  round(mean(comp$obs),2))

write.table(stats,file=paste0(dirbase,"/stats_bra.csv"),sep=";")


comp_t <- gather(comp,"wp","wpc","obs",key="type",value="wp_GWh")

# boxplots
ggplot(data=comp_t,aes(x=type,y=wp_GWh)) +
  geom_boxplot() +
  xlab("data (observed/simulated corrected/uncorrected)") +
  scale_y_continuous(name="daily wind power generation [GWh]")
ggsave(paste(dirbase,"/daily_comp_bra.png",sep=""), width = 8, height = 5.5)
ggsave(paste(dirbase,"/daily_comp_bra_pres.png",sep=""), width = 4, height = 5.3)





