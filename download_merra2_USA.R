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

source("C:/Users/KatharinaG/andereDokumente/Masterarbeit2/daten/MERRA/skripts/MERRA_data.R")



lon1 <- -173
lat1 <- 12
lon2 <- -64
lat2 <- 68


date_seq<-seq(as.POSIXct("2000-12-01",tz="UTC"),as.POSIXct("2018-12-31",tz="UTC"),by="d")

setwd("C:/Users/KatharinaG/Data/MERRA/USA")

getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("U10M"),
                "RE_EXTREME",
                "Re_extreme666!",
                TRUE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("U50M"),
                "RE_EXTREME",
                "Re_extreme666!",
                TRUE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("V10M"),
                "RE_EXTREME",
                "Re_extreme666!",
                TRUE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("V50M"),
                "RE_EXTREME",
                "Re_extreme666!",
                TRUE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("DISPH"),
                "RE_EXTREME",
                "Re_extreme666!",
                TRUE)
