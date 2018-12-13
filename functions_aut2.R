


# function for download of wind park data from ig windkraft
###remove / correct wrong data points in original file
remove_erroneous_data<-function(wind_turbines){
  
  ##zurndorf: Lat/Long vertauscht
  ##zurndorf: Lat/Long vertauscht
  selector<-wind_turbines$Lat<40
  helper<-wind_turbines[selector,]$Lat
  wind_turbines[selector,]$Lat<-wind_turbines[selector,]$Long
  wind_turbines[selector,]$Long<-helper
  
  ##nickelsdorf anlage 9: falsch
  selector<-wind_turbines$Long>17.39
  wind_turbines[selector,]$Long<-17+(wind_turbines[selector,]$Long-17)/10
  
  ##anlagen munderfing 5 und 6 eigenartig
  wind_turbines[wind_turbines$Park=="Munderfing"&wind_turbines$Long>14,]$Long<-
    wind_turbines[wind_turbines$Park=="Munderfing"&wind_turbines$Long>14,]$Long-3.5
  
  ##Correct very large latitude values (e.g. Pottendorf)
  wind_turbines[wind_turbines$Lat>1000,]$Lat<-wind_turbines[wind_turbines$Lat>1000,]$Lat/10^5
  
  ##Gro?hofen: one turbine far away from others
  wind_turbines[wind_turbines$Park=="Groszhofen"&wind_turbines$Long<16.6,]$Long<-
    wind_turbines[wind_turbines$Park=="Groszhofen"&wind_turbines$Long<16.6,]$Long+0.52
  
  return(wind_turbines)
} 




calcstatpower <- function(ratedpower,windspeed,powercurve,useWA){
  # get data of windparks: capacities and start dates and sort by start dates for each location
  load(paste0(dirbase,"/wind_turbines_AUT_complete.RData"))
  # exact commissioning date is not known, therefore set to middle of year (1st of july)
  wind_turbines <- data.frame(wind_turbines,comdate=as.POSIXct(paste(wind_turbines$Jahr,"-07-01 00:00:00",sep=""),tz="UTC",format="%Y-%m-%d %H:%M:%OS"))
  wind_turbines <- wind_turbines[which(wind_turbines$comdate < as.POSIXct("2018-01-01 00:00:00",tz="UTC")),]
  # remove NAcomdates
  wind_turbines <- wind_turbines[which(!is.na(wind_turbines$comdate)),]
  wind_turbines$comdate[which(wind_turbines$comdate<date.start)] <- date.start
  statpowlist <- list()
  # list for saving correction factors
  cfs_mean <<-list()
  for(i in c(1:length(wind_turbines$Long))){
    print(i)
    # find nearest neightbour ERA5 and extrapolate to hubheight
    long <<- wind_turbines$Long[i]
    lat <<- wind_turbines$Lat[i]
    lldo <<- distanceorder()
    # some hub heights are missing, then set to 108m
    height <- wind_turbines$Nabenhoehe[i]
    if(is.na(height)){
      height <- 108
    }
    NNera <- NNdf(height)
    
    # if wind atlas wind speed correction shall be applied (useWA = 1)
    if(useWA > 0){
      # WIND ATLAS wind speed correction
      load(paste(dirwindatlas,"/windatlas.RData",sep=""))
      ppWAdistance <- 6378.388*acos(sin(rad*wind_turbines$Lat[i]) * sin(rad*windatlas[,2]) + cos(rad*wind_turbines$Lat[i]) * cos(rad*windatlas[,2]) * cos(rad*windatlas[,1]-rad*wind_turbines$Long[i]))
      
      # find data of nearest station
      pointn <- which(ppWAdistance==min(ppWAdistance))
      long <<- windatlas[pointn,1]
      lat <<- windatlas[pointn,2]
      # get wind speed data of nearest ERA5 point (at 50m height! as wind atlas data)
      lldo <<- distanceorder()
      windERA100m <- NNdf(100)
      cf <- as.numeric(windatlas[pointn,3])/mean(windERA100m[,2])
      cfs_mean[[i]] <<- cf
      # adapt mean wind speed
      NNera[,2] <- NNera[,2]*cf
    }
    
    
    # make a list of capacities for all dates
    caplist <- data.frame(NNera[,1],rep(0,length(NNera[,1])))
    match <- match(wind_turbines$comdate[i],caplist[,1])
    caplist[match:length(caplist[,1]),2] <- wind_turbines$KW[i]
    
    # calculate power output for all hours from power curve in kWh
    # values are interpolated linearly betweer points of power curve
    whichs <- findInterval(NNera[,2],windspeed)
    whichs.1 <- whichs+1
    whichs.1[which(whichs.1>length(powercurve))] <- length(powercurve)
    statpower <- caplist[,2]/ratedpower*((powercurve[whichs]-powercurve[whichs.1])/(windspeed[whichs]-windspeed[whichs.1])*(NNera[,2]-windspeed[whichs.1])+powercurve[whichs.1])
    # replace NAs created where which is last element of powercurve with 0 because above cut out wind speed
    statpower[which(whichs==length(powercurve))] <- 0
    
    statpowlist[[i]] <- data.frame(NNera[,1],statpower)
    
  }
  
  return(statpowlist)
}



#function for calculating distances and order by distances for ERA5
distanceorder <- function(){
  distance <- 6378.388*acos(sin(rad*lat) * sin(rad*LonLat$lat) + cos(rad*lat) * cos(rad*LonLat$lat) * cos(rad*LonLat$long-rad*long))
  lonlatdistao <- data.frame(LonLat,distance,c(1:length(distance)))
  names(lonlatdistao) <- c("Longitude","Latitude","distance","ERAnum")
  lonlatdistao <- lonlatdistao[order(lonlatdistao[,3]),]
  return(lonlatdistao)
}



NNdf <- function(hubheight=10){
  setwd(direra_aut)
  ########## 1. Nearest Neighbour ##########
  # first row (nearest neighbour) is extracted
  # first column of list is taken (distance to station)
  # and columns 4 to 27, long and lat are excluded
  ERAdf <- getERAPoint(lldo$Longitude[1],lldo$Latitude[1])
  
  # Wind speeds Nearest Neighbor
  WHuv100 <- sqrt(ERAdf$u100^2+ERAdf$v100^2)
  WHuv10 <- sqrt(ERAdf$u10^2+ERAdf$v10^2)
  
  # WH ERA5-DF
  MWH <- data.frame(ERAdf$ERADate,WHuv100,WHuv10)
  MWH1 <- MWH[which(MWH[,1]>=date.start),]
  
  MWH1ext <-data.frame(MWH1[,1],extrap(MWH1,hubheight))
  names(MWH1ext) <- c("date","vext")
  
  return(MWH1ext)
}


#Extrapolation Höhe
#Höhe in der INMET Daten gemessen werden: 10m
extrap <- function(MWH1,hIN=10){
  #alpha friction coefficient
  alpha <- (log(MWH1[,2])-log(MWH1[,3]))/(log(100)-log(10))
  #wind speed with power law
  vext <- MWH1[,2]*(hIN/100)^alpha
  return(vext)
}
