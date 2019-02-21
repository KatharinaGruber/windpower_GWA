# make sum of wind power generation for all of brazil
sum_brasil <- function(complist){
  df <- complist[[1]]
  for(i in c(2:length(complist))){
    df[,2] <- df[,2]+complist[[i]][,2]
  }
  return(df)
}




# function for loading measured wind power for subsystems or brazil
getprodSUBBRA <- function(area){
  a <- data.frame(a=c("NE","S","BRASIL"),b=c("nordeste","sul","brasil"))
  prod <- read.table(paste(dirwindprodsubbra,"/",a$b[match(area,a$a)],"_dia.csv",sep=""),sep=";",header=T,stringsAsFactors=F)
  # extract date and generation in GWh
  prod <- data.frame(date=as.numeric(as.vector(paste(substr(prod[,1],7,10),substr(prod[,1],4,5),substr(prod[,1],1,2),sep=""))),prod_GWh=as.numeric(gsub(",",".",prod[,8],fixed=T)))
  # cut last line because useless
  prod <- prod[2:length(prod[,1]),]
  return(prod)
}


# function that cuts two data frames to same length
# data frames with two columns, first column has dates, second has data
csl <- function(df1,df2){
  cut1 <- max(df1[1,1],df2[1,1])
  cut2 <- min(df1[nrow(df1),1],df2[nrow(df2),1])
  df1.1 <- df1[which(df1[,1]==cut1):which(df1[,1]==cut2),]
  df2.1 <- df2[which(df2[,1]==cut1):which(df2[,1]==cut2),]
  df <- data.frame(df1.1[,1],df1.1[,2],df2.1[,2])
  return(df)
}





# calculate wind power generation per location (wind turbine, wind park or several wind parks at same locations)
# ratedpower and height refer to wind turbines, windspeed and powercurve contain the power curve data
# if useWA is 1, the wind atlas bias correction is performed
calcstatpower <- function(ratedpower,height,windspeed,powercurve,useWA){
  # get data of windparks: capacities and start dates and sort by start dates for each location
  load(paste(dirwindparks,"/windparks_complete.RData",sep=""))
  windparks <- windparks[order(windparks$long),]
  windparks <- data.frame(windparks,comdate=as.POSIXct(paste(windparks$year,"-",windparks$month,"-",windparks$day," 00:00:00",sep=""),tz="UTC"))
  windparks <- windparks[which(windparks$comdate < as.POSIXct("2017-08-31 00:00:00",tz="UTC")),]
  windparks$comdate[which(windparks$comdate<date.start)] <- date.start
  numlon <- as.vector(unlist(rle(windparks$long)[1]))
  # counter for locations
  pp <- 1
  statpowlist <- list()
  powlistind <- 1
  # list for saving correction factors
  cfs_mean <<-list()
  while(pp<=length(windparks$long)){
    print(powlistind)
    numstat <- numlon[powlistind]
    pplon <- windparks$long[pp]
    pplat <- windparks$lat[pp]
    # find nearest neightbour ERA5 and extrapolate to hubheight
    long <<- pplon
    lat <<- pplat
    lldo <<- distanceorder()
    NNera <- NNdf(height)
    if(useWA>0){
      # WIND ATLAS wind speed correction
      load(paste(dirwindatlas,"/windatlas.RData",sep=""))
      ppWAdistance <- 6378.388*acos(sin(rad*pplat) * sin(rad*windatlas[,2]) + cos(rad*pplat) * cos(rad*windatlas[,2]) * cos(rad*windatlas[,1]-rad*pplon))
      
      # find data of nearest station
      pointn <- which(ppWAdistance==min(ppWAdistance))
      long <<- windatlas[pointn,1]
      lat <<- windatlas[pointn,2]
      # get wind speed data of nearest ERA5 point (at 100m height! as wind atlas data)
      lldo <<- distanceorder()
      windERA100m <- NNdf(100)
      cf <- as.numeric(windatlas[pointn,3])/mean(windERA100m[,2])
      cfs_mean[[powlistind]] <<- cf
      # adapt mean wind speed
      NNera[,2] <- NNera[,2]*cf
    }
    
    # get startdates and capacities from municipios
    capdate <- data.frame(windparks$comdate[pp:(pp+numstat-1)],windparks$cap[pp:(pp+numstat-1)],rep(NA,numstat))
    names(capdate) <- c("commissioning","capacity","capacitysum")
    capdate <- capdate[order(capdate$commissioning),]
    capdate$capacitysum <- cumsum(capdate$capacity)
    # make a list of capacities for all dates
    caplist <- data.frame(NNera[,1],rep(0,length(NNera[,1])))
    match <- match(capdate$commissioning,caplist[,1])
    for(i in c(1:length(match))){
      caplist[match[i]:length(caplist[,1]),2] <- capdate$capacitysum[i]
    }
    
    
    # calculate power output for all hours from power curve in kWh
    # values are interpolated linearly betweer points of power curve
    whichs <- findInterval(NNera[,2],windspeed)
    whichs.1 <- whichs+1
    whichs.1[which(whichs.1>length(powercurve))] <- length(powercurve)
    statpower <- caplist[,2]/ratedpower*((powercurve[whichs]-powercurve[whichs.1])/(windspeed[whichs]-windspeed[whichs.1])*(NNera[,2]-windspeed[whichs.1])+powercurve[whichs.1])
    # replace NAs created where which is last element of powercurve with 0 because above cut out wind speed
    # because powercurve becomes flat and no higher power is generated
    statpower[which(whichs==length(powercurve))] <- 0
    
    statpowlist[[powlistind]] <- data.frame(NNera[,1],statpower)
    
    pp <- pp + numstat
    powlistind <- powlistind +1
    
  }
  
  return(statpowlist)
}

# nearest neighbour interpolation and extrapolation to hubheight
NNdf <- function(hubheight){
  setwd(direra_bra)
  # first row (nearest neighbour) is extracted
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


# Extrapolation Height
extrap <- function(MWH1,height){
  # alpha friction coefficient
  alpha <- (log(MWH1[,2])-log(MWH1[,3]))/(log(100)-log(10))
  # wind speed with power law
  vext <- MWH1[,2]*(height/100)^alpha
  return(vext)
}


#function for calculating distances and order by distances for ERA5
distanceorder <- function(){
  distance <- 6378.388*acos(sin(rad*lat) * sin(rad*LonLat$lat) + cos(rad*lat) * cos(rad*LonLat$lat) * cos(rad*LonLat$long-rad*long))
  lonlatdistao <- data.frame(LonLat,distance,c(1:length(distance)))
  names(lonlatdistao) <- c("Longitude","Latitude","distance","ERAnum")
  lonlatdistao <- lonlatdistao[order(lonlatdistao[,3]),]
  return(lonlatdistao)
}