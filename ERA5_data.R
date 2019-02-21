###########################################################################################
######################################FUNCTIONS############################################
###########################################################################################


####reads era5 file from disk
####creates file of lons lats
####saves all points into single feather files
####joining the respective time-period together
convertERAFeather<- function(date_seq,pname,pname1) {

  ncname_head<-paste(direrabase,"/era5_",pname1,"_component_of_wind_",sep="")
  
    
  dates<-format(date_seq,"%Y%m%d")
  ncfiles<-paste(ncname_head,dates,".nc",sep="")
  ncfile <- nc_open(ncfiles[1])
  
  
  #Longitude
  longitude <- ncvar_get(ncfile, "longitude", verbose = F)
  nlon <- dim(longitude)
  
  #Latitude
  latitude <- ncvar_get(ncfile, "latitude", verbose = F)
  nlat <- dim(latitude)
  
  lonlat <- expand.grid(longitude,latitude)
  nc_close(ncfile)
  
  res<-sapply(ncfiles,readSingleParam,pname,nlon,nlat,simplify=FALSE)
  df<-bind_rows(res)
  
  ###write era5 points for whole period
  dir<-paste0(direrabase,"/feather")
  dir.create(dir, showWarnings = FALSE)
  dir<-paste0(direrabase,"/feather/",paste(pname,
             format(date_seq[1],"%Y%m%d"),
             format(date_seq[length(date_seq)],"%Y%m%d"),sep="_"))
  dir.create(dir, showWarnings = FALSE)

  ###write files
  for(i in 1:ncol(df)){
    outfile<-paste(dir,"/",lonlat[i,1],"_",lonlat[i,2],".feather",sep="")
    write_feather(df[,i],outfile)
  }

  write_feather(lonlat,paste(dir,"/lonlat.feather",sep=""))
}

######read a single parameter from a single ERA5-ncfile
readSingleParam<-function(ncfileIn,pname,nlon,nlat){
  print(ncfileIn)
  ncfile<-nc_open(ncfileIn)
  #Time
  time <- ncvar_get(ncfile, "time")
  tunits <- ncatt_get(ncfile, "time", "units")
  ntime <- dim(time)
  m.array <- ncvar_get(ncfile, pname)
  m.vec.long <- as.vector(m.array)
  m.mat <- t(matrix(m.vec.long, nrow = nlon * nlat, ncol = ntime))
  nc_close(ncfile)
  return(as_tibble(m.mat))
}


# merge era variables from different files into one file
mergeERAvars <- function(date_seq,pnames,long,lat){
  dir<-paste0(direrabase,"/feather/",paste(pnames,
             format(date_seq[1],"%Y%m%d"),
             format(date_seq[length(date_seq)],"%Y%m%d"),sep="_"))
  
  files <- paste0(dir,"/",long,"_",lat,".feather")
  
  out<-NULL
  for(i in files){
    t<-read_feather(i)
    out<-bind_cols(out,t)
  }
  names(out) <- pnames
  return(out)
}

# save era5 data from per-day files to per-point files
saveERAPoint <- function(lonlat,pnames,date_seq){
  long <- as.vector(unlist(lonlat[1]))
  lat <- as.vector(unlist(lonlat[2]))
  eraPoint <- lapply(as.list(date_seq),mergeERAvars,pnames,long,lat)
  eraPoint <- bind_rows(eraPoint)
  dir <- paste0(direrabase,"/ERApoints")
  dir.create(dir,showWarnings = FALSE)
  write_feather(eraPoint,paste(dir,"/",as.character(long),"_",as.character(lat),".feather",sep=""))
}

# read era5 data on one point
getERAPoint <- function(long,lat){
  x <- read_feather(paste(direrabase,"/ERAPoints/",long,"_",lat,".feather",sep=""))
  time <- read_feather(paste(direrabase,"/ERADate.feather",sep=""))
  return(cbind(time,x))
}


