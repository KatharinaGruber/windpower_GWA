# Script for dealing with ERA5 data, similar to Script ERA5_data.R
# but for larger datasets
# area of Brazil was too large to handle it with the existing functions, especially bind_rows led to problems
# either not working at all or being really really slow
# fixed problem by using loops, which is not very efficient, but faster than bind_rows and, most importantly, WORKS!
# data are first aggregated by year, as larger datasets will not fit in memory
# then yearly files are aggregated to all downloaded years and saved per point in .feather format for faster loading



###########################################################################################
######################################FUNCTIONS############################################
###########################################################################################


####reads era5 file from disk
####creates file of lons lats
####saves all points into single feather files
####joining the respective time-period together
####date_seq is a list of vectors of time steps (resolution of downloaded data)
convertERAFeather<- function(date_seq,pname,pname1) {

  ncname_head<-paste(direrabase,"/era5_",pname1,"_component_of_wind_",sep="")
  
  ncfile <- nc_open(paste0(ncname_head,format(date_seq[[1]][1],"%Y%m%d"),".nc"))
  
  
  #Longitude
  longitude <- ncvar_get(ncfile, "longitude", verbose = F)
  nlon <- dim(longitude)
  
  #Latitude
  latitude <- ncvar_get(ncfile, "latitude", verbose = F)
  nlat <- dim(latitude)
  
  lonlat <- expand.grid(longitude,latitude)
  nc_close(ncfile)
  
  r <- vector("list",dim(lonlat)[1])
  
  for(month in c(1:length(date_seq))){
    
    print(format(date_seq[[month]][1],"%Y%m"))  
    dates<-format(date_seq[[month]],"%Y%m%d")
    ncfiles<-paste(ncname_head,dates,".nc",sep="")
    
    res <- sapply(ncfiles,readSingleParam,pname,nlon,nlat,simplify=FALSE)
    
    
    for(j in c(1:dim(res[[1]])[2])){
      for(i in c(1:length(res))){
        r[[j]] <- c(r[[j]],res[[i]][,j])
      }
    }
  }
  
  r1 <- lapply(lapply(r,unlist),as.vector)
  
  
  ### create directory for era5 data points
  dir<-paste0(direrabase,"/feather")
  dir.create(dir, showWarnings = FALSE)
  dir<-paste0(direrabase,"/feather/",paste(pname,
                                           format(date_seq[[1]][1],"%Y%m%d"),
                                           format(last(date_seq[[length(date_seq)]]),"%Y%m%d"),sep="_"))
  dir.create(dir, showWarnings = FALSE)
  
  ###write era5 points for whole period
  for(i in c(1:length(r1))){
    outfile<-paste(dir,"/",lonlat[i,1],"_",lonlat[i,2],".feather",sep="")
    write_feather(data.frame(era=r1[[i]]),outfile)
  }
  
  write_feather(lonlat,paste(dir,"/lonlat.feather",sep=""))
}

######read a single parameter from a single ERA5-ncfile
readSingleParam<-function(ncfileIn,pname,nlon,nlat){
  #print(ncfileIn)
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


