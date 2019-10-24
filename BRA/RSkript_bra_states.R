library(Metrics)
library(zoo)

dirmerrabase <- "C:/..."
direrabase <- "C:/..."
dirons <- "C:/..."
dirwindproddaily <- paste0(dirons,"/daily_data")
dirwindparks <- "C:/..."
dirresults <- "C:/..."

source(paste0(dirresults,"/functions_bra2.R"))


##### prepare results for states #####

# MERRA
load(paste(dirresults,"/results/STATEpowlist_NN.RData",sep=""))
STATEpowlist_NNM <- STATEpowlist
load(paste(dirresults,"/results/STATEpowlist_wsmaWA.RData",sep=""))
STATEpowlist_WAM <- STATEpowlist

# ERA
load(paste0(direrabase,"/results/statpowerWA_bra.RData"))
load(paste0(direrabase,"/results/statpower_bra.RData"))

##### sum up per state for ERA, has not been done before

# load windparkdata to know which state they are in
load(paste0(dirwindparks,"/windparks_complete.RData"))
windparks <- windparks[order(windparks$long),]
windparks <- data.frame(windparks,comdate=as.POSIXct(paste(windparks$year,"-",windparks$month,"-",windparks$day," 00:00:00",sep=""),tz="UTC"))
windparks <- windparks[which(windparks$comdate < as.POSIXct("2017-08-31 00:00:00",tz="UTC")),]
states <- data.frame(num=c(1:length(rle(windparks$long)$lengths)),states=windparks$state[cumsum(rle(windparks$long)$lengths)])
states<- states[order(states$states),]
cumnum = c(0,cumsum(rle(as.vector(states$states))$lengths))
statename = rle(as.vector(states$states))$values

STATEpowlist_NNE <- list()
STATEpowlist_WAE <- list()
for(i in c(1:length(statename))){
  print(statename[i])
  statepowNN <- NULL
  statepowWA <- NULL
  for(j in c((cumnum[i]+1):cumnum[i+1])){
    if(length(statepowNN>0)){
      statepowNN[,2]=statepowNN[,2]+statpower[[states[j,1]]][,2]
      statepowWA[,2]=statepowWA[,2]+statpower_meanAPT[[states[j,1]]][,2]
    }else{
      statepowNN=statpower[[states[j,1]]]
      statepowWA=statpower_meanAPT[[states[j,1]]]
    }
  }
  STATEpowlist_NNE[[i]] <- statepowNN
  STATEpowlist_WAE[[i]] <- statepowWA
}
names(STATEpowlist_NNE) <- statename
names(STATEpowlist_WAE) <- statename

####



# remove states which have no comparison data
no_comp <- c(3,4,5,6,9,13)
for(i in c(no_comp[length(no_comp):1])){
  STATEpowlist_NNM[[i]] <- NULL
  STATEpowlist_WAM[[i]] <- NULL
  STATEpowlist_NNE[[i]] <- NULL
  STATEpowlist_WAE[[i]] <- NULL
}



# aggregate daily
STATEpowlist_NNMd <- list()
STATEpowlist_WAMd <- list()
STATEpowlist_NNEd <- list()
STATEpowlist_WAEd <- list()
for(i in c(1:length(STATEpowlist_NNM))){
  STATEpowlist_NNMd[[i]] <- aggregate(STATEpowlist_NNM[[i]][,2],by=list(format(STATEpowlist_NNM[[i]][,1],"%Y%m%d")),sum)
  STATEpowlist_WAMd[[i]] <- aggregate(STATEpowlist_WAM[[i]][,2],by=list(format(STATEpowlist_WAM[[i]][,1],"%Y%m%d")),sum)
  STATEpowlist_NNEd[[i]] <- aggregate(STATEpowlist_NNE[[i]][,2],by=list(format(STATEpowlist_NNE[[i]][,1],"%Y%m%d")),sum)
  STATEpowlist_WAEd[[i]] <- aggregate(STATEpowlist_WAE[[i]][,2],by=list(format(STATEpowlist_WAE[[i]][,1],"%Y%m%d")),sum)
}

names(STATEpowlist_NNMd) <- names(STATEpowlist_NNM)
names(STATEpowlist_WAMd) <- names(STATEpowlist_WAM)
names(STATEpowlist_NNEd) <- names(STATEpowlist_NNE)
names(STATEpowlist_WAEd) <- names(STATEpowlist_WAE)

save(STATEpowlist_NNMd,STATEpowlist_WAMd,STATEpowlist_NNEd,STATEpowlist_WAEd,file=paste0(dirresults,"/results/simulated_BRA_daily.RData"))




# convert calculated wind power from kWh to GWh
for(i in c(1:length(STATEpowlist_NNEd))){
  STATEpowlist_NNMd[[i]][,2] <- STATEpowlist_NNMd[[i]][,2]/10^6
  STATEpowlist_WAMd[[i]][,2] <- STATEpowlist_WAMd[[i]][,2]/10^6
  STATEpowlist_NNEd[[i]][,2] <- STATEpowlist_NNEd[[i]][,2]/10^6
  STATEpowlist_WAEd[[i]][,2] <- STATEpowlist_WAEd[[i]][,2]/10^6
}




##### load measured wind power ####
# states
STATEprod <- list()

states <- c("Bahia","Ceará","Paraíba","Paraná","Pernambuco","Piaui","Rio de Janeiro","Rio Grande do Norte","Rio Grande do Sul","Santa Catarina","Sergipe")

for(i in c(1:length(STATEpowlist_NNEd))){
  ind <- which(states==names(STATEpowlist_NNEd)[i])
  STATEprod[[i]] = getSTATEproddaily(ind)
  
}

names(STATEprod) <- names(STATEpowlist_NNEd)

########




# cut simulated and observed wind power to same lengths
STATEcomp_NNMd <- list()
STATEcomp_WAMd <- list()
STATEcomp_NNEd <- list()
STATEcomp_WAEd <- list()
for(i in c(1:length(STATEpowlist_NNEd))){
  STATEcomp_NNEd[[i]] <- csl(STATEpowlist_NNEd[[i]],STATEprod[[i]])
  STATEcomp_WAEd[[i]] <- csl(STATEpowlist_WAEd[[i]],STATEprod[[i]])
  STATEcomp_NNMd[[i]] <- csl(STATEpowlist_NNMd[[i]],STATEprod[[i]])
  STATEcomp_WAMd[[i]] <- csl(STATEpowlist_WAMd[[i]],STATEprod[[i]])
}






save(STATEcomp_NNEd,STATEcomp_NNMd,STATEcomp_WAEd,STATEcomp_WAMd,file=paste0(dirresults,"/results/comp_BRA_daily.RData"))







##### analyse results for each state
stats <- data.frame('Country'='Brazil',
                    'State'=rep(names(STATEpowlist_NNEd),each=4),
                    'Data'=rep(c('ERA5','ERA5','MERRA2','MERRA2'),length(STATEpowlist_NNEd)),
                    'GWA'=rep(c(0,1,0,1),length(STATEpowlist_NNEd)),
                    'RMSE'=NA,
                    'MBE'=NA)

for(i in c(1:length(STATEpowlist_NNEd))){
  stats$RMSE[(i-1)*4+1] <- rmse(STATEcomp_NNEd[[i]][,2],STATEcomp_NNEd[[i]][,3]) 
  stats$RMSE[(i-1)*4+2] <- rmse(STATEcomp_WAEd[[i]][,2],STATEcomp_WAEd[[i]][,3]) 
  stats$RMSE[(i-1)*4+3] <- rmse(STATEcomp_NNMd[[i]][,2],STATEcomp_NNMd[[i]][,3]) 
  stats$RMSE[(i-1)*4+4] <- rmse(STATEcomp_WAMd[[i]][,2],STATEcomp_WAMd[[i]][,3])
  
  stats$MBE[(i-1)*4+1] <- mean(STATEcomp_NNEd[[i]][,2]-STATEcomp_NNEd[[i]][,3]) 
  stats$MBE[(i-1)*4+2] <- mean(STATEcomp_WAEd[[i]][,2]-STATEcomp_WAEd[[i]][,3]) 
  stats$MBE[(i-1)*4+3] <- mean(STATEcomp_NNMd[[i]][,2]-STATEcomp_NNMd[[i]][,3]) 
  stats$MBE[(i-1)*4+4] <- mean(STATEcomp_WAMd[[i]][,2]-STATEcomp_WAMd[[i]][,3])
}


write.csv(stats,file=paste0(dirresults,"/results/stats_BRAstates_daily.csv"))








###### Calculate mean capacities brazilian states
meancaps <- data.frame('state' = names(STATEprod),
                       'cap' = NA)
for(i in c(1:length(STATEcomp_NNEd))){
  wp_state <- windparks[which(windparks$state==names(STATEprod)[i]),]
  caps = data.frame('time'=format(wp_state$comdate[order(wp_state$comdate)],"%Y%m%d"),
                    'cap'=wp_state$cap[order(wp_state$comdate)])
  caps = aggregate(caps$cap,by=list(caps$time),sum)
  names(caps) = c("time","cap")
  capsum = data.frame('comdate'=caps$time,
                      'cap'=cumsum(caps$cap)/10^6)
  time = STATEcomp_NNEd[[i]][,1]
  cap = data.frame('time'=time,'cap'=NA)
  # is there installed capacity before first compdate?
  if(as.vector(capsum[1,1])<=as.vector(cap[1,1])){
    # if yes set first date to last date before that date
    cap[1,2] = tail(capsum$cap[which(as.vector(capsum$comdate)<=as.vector(cap$time[1]))],1)
    # and then cut the remaining capacities
    capsum = capsum[which(as.vector(capsum$comdate)>as.vector(cap$time[1])),]
  }else{
    # starting capacity is 0
    cap[1,2] = 0
  }
  
  # fill in capacities
  cap[match(capsum$comdate,cap$time),2] = capsum[,2]
  # fill in rest of days with last capapcity
  capall <- cap
  capall[,2] <- na.locf(cap[,2])
  meancaps$cap[i] <- mean(capall[,2])
}

save(meancaps,file=paste0(dirresults,"/results/meancaps_states_bra.RData"))





#### RELATIVE RESULTS

stats_rel <- stats
# *24 to account for 24 hours of day and *100 to get in %
stats_rel$RMSE <- stats$RMSE/(meancaps$cap[match(stats$State,meancaps$state)]*24)*100
stats_rel$MBE <- stats$MBE/(meancaps$cap[match(stats$State,meancaps$state)]*24)*100


write.csv(stats_rel,file=paste0(dirresults,"/results/stats_BRAstates_daily_rel.csv"))









####################################################################################
# Monthly results
####################################################################################

load(paste0(dirresults,"/results/comp_BRA_daily.RData"))

STATEcomp_NNEm <- list()
STATEcomp_NNMm <- list()
STATEcomp_WAEm <- list()
STATEcomp_WAMm <- list()

for(i in c(1:length(STATEcomp_NNEd))){
  STATEcomp_NNEm[[i]] <- aggregate(STATEcomp_NNEd[[i]][,2:3],by=list(substr(STATEcomp_NNEd[[i]][,1],1,6)),sum)
  names(STATEcomp_NNEm[[i]]) <- c('time','sim','obs')
  STATEcomp_NNMm[[i]] <- aggregate(STATEcomp_NNMd[[i]][,2:3],by=list(substr(STATEcomp_NNMd[[i]][,1],1,6)),sum)
  names(STATEcomp_NNMm[[i]]) <- c('time','sim','obs')
  STATEcomp_WAEm[[i]] <- aggregate(STATEcomp_WAEd[[i]][,2:3],by=list(substr(STATEcomp_WAEd[[i]][,1],1,6)),sum)
  names(STATEcomp_WAEm[[i]]) <- c('time','sim','obs')
  STATEcomp_WAMm[[i]] <- aggregate(STATEcomp_WAMd[[i]][,2:3],by=list(substr(STATEcomp_WAMd[[i]][,1],1,6)),sum)
  names(STATEcomp_WAMm[[i]]) <- c('time','sim','obs')
}

names(STATEcomp_NNEm) <- as.vector(meancaps$state)
names(STATEcomp_NNMm) <- as.vector(meancaps$state)
names(STATEcomp_WAEm) <- as.vector(meancaps$state)
names(STATEcomp_WAMm) <- as.vector(meancaps$state)

save(STATEcomp_NNEm,STATEcomp_NNMm,STATEcomp_WAEm,STATEcomp_WAMm,file=paste0(dirresults,"/results/comp_BRA_monthly.RData"))


##### analyse results for each state
stats_m <- data.frame('Country'='Brazil',
                      'State'=rep(names(STATEcomp_NNEm),each=4),
                      'Data'=rep(c('ERA5','ERA5','MERRA2','MERRA2'),length(STATEcomp_NNEm)),
                      'GWA'=rep(c(0,1,0,1),length(STATEcomp_NNEm)),
                      'RMSE'=NA,
                      'MBE'=NA)

for(i in c(1:length(STATEcomp_NNEm))){
  stats_m$RMSE[(i-1)*4+1] <- rmse(STATEcomp_NNEm[[i]][,2],STATEcomp_NNEm[[i]][,3]) 
  stats_m$RMSE[(i-1)*4+2] <- rmse(STATEcomp_WAEm[[i]][,2],STATEcomp_WAEm[[i]][,3]) 
  stats_m$RMSE[(i-1)*4+3] <- rmse(STATEcomp_NNMm[[i]][,2],STATEcomp_NNMm[[i]][,3]) 
  stats_m$RMSE[(i-1)*4+4] <- rmse(STATEcomp_WAMm[[i]][,2],STATEcomp_WAMm[[i]][,3])
  
  stats_m$MBE[(i-1)*4+1] <- mean(STATEcomp_NNEm[[i]][,2]-STATEcomp_NNEm[[i]][,3]) 
  stats_m$MBE[(i-1)*4+2] <- mean(STATEcomp_WAEm[[i]][,2]-STATEcomp_WAEm[[i]][,3]) 
  stats_m$MBE[(i-1)*4+3] <- mean(STATEcomp_NNMm[[i]][,2]-STATEcomp_NNMm[[i]][,3]) 
  stats_m$MBE[(i-1)*4+4] <- mean(STATEcomp_WAMm[[i]][,2]-STATEcomp_WAMm[[i]][,3])
}


write.csv(stats_m,file=paste0(dirresults,"/results/stats_BRAstates_monthly.csv"))


#### RELATIVE RESULTS

stats_mrel <- stats_m
# *24 to account for 24 hours of day, *30 to account for 30 days in month (approx.) and *100 to get in %
stats_mrel$RMSE <- stats_m$RMSE/(meancaps$cap[match(stats_m$State,meancaps$state)]*24*30)*100
stats_mrel$MBE <- stats_m$MBE/(meancaps$cap[match(stats_m$State,meancaps$state)]*24*30)*100


write.csv(stats_mrel,file=paste0(dirresults,"/results/stats_BRAstates_monthly_rel.csv"))
