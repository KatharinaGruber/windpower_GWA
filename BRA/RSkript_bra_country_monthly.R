
# compare wind power generation in Brazil between ERA5 and MERRA-2 with and without GWA correction


# base directory for simulation
direrabase <- "C:/..."
dirmerrabase <- "C:/..."
# directory where observed wind power generation data from ONS are stored
dirons <- paste0("C:/...")
# directory for ICEM
dirresults <- "C:/..."


# load results ERA5
load(paste0(direrabase,"/results/wp_brasil.RData"))
wpbra_ERA5 <- wpbra
wpbra_ERA5gwa <- wpbra_meanAPT
# load results MERRA2
load(paste0(dirmerrabase,"/results/comp_wsma.RData"))
wpbra_MERRA2 <- Bpowlist_NNd
wpbra_MERRA2gwa <- Bpowlist_WAd
names(wpbra_MERRA2) <- c("time","wp_kwh")
names(wpbra_MERRA2gwa) <- c("time","wp_kwh")
prod_ons <- Bprod
rm(wpbra,wpbra_meanAPT,statpowlist_NNd,statpowlist_INd,statpowlist_WAd,STATEpowlist_NNd,STATEpowlist_INd,STATEpowlist_WAd,SUBpowlist_NNd,SUBpowlist_INd,SUBpowlist_WAd,Bpowlist_NNd,Bpowlist_INd,Bpowlist_WAd,statprod,STATEprod,SUBprod,Bprod)


# correct capacities
load(paste0(dirmerrabase,"/capacidades_subsistemas_ONS/cap_cfs.RData"))
rm(cfNE,cfS)
wpbra_ERA5$wp_kwh <- wpbra_ERA5$wp_kwh*cfB/10^6
wpbra_ERA5gwa$wp_kwh <- wpbra_ERA5gwa$wp_kwh*cfB/10^6
wpbra_MERRA2$wp_kwh <- wpbra_MERRA2$wp_kwh*cfB/10^6
wpbra_MERRA2gwa$wp_kwh <- wpbra_MERRA2gwa$wp_kwh*cfB/10^6
names(wpbra_ERA5) <- c("time","wp_GWh")
names(wpbra_ERA5gwa) <- c("time","wp_GWh")
names(wpbra_MERRA2) <- c("time","wp_GWh")
names(wpbra_MERRA2gwa) <- c("time","wp_GWh")




# MONTHLY
# aggregate monthly for a monthly comparison
wpbra_ERA5_m <- aggregate(wpbra_ERA5$wp_GWh,by=list(format(wpbra_ERA5$time,"%Y%m")),sum)
wpbra_ERA5gwa_m <- aggregate(wpbra_ERA5gwa$wp_GWh,by=list(format(wpbra_ERA5gwa$time,"%Y%m")),sum)
wpbra_MERRA2_m <- aggregate(wpbra_MERRA2$wp_GWh,by=list(substr(wpbra_MERRA2$time,1,6)),sum)
wpbra_MERRA2gwa_m <- aggregate(wpbra_MERRA2gwa$wp_GWh,by=list(substr(wpbra_MERRA2gwa$time,1,6)),sum)
prod_ons_m <- aggregate(prod_ons$prod_GWh,by=list(substr(prod_ons$date,1,6)),sum)


# cut to right time spans
wpbra_ERA5_m <- wpbra_ERA5_m[which((wpbra_ERA5_m[,1]>=prod_ons_m[1,1])&(wpbra_ERA5_m[,1]<=tail(prod_ons_m[,1],1))),]
wpbra_ERA5gwa_m <- wpbra_ERA5gwa_m[which((wpbra_ERA5gwa_m[,1]>=prod_ons_m[1,1])&(wpbra_ERA5gwa_m[,1]<=tail(prod_ons_m[,1],1))),]
wpbra_MERRA2_m <- wpbra_MERRA2_m[which((wpbra_MERRA2_m[,1]>=prod_ons_m[1,1])&(wpbra_MERRA2_m[,1]<=tail(prod_ons_m[,1],1))),]
wpbra_MERRA2gwa_m <- wpbra_MERRA2gwa_m[which((wpbra_MERRA2gwa_m[,1]>=prod_ons_m[1,1])&(wpbra_MERRA2gwa_m[,1]<=tail(prod_ons_m[,1],1))),]

# remove first month because just two days of production
comp_monthly <- data.frame(time=seq(as.POSIXct("2006-04-01",tz="UTC"),as.POSIXct("2017-08-01",tz="UTC"),by="month"),
                           obs=prod_ons_m[-1,2],
                           era5 = wpbra_ERA5_m[-1,2],
                           era5gwa = wpbra_ERA5gwa_m[-1,2],
                           merra2 = wpbra_MERRA2_m[-1,2],
                           merra2gwa = wpbra_MERRA2gwa_m[-1,2])

save(comp_monthly, file=paste0(dirresults,"/monthly_comp_bra.RData"))




# calculate statistics
stats <- data.frame(type=c("RMSE","MBE"),
                    era5=NA,
                    era5gwa=NA,
                    merra2=NA,
                    merra2gwa=NA)

# RMSE
stats[1,2:5] <- c(rmse(comp_monthly$obs,comp_monthly$era5),
                  rmse(comp_monthly$obs,comp_monthly$era5gwa),
                  rmse(comp_monthly$obs,comp_monthly$merra2),
                  rmse(comp_monthly$obs,comp_monthly$merra2gwa))
#MBE
stats[2,2:5] <- c(mean(comp_monthly$era5-comp_monthly$obs),
                  mean(comp_monthly$era5gwa-comp_monthly$obs),
                  mean(comp_monthly$merra2-comp_monthly$obs),
                  mean(comp_monthly$merra2gwa-comp_monthly$obs))


write.table(stats,file=paste0(dirresults,"/stats_BRAm_abs.csv"),sep=";")



#######



# mean installed capacity
load(paste0(dirwindparks,"/windparks_complete.RData"))
windparks <- windparks[order(windparks$long),]
windparks <- data.frame(windparks,comdate=as.POSIXct(paste(windparks$year,"-",windparks$month,"-",windparks$day," 00:00:00",sep=""),tz="UTC"))
windparks <- windparks[which(windparks$comdate < as.POSIXct("2017-08-31 00:00:00",tz="UTC")),]



cap = aggregate(windparks$cap/10^6,by=list(format(windparks$comdate,"%Y%m")),sum)

cap$date = as.POSIXct(paste0(substr(cap[,1],1,4),"-",substr(cap[,1],5,6),"-01"),tz="UTC")

cap$capsum = cumsum(cap[,2])                

capslist = data.frame(date = comp_monthly$time,
                      cap = NA)
# get first capacity in period
capslist$cap[1] = tail(cap$capsum[cap$date<comp_monthly$time[1]],1)
# cut capacities after
cap <- cap[cap$date>=comp_monthly$time[1],]
# fill in capacities
capslist$cap[match(cap$date,capslist$date)] = cap$capsum
# fill nas
capslist$cap = na.locf(capslist$cap)

# calculate mean capacity
meancap = mean(capslist$cap)







# relative results
stats_r <- stats
stats_r[,2:5] <- round(stats[,2:5]/(meancap*24*30)*100,1)


write.table(stats_r,file=paste0(dirresults,"/stats_BRAm_rel.csv"),sep=";")
