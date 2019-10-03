##arthurlutz 20190729
##Script to calculate the uncertainty in water gap calculation. Calculate water gap monthly. Also write annual total demands for uncertainty analysis
rm(list=ls(all=TRUE))
library(raster)
library(maptools)
library(rgdal)
library(ncdf4)
library(foreach)
library(doParallel)

##SETTINGS
base <- "e:\\Dropbox (FutureWater)\\Team\\Projects\\Active\\2019001_NGS_WaterTowers\\data\\"
WTU <- raster(paste(base,"index\\units\\WTU.tif",sep=""))
basins <- raster(paste(base,"index\\units\\basins.tif",sep=""))
downstream <- raster(paste(base,"index\\units\\basins_downstream.tif",sep=""))
WTU_specs <- read.csv(paste(base,"index\\units\\WTU_specs.csv",sep=""),colClasses=c("NULL",NA,NA,NA,NA))
outdir <- paste("c:\\workdir\\",sep="")
resolution <- 0.05
demand_threshold <- 0.000001 #km3 per cell

##Monthly input
#demand input (demand source files are million m3)
Dom_use_m <- brick(paste(base,"WaterDemandsWRI\\wri_waterdemand_historical_5min_jul2016\\global_historical_PDomUse_month_millionm3_5min_1960_2014.nc4",sep=""))
Ind_use_m <- brick(paste(base,"WaterDemandsWRI\\wri_waterdemand_historical_5min_jul2016\\global_historical_PIndUse_month_millionm3_5min_1960_2014.nc4",sep=""))
Irr_use_m <- brick(paste(base,"WaterDemandsWRI\\wri_waterdemand_historical_5min_jul2016\\global_historical_PIrrWN_month_millionm3_5min_1960_2014.nc4",sep=""))

#precipitation input
Pdir <- paste(base,"ERA5\\",sep="")
Pymonmean <- brick(paste(Pdir,"monthly\\era5_total-precipitation_ymonmean_2001-2017_global.tif",sep=""))

#evaporation input
ETymonmean <- brick(paste(Pdir,"evaporation\\ERA5_evaporation_ymonmean_2001_2017.nc",sep=""))

#EFR input
EFR <- brick(paste(base,"EnvironmentalFlow\\global_historical_riverdischarge_ymonmean_m3second_5min_2001_2014.nc4",sep=""))

###Annual demand input
demands <- read.csv(paste(base,"WaterDemandsWRI\\processed\\WTU_Demand_indicators_v3.csv",sep=""))

#uncertainty input
P_un_sd_wtu <- read.csv(paste(base,"uncertainty\\precipitation\\P_uncertainty_per_WTU.csv",sep=""))
P_un_sd_ds <- read.csv(paste(base,"uncertainty\\precipitation\\P_uncertainty_per_downstream.csv",sep=""))
ET_un_sd_ds <- read.csv(paste(base,"uncertainty\\evaporation\\ET_uncertainty_per_downstream.csv",sep=""))
dem_dom_un_sd <- 0.063
dem_ind_un_sd <- 0.052
dem_irr_un_sd <- 0.198

#number of iteration for uncertainty analysis
nstart <- 1
nend <- 100
cores_unused <- 2
##SETTINGS END

##increase speed of raster operations
rasterOptions(maxmemory = 1e+09)

#create raster template at project resolution
template <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,resolution=c(resolution,resolution))

#create raster with area per cell (km2)
area_km2 <- area(template)

#turn off scientific notation
options(scipen=999)

##read ERA 5 monthly ET climatology at project resolution
ETymonmean_res <- brick(paste(Pdir,"evaporation\\era5_total-evaporation_ymonmean_2001-2017_global_005.tif",sep=""))

#initiate dataframe to store water gap values per WTU
no <- length(unique(WTU))
df<-as.data.frame(matrix(ncol=14,nrow=no))
df[,1]<-1:no
colnames(df)<-c("WTU_ID","Water_gap_total","Water_gap_Jan","Water_gap_Feb","Water_gap_Mar","Water_gap_Apr","Water_gap_May","Water_gap_Jun","Water_gap_Jul","Water_gap_Aug","Water_gap_Sep","Water_gap_Oct","Water_gap_Nov","Water_gap_Dec")

##load demand grids at project resolution
Tot_demand_ymonmean <- brick(paste(base,"WaterDemandsWRI\\processed\\demand_monthly_climatology\\Tot_use_ymonmean_2001_2014_005.tif",sep=""))
Dom_use_ymonmean <- brick(paste(base,"WaterDemandsWRI\\processed\\demand_monthly_climatology\\Dom_use_ymonmean_2001_2014_005.tif",sep=""))
Ind_use_ymonmean <- brick(paste(base,"WaterDemandsWRI\\processed\\demand_monthly_climatology\\Ind_use_ymonmean_2001_2014_005.tif",sep=""))
Irr_use_ymonmean <- brick(paste(base,"WaterDemandsWRI\\processed\\demand_monthly_climatology\\Irr_use_ymonmean_2001_2014_005.tif",sep=""))

##read ERA 5 monthly P climatology at project resolution
Pymonmean2 <- brick(paste(Pdir,"monthly\\era5_total-precipitation_ymonmean_2001-2017_global_005.tif",sep=""))

##initiate dataframes for sectoral water gaps. 
no <- length(unique(WTU))
df_dom<-as.data.frame(matrix(ncol=14,nrow=no))
df_dom[,1]<-1:no
colnames(df_dom)<-c("WTU_ID","Water_gap_total","Water_gap_Jan","Water_gap_Feb","Water_gap_Mar","Water_gap_Apr","Water_gap_May","Water_gap_Jun","Water_gap_Jul","Water_gap_Aug","Water_gap_Sep","Water_gap_Oct","Water_gap_Nov","Water_gap_Dec")

df_ind<-as.data.frame(matrix(ncol=14,nrow=no))
df_ind[,1]<-1:no
colnames(df_ind)<-c("WTU_ID","Water_gap_total","Water_gap_Jan","Water_gap_Feb","Water_gap_Mar","Water_gap_Apr","Water_gap_May","Water_gap_Jun","Water_gap_Jul","Water_gap_Aug","Water_gap_Sep","Water_gap_Oct","Water_gap_Nov","Water_gap_Dec")

df_irr<-as.data.frame(matrix(ncol=14,nrow=no))
df_irr[,1]<-1:no
colnames(df_irr)<-c("WTU_ID","Water_gap_total","Water_gap_Jan","Water_gap_Feb","Water_gap_Mar","Water_gap_Apr","Water_gap_May","Water_gap_Jun","Water_gap_Jul","Water_gap_Aug","Water_gap_Sep","Water_gap_Oct","Water_gap_Nov","Water_gap_Dec")

df_nat<-as.data.frame(matrix(ncol=14,nrow=no))
df_nat[,1]<-1:no
colnames(df_nat)<-c("WTU_ID","Water_gap_total","Water_gap_Jan","Water_gap_Feb","Water_gap_Mar","Water_gap_Apr","Water_gap_May","Water_gap_Jun","Water_gap_Jul","Water_gap_Aug","Water_gap_Sep","Water_gap_Oct","Water_gap_Nov","Water_gap_Dec")

##get natural demand per month (stays constant for all water gap runs)
df_nat_dem <- as.data.frame(matrix(ncol=13,nrow=no))
df_nat_dem[,1]<- 1:no
for(m in 1:12)
{
  EFR_m <- EFR[[m]] * 3600*24*30.42 * 0.001*0.001*0.001
  #resample EFR to model resolution
  EFR_m_res <- resample(EFR_m,basins,method="ngb")
  #extract maximum values
  EFR_max <- zonal(EFR_m_res,basins,fun='max')
  df_nat_dem[,m+1]<-EFR_max[,2]
}
colnames(df_nat_dem)<-c('WTU_ID',1,2,3,4,5,6,7,8,9,10,11,12)

start_time <- Sys.time()

#loop n times for uncertainty analysis
cores = detectCores()
cl = makeCluster(cores - cores_unused)
registerDoParallel(cl)
foreach(i = nstart:nend,.packages = c('base','raster','maptools','rgdal','ncdf4','doParallel','foreach'),.verbose=F) %dopar% {
  
  #creates unique filepath for temp directory
  dir.create(file.path(paste("c:/temp/temp",sprintf("%05d",i),sep="")))
  #sets temp directory
  rasterOptions(tmpdir=file.path(paste("c:/temp/temp",sprintf("%05d",i),sep="")))
             
  ##prepare maps and constants with multiplication factors for uncertainty
  mat <- as.matrix(cbind(P_un_sd_ds$WTU_ID,rnorm(nrow(P_un_sd_ds),1,P_un_sd_ds$sd_normalized[1])))
  #write Precipitation multiplication factor downstream for later use in uncertainty_analysis.R
  colnames(mat)<-c("WTU_ID","P_downstream_mult")
  write.csv(mat,paste(outdir,sprintf("%05d",i),"_P_downstream_multiplication_factor.csv",sep=""))
  
  P_mulfactor_ds <- raster::reclassify(downstream,mat)
  mat <- as.matrix(cbind(ET_un_sd_ds$WTU_ID,rnorm(nrow(ET_un_sd_ds),1,ET_un_sd_ds$sd_normalized[1])))
  ET_mulfactor_ds <- raster::reclassify(downstream,mat)
  dem_dom_mulfactor <- rnorm(1,1,dem_dom_un_sd)
  dem_ind_mulfactor <- rnorm(1,1,dem_ind_un_sd)
  dem_irr_mulfactor <- rnorm(1,1,dem_irr_un_sd)
  EFR_mulfactor <- 1

  ##multiply inputs with factors
  Pymonmean_mult <- Pymonmean2 * P_mulfactor_ds
  ETymonmean_mult <- ETymonmean_res * ET_mulfactor_ds
  Dom_use_ymonmean_mult <- Dom_use_ymonmean * dem_dom_mulfactor
  Ind_use_ymonmean_mult <- Ind_use_ymonmean * dem_ind_mulfactor
  Irr_use_ymonmean_mult <- Irr_use_ymonmean * dem_irr_mulfactor
  
  ##calculate file with annual total demands and write table
  demands2 <- demands
  demands2$Dom_demand_basin_km3yr.1 <- demands$Dom_demand_basin_km3yr.1 * dem_dom_mulfactor
  demands2$Ind_demand_basin_km3yr.1 <- demands$Ind_demand_basin_km3yr.1 * dem_ind_mulfactor
  demands2$Irr_demand_basin_km3yr.1 <- demands$Irr_demand_basin_km3yr.1 * dem_irr_mulfactor
  demands2$Nat_demand_basin_km3yr.1 <- demands$Nat_demand_basin_km3yr.1 * EFR_mulfactor
  demands2$Total_demand_basin_km3yr.1 <- demands2$Dom_demand_basin_km3yr.1 + demands2$Ind_demand_basin_km3yr.1 + demands2$Irr_demand_basin_km3yr.1
  write.csv(demands2,paste(outdir,sprintf("%05d",i),"_WTU_Demand_indicators.csv",sep=""))
  
  #loop over months, get basin demand, downstream P and water gap per sector per WTU.
  print("Calculating water gap per WTU per month, adding to table and writing as grid...")
  for(m in 1:12)
  {
    #calculate P available for demand by subtracting ET
    print("Calculating downstream P - ET...")
    #convert P and ET to km3 (from mm)
    P_km3 <- subset(Pymonmean_mult,m) * area_km2*0.000001
    ET_km3 <- subset(ETymonmean_mult,m) * area_km2 * 0.000001
    Pnet_km3 <- min(P_km3,max((P_km3+ET_km3),0))

    ###DOMESTIC DEMAND INDICATORS
    print("Calculating domestic demand indicators...")
    ##get ymonmean for month m and convert to km3
    Dom_use_ymonmean_m <- subset(Dom_use_ymonmean_mult,m)*0.001

    #calculate water gap at grid cell basis
    Dom_gap_m <- max(Dom_use_ymonmean_m - Pnet_km3,0)

    #mask for cells with demand above threshold
    Dom_gap_m_masked <- Dom_gap_m * (Dom_use_ymonmean_m > demand_threshold)

    #aggregate gap to basin
    Dom_gap_m_basin <- zonal(Dom_gap_m_masked,downstream,fun='sum')
    colnames(Dom_gap_m_basin)<-c("WTU_ID","P")
    temp <- merge(WTU_specs, Dom_gap_m_basin, by = "WTU_ID",all.x=T)
    df_dom[,m+2] <- temp$P
    print("Done")

    ###INDUSTRIAL DEMAND INDICATORS
    print("Calculating industrial demand indicators...")
    ##get ymonmean for month m and convert to km3
    Ind_use_ymonmean_m <- subset(Ind_use_ymonmean_mult,m)*0.001

    #calculate water gap at grid cell basis
    Ind_gap_m <- max(Ind_use_ymonmean_m - Pnet_km3,0)

    #mask for cells with demand above threshold
    Ind_gap_m_masked <- Ind_gap_m * (Ind_use_ymonmean_m > demand_threshold)

    #aggregate gap to basin
    Ind_gap_m_basin <- zonal(Ind_gap_m_masked,downstream,fun='sum')
    colnames(Ind_gap_m_basin)<-c("WTU_ID","P")
    temp <- merge(WTU_specs, Ind_gap_m_basin, by = "WTU_ID",all.x=T)
    df_ind[,m+2] <- temp$P
    print("Done")

    ###IRRIGATION DEMAND INDICATORS
    print("Calculating irrigation demand indicators...")
    ##get ymonmean for month m and convert to km3
    Irr_use_ymonmean_m <- subset(Irr_use_ymonmean_mult,m)*0.001

    #calculate water gap at grid cell basis
    Irr_gap_m <- max(Irr_use_ymonmean_m - Pnet_km3,0)

    #mask for cells with demand above threshold
    Irr_gap_m_masked <- Irr_gap_m * (Irr_use_ymonmean_m > demand_threshold)

    #aggregate gap to basin
    Irr_gap_m_basin <- zonal(Irr_gap_m_masked,downstream,fun='sum')
    colnames(Irr_gap_m_basin)<-c("WTU_ID","P")
    temp <- merge(WTU_specs, Irr_gap_m_basin, by = "WTU_ID",all.x=T)
    df_irr[,m+2] <- temp$P
    print("Done")

    ##NATURAL DEMAND INDICATORS
    #calculate water gap aggregated for basin
    Pnet_basin_ds <- zonal(Pnet_km3,downstream,fun='sum')
    colnames(Pnet_basin_ds)<-c("WTU_ID","Pnet")
    temp <- merge(WTU_specs, Pnet_basin_ds, by = "WTU_ID",all.x=T)
    EFR <- cbind(df_nat_dem[,1],df_nat_dem[,m+1])
    colnames(EFR)<-c("WTU_ID","EFR")
    temp2 <- merge(temp,EFR,by = "WTU_ID",all.x=T)
    Nat_gap_m <- temp2$EFR - temp2$Pnet
    Nat_gap_m[Nat_gap_m<0] <- 0
    df_nat[,m+2] <- Nat_gap_m
    print("Done")

    #calculate total water gap per month
    df[,m+2] <- df_dom[,m+2]+df_ind[,m+2]+df_irr[,m+2]+df_nat[,m+2]

    print(paste("Month ",m," done",sep=""))

  }
  print("Done")

  ##sum gaps to annual
  df[,2]<- df[,3]+df[,4]+df[,5]+df[,6]+df[,7]+df[,8]+df[,9]+df[,10]+df[,11]+df[,12]+df[,13]+df[,14]
  df_dom[,2]<- df_dom[,3]+df_dom[,4]+df_dom[,5]+df_dom[,6]+df_dom[,7]+df_dom[,8]+df_dom[,9]+df_dom[,10]+df_dom[,11]+df_dom[,12]+df_dom[,13]+df_dom[,14]
  df_ind[,2]<- df_ind[,3]+df_ind[,4]+df_ind[,5]+df_ind[,6]+df_ind[,7]+df_ind[,8]+df_ind[,9]+df_ind[,10]+df_ind[,11]+df_ind[,12]+df_ind[,13]+df_ind[,14]
  df_irr[,2]<- df_irr[,3]+df_irr[,4]+df_irr[,5]+df_irr[,6]+df_irr[,7]+df_irr[,8]+df_irr[,9]+df_irr[,10]+df_irr[,11]+df_irr[,12]+df_irr[,13]+df_irr[,14]
  df_nat[,2]<- df_nat[,3]+df_nat[,4]+df_nat[,5]+df_nat[,6]+df_nat[,7]+df_nat[,8]+df_nat[,9]+df_nat[,10]+df_nat[,11]+df_nat[,12]+df_nat[,13]+df_nat[,14]

  #write tables
  write.csv(df,paste(outdir,sprintf("%05d",i),"_WTU_Total_Water_Gap_monthly.csv",sep=""))
  write.csv(df_dom,paste(outdir,sprintf("%05d",i),"_WTU_Domestic_Water_Gap_monthly.csv",sep=""))
  write.csv(df_ind,paste(outdir,sprintf("%05d",i),"_WTU_Industrial_Water_Gap_monthly.csv",sep=""))
  write.csv(df_irr,paste(outdir,sprintf("%05d",i),"_WTU_Irrigation_Water_Gap_monthly.csv",sep=""))
  write.csv(df_nat,paste(outdir,sprintf("%05d",i),"_WTU_Natural_Water_Gap_monthly.csv",sep=""))

  #removes entire temp directory without affecting other running processes
  system(paste0("rm -r c:/temp/temp",sprintf("%05d",i)))
  }

stopCluster(cl)
end_time <- Sys.time()
print(end_time-start_time)