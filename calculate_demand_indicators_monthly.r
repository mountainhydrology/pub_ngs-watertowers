##arthurlutz 20190313
##Script to calculate demand indicators at monthly time scale for NGS-Water Tower Index project
rm(list=ls(all=TRUE))
library(raster)
library(maptools)
library(rgdal)
library(ncdf4)

##SETTINGS
base <- "d:\\Dropbox (FutureWater)\\Team\\Projects\\Active\\2019001_NGS_WaterTowers\\data\\"
WTU <- raster(paste(base,"index\\units\\WTU.tif",sep=""))
basins <- raster(paste(base,"index\\units\\basins.tif",sep=""))
downstream <- raster(paste(base,"index\\units\\basins_downstream.tif",sep=""))
WTU_specs <- read.csv(paste(base,"index\\units\\WTU_specs.csv",sep=""),colClasses=c("NULL",NA,NA,NA,NA))
outdir <- paste(base,"WaterDemandsWRI\\processed\\",sep="")
resolution <- 0.05
calculate_demand_climatologies <- F
resample_era5_climatology <- F
resample_era5_ET_climatology <- F
demand_threshold <- 0.000001 #km3 per cell

#demand input (demand source files are million m3)
Dom_use_m <- brick(paste(base,"WaterDemandsWRI\\wri_waterdemand_historical_5min_jul2016\\global_historical_PDomUse_month_millionm3_5min_1960_2014.nc4",sep=""))
Ind_use_m <- brick(paste(base,"WaterDemandsWRI\\wri_waterdemand_historical_5min_jul2016\\global_historical_PIndUse_month_millionm3_5min_1960_2014.nc4",sep=""))
Irr_use_m <- brick(paste(base,"WaterDemandsWRI\\wri_waterdemand_historical_5min_jul2016\\global_historical_PIrrWN_month_millionm3_5min_1960_2014.nc4",sep=""))

#precipitation input
Pdir <- paste(base,"ERA5\\",sep="")
Pymonmean <- brick(paste(Pdir,"monthly\\era5_total-precipitation_ymonmean_2001-2017_global.tif",sep=""))

#evaporation input
ETymonmean <- brick(paste(base,"ERA5\\evaporation\\ERA5_evaporation_ymonmean_2001_2017.nc",sep=""))

#EFR input
EFR <- brick(paste(base,"EnvironmentalFlow\\global_historical_riverdischarge_ymonmean_m3second_5min_2001_2014.nc4",sep=""))
##SETTINGS END

##increase speed of raster operations
rasterOptions(maxmemory = 1e+09)

#create raster template at project resolution
template <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,resolution=c(resolution,resolution))

#create raster with area per cell (km2)
area_km2 <- area(template)

#turn off scientific notation
options(scipen=999)

if(resample_era5_ET_climatology == T)
{
  print("Resampling ERA 5 monthly climatology to project resolution and writing grids...")
  #shift coordinates ET grid and convert to mm/month and resample
  ET_ymonmean_shift <- shift(ETymonmean,x=0.125)
  ET1 <- crop(ET_ymonmean_shift , c(xmn=0, xmx=180,ymn=-90,ymx=90))
  ET2 <- crop(ET_ymonmean_shift , c(xmn=180, xmx=360,ymn=-90,ymx=90))
  ET2_s <- shift(ET2,x=-360)
  ETymonmean2 <- merge(ET1,ET2_s)*1000
  ETymonmean_res <- resample(ETymonmean2,template,method="bilinear")
  writeRaster(ETymonmean_res,paste(Pdir,"evaporation\\era5_total-evaporation_ymonmean_2001-2017_global_005.tif",sep=""),overwrite=T)
}
##read ERA 5 monthly ET climatology at project resolution
ETymonmean_res <- brick(paste(Pdir,"evaporation\\era5_total-evaporation_ymonmean_2001-2017_global_005.tif",sep=""))

##initiate dataframe to store water gap values per WTU
no <- length(unique(WTU))
df<-as.data.frame(matrix(ncol=14,nrow=no))
df[,1]<-1:no
colnames(df)<-c("WTU_ID","Water_gap_total","Water_gap_Jan","Water_gap_Feb","Water_gap_Mar","Water_gap_Apr","Water_gap_May","Water_gap_Jun","Water_gap_Jul","Water_gap_Aug","Water_gap_Sep","Water_gap_Oct","Water_gap_Nov","Water_gap_Dec")

##create monthly climatologies for demands, if TRUE in settings
if(calculate_demand_climatologies == T)
{
  ##create monthly climatology of demand files, sum demands, resample to project resolution
  print("Create monthly climatology of demand files, sum demands, and resample to project resolution..")
  #climatology per month for each demand (jan 2001 is layer 493, dec 2014 is layer 660)
  Dom_clim <- stack()
  Ind_clim <- stack()
  Irr_clim <- stack()
  for(m in 1:12)
  {
    timesteps <- c()
    for(y in 1:14)
    {
      x <- (492+y*12)-12+m
      timesteps[y]<-x
    }
    timesteps
    Dom_use_m_sub <- subset(Dom_use_m,timesteps)
    Dom_clim <- stack(Dom_clim,mean(Dom_use_m_sub))
    print(paste("Domestic month ",m," done...",sep=""))
    Ind_use_m_sub <- subset(Ind_use_m,timesteps)
    Ind_clim <- stack(Ind_clim,mean(Ind_use_m_sub))
    print(paste("Industrial month ",m," done...",sep=""))
    Irr_use_m_sub <- subset(Irr_use_m,timesteps)
    Irr_clim <- stack(Irr_clim,mean(Irr_use_m_sub))
    print(paste("Irrigation month ",m," done...",sep=""))
  }
  
  #calculate total human demand
  print("Sum to total demand...")
  Tot_clim <- Dom_clim + Ind_clim + Irr_clim
  print("Done")
  
  #correct values for change in resolution, write grids, use gdal_translate to resample to project resolution
  agf <- res(Dom_use_m)[1] / res(template)[1]
  
  print("Resample domestic and write grids...")
  Dom_clim_cor <- Dom_clim/(agf^2)
  writeRaster(Dom_clim_cor,paste(outdir,"demand_monthly_climatology\\Dom_use_ymonmean_2001_2014_5min_temp.tif",sep=""))
  setwd(outdir)
  command <- paste("gdal_translate -tr 0.05 0.05 -r nearest demand_monthly_climatology\\Dom_use_ymonmean_2001_2014_5min.tif demand_monthly_climatology\\Dom_use_ymonmean_2001_2014_005.tif",sep="")
  system(command)
  file.remove(paste(outdir,"demand_monthly_climatology\\Dom_use_ymonmean_2001_2014_5min_temp.tif",sep=""))
  print("Done")
  
  print("Resample industrial and write grids...")
  Ind_clim_cor <- Ind_clim/(agf^2)
  writeRaster(Ind_clim_cor,paste(outdir,"demand_monthly_climatology\\Ind_use_ymonmean_2001_2014_5min_temp.tif",sep=""))
  setwd(outdir)
  command <- paste("gdal_translate -tr 0.05 0.05 -r nearest demand_monthly_climatology\\Ind_use_ymonmean_2001_2014_5min_temp.tif demand_monthly_climatology\\Ind_use_ymonmean_2001_2014_005.tif",sep="")
  system(command)
  file.remove(paste(outdir,"demand_monthly_climatology\\Ind_use_ymonmean_2001_2014_5min_temp.tif",sep=""))
  print("Done")
  
  print("Resample irrigation and write grids...")
  Irr_clim_cor <- Irr_clim/(agf^2)
  writeRaster(Irr_clim_cor,paste(outdir,"demand_monthly_climatology\\Irr_use_ymonmean_2001_2014_5min_temp.tif",sep=""))
  setwd(outdir)
  command <- paste("gdal_translate -tr 0.05 0.05 -r nearest demand_monthly_climatology\\Irr_use_ymonmean_2001_2014_5min_temp.tif demand_monthly_climatology\\Irr_use_ymonmean_2001_2014_005.tif",sep="")
  system(command)
  file.remove(paste(outdir,"demand_monthly_climatology\\Irr_use_ymonmean_2001_2014_5min_temp.tif",sep=""))
  print("Done")
  
  print("Resample total demand and write grids...")
  Tot_clim_cor <- Tot_clim/(agf^2)
  writeRaster(Tot_clim_cor,paste(outdir,"demand_monthly_climatology\\Tot_use_ymonmean_2001_2014_5min_temp.tif",sep=""))
  setwd(outdir)
  command <- paste("gdal_translate -tr 0.05 0.05 -r nearest demand_monthly_climatology\\Tot_use_ymonmean_2001_2014_5min_temp.tif demand_monthly_climatology\\Tot_use_ymonmean_2001_2014_005.tif",sep="")
  system(command)
  file.remove(paste(outdir,"demand_monthly_climatology\\Tot_use_ymonmean_2001_2014_5min_temp.tif",sep=""))
  print("Done")
}

##load demand grids at project resolution
Tot_demand_ymonmean <- brick(paste(outdir,"demand_monthly_climatology\\Tot_use_ymonmean_2001_2014_005.tif",sep=""))
Dom_use_ymonmean <- brick(paste(outdir,"demand_monthly_climatology\\Dom_use_ymonmean_2001_2014_005.tif",sep=""))
Ind_use_ymonmean <- brick(paste(outdir,"demand_monthly_climatology\\Ind_use_ymonmean_2001_2014_005.tif",sep=""))
Irr_use_ymonmean <- brick(paste(outdir,"demand_monthly_climatology\\Irr_use_ymonmean_2001_2014_005.tif",sep=""))

##resample ERA5 P files to project resolution, if needed
if(resample_era5_climatology == T)
{
  print("Resampling ERA 5 monthly climatology to project resolution and writing grids...")
  #shift coordinates
  Pymonmean_shift <- shift(Pymonmean,x=0.125)
  P1 <- crop(Pymonmean_shift, c(xmn=0, xmx=180,ymn=-90,ymx=90))
  P2 <- crop(Pymonmean_shift, c(xmn=180, xmx=360,ymn=-90,ymx=90))
  P2_s <- shift(P2,x=-360)
  Pymonmean2 <- merge(P1,P2_s)
  #resample and convert to mm
  Pymonmean_res <- resample(Pymonmean2*1000,template,method="bilinear")
  writeRaster(Pymonmean_res,paste(Pdir,"monthly\\era5_total-precipitation_ymonmean_2001-2017_global_005.tif",sep=""),overwrite=T)
  print("Done")
}

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

#loop over months, get basin demand, downstream P and water gap per sector per WTU.
print("Calculating water gap per WTU per month, adding to table and writing as grid...")
for(m in 1:12)
{
  #calculate P available for demand by subtracting ET
  print("Calculating downstream P - ET...")
  #convert P and ET to km3 (from mm)
  P_km3 <- subset(Pymonmean2,m) * area_km2*0.000001
  ET_km3 <- subset(ETymonmean_res,m) * area_km2 * 0.000001
  Pnet_km3 <- min(P_km3,max((P_km3+ET_km3),0))
  
  ###DOMESTIC DEMAND INDICATORS
  print("Calculating domestic demand indicators...")
  ##get ymonmean for month m and convert to km3
  Dom_use_ymonmean_m <- subset(Dom_use_ymonmean,m)*0.001
 
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
  Ind_use_ymonmean_m <- subset(Ind_use_ymonmean,m)*0.001
  
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
  Irr_use_ymonmean_m <- subset(Irr_use_ymonmean,m)*0.001
  
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
  ##get natural demand per month
  EFR_m <- EFR[[m]] * 3600*24*30.42 * 0.001*0.001*0.001
  #resample EFR to model resolution
  EFR_m_res <- resample(EFR_m,basins,method="ngb")
  #extract maximum values
  EFR_max <- zonal(EFR_m_res,basins,fun='max')
  colnames(EFR_max) <- c("WTU_ID","EFR")
  #calculate water gap aggregated for basin
  Pnet_basin_ds <- zonal(Pnet_km3,downstream,fun='sum')
  colnames(Pnet_basin_ds)<-c("WTU_ID","Pnet")
  temp <- merge(WTU_specs, Pnet_basin_ds, by = "WTU_ID",all.x=T)
  temp2 <- merge(temp,EFR_max,by = "WTU_ID",all.x=T)
  Nat_gap_m <- temp2$EFR - temp2$Pnet
  Nat_gap_m[Nat_gap_m<0] <- 0
  df_nat[,m+2] <- Nat_gap_m
  print("Done")
  
   #caclulate total water gap per month
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
write.csv(df,paste(outdir,"WTU_Total_Water_Gap_monthly.csv",sep=""))
write.csv(df_dom,paste(outdir,"WTU_Domestic_Water_Gap_monthly.csv",sep=""))
write.csv(df_ind,paste(outdir,"WTU_Industrial_Water_Gap_monthly.csv",sep=""))
write.csv(df_irr,paste(outdir,"WTU_Irrigation_Water_Gap_monthly.csv",sep=""))
write.csv(df_nat,paste(outdir,"WTU_Natural_Water_Gap_monthly.csv",sep=""))
print("Done")