##arthurlutz 20190205
##Script to preprocess grids for demand indicators of water tower units for NGS-Water Tower Index project
##Also calculates downstream P for cells with demand above threshold
rm(list=ls(all=TRUE))
library(raster)
library(maptools)
library(rgdal)
library(ncdf4)
library(sf)
library(tibble)
library(rgeos)
library(lwgeom)

##SETTINGS
base <- "d:\\Dropbox (FutureWater)\\Team\\Projects\\Active\\2019001_NGS_WaterTowers\\data\\"
WTU <- raster(paste(base,"index\\units\\WTU.tif",sep=""))
basins <- raster(paste(base,"index\\units\\basins.tif",sep=""))
downstream <- raster(paste(base,"index\\units\\basins_downstream.tif",sep=""))
WTU_specs <- read.csv(paste(base,"index\\units\\WTU_specs.csv",sep=""),colClasses=c("NULL",NA,NA,NA,NA))
outdir <- paste(base,"WaterDemandsWRI\\processed\\",sep="")
resolution <- 0.05
demand_threshold <- 0.00001 #km3 per cell

#demand input
Dom_use_y <- brick(paste(base,"WaterDemandsWRI\\wri_waterdemand_historical_5min_jul2016\\global_historical_PDomUse_year_millionm3_5min_1960_2014.nc4",sep=""))
Ind_use_y <- brick(paste(base,"WaterDemandsWRI\\wri_waterdemand_historical_5min_jul2016\\global_historical_PIndUse_year_millionm3_5min_1960_2014.nc4",sep=""))
Irr_use_y <- brick(paste(base,"WaterDemandsWRI\\wri_waterdemand_historical_5min_jul2016\\global_historical_PIrrWN_year_millionm3_5min_1960_2014.nc4",sep=""))
EFR <- brick(paste(base,"EnvironmentalFlow\\global_historical_riverdischarge_ymonmean_m3second_5min_2001_2014.nc4",sep=""))

#precipitation input
Pannual <- raster(paste(base,"ERA5\\processed\\P_avg_annual_mm.tif",sep=""))

#evaporation input
ET_avg_annual <- raster(paste(base,"ERA5\\evaporation\\ERA5_evaporation_avgannual_2001_2017.nc",sep=""))

##SETTINGS END

##increase speed of raster operations
rasterOptions(maxmemory = 1e+09)

#create raster template at project resolution
template <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,resolution=c(resolution,resolution))

#create raster with area per cell (km2)
area_km2 <- area(template)

#turn off scientific notation
options(scipen=999)

#shift coordinates ET grid and convert to mm/yr and resample
ET_avg_annual_shift <- shift(ET_avg_annual,x=0.125)
ET1 <- crop(ET_avg_annual_shift , c(xmn=0, xmx=180,ymn=-90,ymx=90))
ET2 <- crop(ET_avg_annual_shift , c(xmn=180, xmx=360,ymn=-90,ymx=90))
ET2_s <- shift(ET2,x=-360)
ET_avg_annual2 <- merge(ET1,ET2_s)*1000
ET_avg_annual_res <- resample(ET_avg_annual2,template,method="bilinear")

#initiate dataframe to store demand values per WTU + dataframe to store available downstream P
no <- length(unique(WTU))
df<-as.data.frame(matrix(ncol=6,nrow=no))
df[,1]<-1:no
colnames(df)<-c("WTU_ID","Dom_demand_basin_km3yr-1","Ind_demand_basin_km3yr-1","Irr_demand_basin_km3yr-1","Nat_demand_basin_km3yr-1","Total_demand_basin_km3yr-1")

no <- length(unique(WTU))
df_P<-as.data.frame(matrix(ncol=6,nrow=no))
df_P[,1]<-1:no
colnames(df_P)<-c("WTU_ID","Dom_P_ds_available_km3yr-1","Ind_P_ds_available_km3yr-1","Irr_P_ds_available_km3yr-1","Nat_P_ds_available_km3yr-1","Total_demand_P_ds_available_km3yr-1")

###DOMESTIC DEMAND INDICATORS
print("Calculating domestic demand indicators...")
print("Calculating and writing average annual grid, resampled to project resolution...")
#subset years (2001-2014)
Dom_use_y_sub <- subset(Dom_use_y,42:55)

##calculate average annual mean and convert to km3
Dom_use_avg_annual <- mean(Dom_use_y_sub)*0.001

#resample and correct totals for change in grid resolution
agf <- res(Dom_use_avg_annual)[1] / res(template)[1]
Dom_use_avg_annual2 <- resample(Dom_use_avg_annual/(agf^2),template,method='ngb')
#writeRaster(Dom_use_avg_annual2,paste(outdir,"Domestic_use_avg_annual_km3.tif",sep=""),overwrite=T)
print("Done")

#aggregate for basin
print("Aggregating for total basin, adding to table and writing as grid...")
Dom_use_sum_basin <- zonal(Dom_use_avg_annual2,basins,fun='sum')
df[,2] <- Dom_use_sum_basin[,2]
Dom_use_basin_grid <- subs(basins,df,by=1,which=2)
writeRaster(Dom_use_basin_grid,paste(outdir,"Domestic_use_avg_annual_basin_km3.tif",sep=""),overwrite=T)
print("Done")

#calculate P available for demand
print("Calculating downstream P available to fulfil demand...")
#convert P and ET to km3 (from mm)
P_km3 <- Pannual * area_km2*0.000001
ET_km3 <- ET_avg_annual_res * area_km2 * 0.000001

#mask for cells with demand above threshold
P_masked <- min(P_km3,max((P_km3+ET_km3),0)) * (Dom_use_avg_annual2 > demand_threshold)
P_available <- zonal(P_masked,downstream,fun='sum')
colnames(P_available)<-c("WTU_ID","P")
temp <- merge(WTU_specs, P_available, by = "WTU_ID",all.x=T)
df_P[,2] <- temp$P
print("Done")

###INDUSTRIAL DEMAND INDICATORS
print("Calculating industrial demand indicators...")
print("Calculating and writing average annual grid, resampled to project resolution...")
#subset years (2001-2014)
Ind_use_y_sub <- subset(Ind_use_y,42:55)

##calculate average annual mean and convert to km3
Ind_use_avg_annual <- mean(Ind_use_y_sub)*0.001

#resample and correct totals for change in grid resolution
agf <- res(Ind_use_avg_annual)[1] / res(template)[1]
Ind_use_avg_annual2 <- resample(Ind_use_avg_annual/(agf^2),template,method='ngb')
#writeRaster(Ind_use_avg_annual2,paste(outdir,"Industrial_use_avg_annual_km3.tif",sep=""),overwrite=T)
print("Done")

#aggregate for basin
print("Aggregating for total basin, adding to table and writing as grid...")
Ind_use_sum_basin <- zonal(Ind_use_avg_annual2,basins,fun='sum')
df[,3] <- Ind_use_sum_basin[,2]
Ind_use_basin_grid <- subs(basins,df,by=1,which=3)
#writeRaster(Ind_use_basin_grid,paste(outdir,"Industrial_use_avg_annual_basin_km3.tif",sep=""),overwrite=T)
print("Done")

#calculate P available for demand
print("Calculating downstream P available to fulfil demand...")
#mask for cells with demand above threshold
P_masked <- min(P_km3,max((P_km3+ET_km3),0)) * (Ind_use_avg_annual2 > demand_threshold)
P_available <- zonal(P_masked,downstream,fun='sum')
colnames(P_available)<-c("WTU_ID","P")
temp <- merge(WTU_specs, P_available, by = "WTU_ID",all.x=T)
df_P[,3] <- temp$P
print("Done")

###IRRIGATION DEMAND INDICATORS
print("Calculating irrigation demand indicators...")
print("Calculating and writing average annual grid, resampled to project resolution...")
#subset years (2001-2014)
Irr_use_y_sub <- subset(Irr_use_y,42:55)

##calculate average annual mean and convert to km3
Irr_use_avg_annual <- mean(Irr_use_y_sub) *0.001

#resample and correct totals for change in grid resolution
agf <- res(Irr_use_avg_annual)[1] / res(template)[1]
Irr_use_avg_annual2 <- resample(Irr_use_avg_annual/(agf^2),template,method='ngb')
writeRaster(Irr_use_avg_annual2,paste(outdir,"Irrigation_use_avg_annual_km3.tif",sep=""),overwrite=T)
print("Done")

#aggregate for basin
print("Aggregating for total basin, adding to table and writing as grid...")
Irr_use_sum_basin <- zonal(Irr_use_avg_annual2,basins,fun='sum')
df[,4] <- Irr_use_sum_basin[,2]
Irr_use_basin_grid <- subs(basins,df,by=1,which=4)
writeRaster(Irr_use_basin_grid,paste(outdir,"Irrigation_use_avg_annual_basin_km3.tif",sep=""),overwrite=T)
print("Done")

#calculate P available for demand
print("Calculating downstream P available to fulfil demand...")
#mask for cells with demand above threshold
P_masked <- min(P_km3,max((P_km3+ET_km3),0)) * (Irr_use_avg_annual2 > demand_threshold)
P_available <- zonal(P_masked,downstream,fun='sum')
colnames(P_available)<-c("WTU_ID","P")
temp <- merge(WTU_specs, P_available, by = "WTU_ID",all.x=T)
df_P[,4] <- temp$P
print("Done")

###NATURAL DEMAND INDICATORS
print("Calculating natural demand indicators...")
print("Calculating average annual grid")

##calculate average annual mean and convert from m3/s to km3 / yr
EFR_avg_annual <- mean(EFR) * 3600*24*365 * 0.001*0.001*0.001

##extract maximum values in downstream basin
#resample basins to EFR resolution
print("Extracting maximum EFR per basin, writing as grid and storing in table...")
EFR_res <- resample(EFR_avg_annual,basins,method="ngb")
#extract maximum values
EFR_max <- zonal(EFR_res,basins,fun='max')
df[,5] <- EFR_max[,2]
Nat_demand_basin_grid <- subs(basins,df,by=1,which=5)
writeRaster(Nat_demand_basin_grid,paste(outdir,"Natural_demand_avg_annual_basin_km3.tif",sep=""),overwrite=T)
print("Done")

#calculate P available for demand
print("Calculating downstream P available to fulfil demand...")
P_available <- zonal(min(P_km3,max((P_km3+ET_km3),0)),downstream,fun='sum')
colnames(P_available)<-c("WTU_ID","P")
temp <- merge(WTU_specs, P_available, by = "WTU_ID",all.x=T)
df_P[,5] <- temp$P
print("Done")

##TOTAL HUMAN DEMAND
print("Sum for total demand and add to table and write as grid...")
#sum for total demand
Total_demand_km3 <- Dom_use_avg_annual2 + Ind_use_avg_annual2 + Irr_use_avg_annual2
df[,6] <- df[,2]+df[,3]+df[,4]
Total_demand_basin_grid <- subs(basins,df,by=1,which=6)
writeRaster(Total_demand_basin_grid,paste(outdir,"Total_human_demand_avg_annual_basin_km3.tif",sep=""),overwrite=T)
print("Done")

#calculate P available for demand
print("Calculating downstream P available to fulfil demand...")
#mask for cells with demand above threshold
P_masked <- P_km3 * (Total_demand_km3 > demand_threshold)
P_available <- zonal(min(P_km3,max((P_km3+ET_km3),0)),downstream,fun='sum')
colnames(P_available)<-c("WTU_ID","P")
temp <- merge(WTU_specs, P_available, by = "WTU_ID",all.x=T)
df_P[,6] <- temp$P
print("Done")

#write tables
print("Writing CSV tables with Demand indicators...")
write.csv(df,paste(outdir,"WTU_Demand_indicators.csv",sep=""))
write.csv(df_P,paste(outdir,"WTU_Demand_DS_P_available.csv",sep=""))
print("Done")