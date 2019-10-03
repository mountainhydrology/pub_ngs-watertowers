##arthurlutz 20190124
##Script to preprocess grids for precipitation indicators of water tower units for NGS-Water Tower Index project
rm(list=ls(all=TRUE))
library(raster)
library(maptools)
library(rgdal)

##SETTINGS
base <- "e:\\Dropbox (FutureWater)\\Team\\Projects\\Active\\2019001_NGS_WaterTowers\\data\\"
WTU <- raster(paste(base,"index\\units\\WTU.tif",sep=""))
basins <- raster(paste(base,"index\\units\\basins.tif",sep=""))
downstream <- raster(paste(base,"index\\units\\basins_downstream.tif",sep=""))
WTU_specs <- read.csv(paste(base,"index\\units_v3\\wTU_specs.csv",sep=""),colClasses=c("NULL",NA,NA,NA,NA))
outdir <- paste(base,"ERA5\\processed\\",sep="")
resolution <- 0.05

#precipitation input
Pannual <- brick(paste(base,"ERA5\\yearly\\era5_total-precipitation_yearsum_2001-2017.tif",sep=""))
Pymonmean <- brick(paste(base,"ERA5\\monthly\\era5_total-precipitation_ymonmean_2001-2017_global.tif",sep=""))

##SETTINGS END

##increase speed of raster operations
rasterOptions(maxmemory = 1e+09)

#create raster template at project resolution
template <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,resolution=c(resolution,resolution))

#create raster with area per cell (km2)
area_km2 <- area(template)

#turn off scientific notation
options(scipen=999)

###PRECIPITATION INDICATORS
#initiate dataframe to store precipitation values per WTU
no <- length(unique(WTU))
df<-as.data.frame(matrix(ncol=13,nrow=no))
df[,1]<-1:no
colnames(df)<-c("WTU_ID","Pannual_WT_mm","Pannual_basin_downstream_mm","Pannual_basin_mm","Pannual_WT_km3","Pannual_basin_downstream_km3","Pannual_basin_km3","Var_inter_upstream","Var_inter_downstream","Var_inter_total","Var_intra_upstream","Var_intra_downstream","Var_intra_total")

#shift coordinates
Pannual_shift <- shift(Pannual,x=0.125)
P1 <- crop(Pannual_shift, c(xmn=0, xmx=180,ymn=-90,ymx=90))
P2 <- crop(Pannual_shift, c(xmn=180, xmx=360,ymn=-90,ymx=90))
P2_s <- shift(P2,x=-360)
Pannual2 <- merge(P1,P2_s)

##CALCULATE AVERAGE ANNUAL TOTAL P
print("Calculating average annual total P...")
Pavg_annual <- mean(Pannual2)*1000

#resample to project resolution
print("Resampling to project resolution and writing global grid...")
Pavg_annual2 <- resample(Pavg_annual,template,method="bilinear")

#write grid
writeRaster(Pavg_annual2,paste(outdir,"P_avg_annual_mm.tif",sep=""),overwrite=T)
print("Done")

#read grid of resampled average annual precipitation
Pavg_annual2 <- raster(paste(outdir,"P_avg_annual_mm.tif",sep=""))

#km3 per grid cell, weighting for area
Pavg_annual_km3 <- Pavg_annual2*area_km2*0.000001

#zonal statistics over upstream
print("Aggregating for upstream WTU, adding to table and writing as grid...")
Pannual_sum_WTU <- zonal(Pavg_annual_km3,WTU,fun='sum')
df[,5] <- Pannual_sum_WTU[,2]
df[,2] <- df[,5]/WTU_specs[,2]*1000000 #km3 to mm
P_WTU_grid <- subs(WTU,df,by=1,which=2)
writeRaster(P_WTU_grid,paste(outdir,"P_avg_annual_WT_mm.tif",sep=""),overwrite=T)
print("Done")

#zonal statistics over downstream
print("Aggregating for downstream WTU, adding to table and writing as grid...")
Pannual_sum_downstream <- zonal(Pavg_annual_km3,downstream,fun='sum')
colnames(Pannual_sum_downstream)<-c("WTU_ID","P")
temp <- merge(df, Pannual_sum_downstream, by = "WTU_ID",all.x=T)
df[,6] <- temp$P
df[,3] <-df[,6]/WTU_specs[,3]*1000000 #km3 to mm
P_downstream_grid <- subs(downstream,df,by=1,which=3)
writeRaster(P_downstream_grid,paste(outdir,"P_avg_annual_DS_mm.tif",sep=""),overwrite=T)
print("Done")

#zonal statistics over total basin
print("Aggregating for total basin, adding to table and writing as grid...")
Pannual_sum_basin <- zonal(Pavg_annual_km3,basins,fun='sum')
colnames(Pannual_sum_basin)<-c("WTU_ID","P")
temp <- merge(df, Pannual_sum_basin, by = "WTU_ID",all.x=T)
df[,7] <- temp$P
df[,4] <- df[,7]/WTU_specs[,4]*1000000 #km3 to mm
P_basin_grid <- subs(basins,df,by=1,which=4)
writeRaster(P_basin_grid,paste(outdir,"P_avg_annual_basin_mm.tif",sep=""),overwrite=T)
print("Done")

##CALCULATE INTERANNUAL VARIABILITY
print("Calculating interannual variability...")
#write global grid of interannual variability
print("Calculating and writing global interannual variability grid...")
Pannual_var <- 1-((max(Pannual2)-min(Pannual2))/max(Pannual2))
writeRaster(Pannual_var,paste(outdir,"P_var_interannual.tif",sep=""),overwrite=T)
print("Done")

#create 3-dimensional array, dim1: upstream,downstream,basin. dim2:years. dim3:WTUs
print("Calculating aggregated annual P per year...")
df_var_inter <- array(NA,dim=c(3,nlayers(Pannual2),nrow(df)))
#loop over years, calculate WTU aggregated P per year
for(i in 1:nlayers(Pannual2))
{
  P1 <- subset(Pannual2,i)
  P2 <- resample(P1,template,method="bilinear")
  #weighting for area
  P3 <- P2*area_km2
  #zonal statistics over upstream
  P_WTU <- zonal(P3,WTU,fun='sum')
  df_var_inter[1,i,] <- (P_WTU[,2] / WTU_specs[,2]) * 1000
  #zonal statistics over downstream
  P_downstream <- zonal(P3,downstream,fun='sum')
  colnames(P_downstream)<-c("WTU_ID","P")
  temp <- merge(df, P_downstream, by = "WTU_ID",all.x=T)
  df_var_inter[2,i,] <- (temp$P / WTU_specs[,3]) * 1000
  #zonal statistics over basin
  P_basin <- zonal(P3,basins,fun='sum')
  colnames(P_basin)<-c("WTU_ID","P")
  temp <- merge(df, P_basin, by = "WTU_ID",all.x=T)
  df_var_inter[3,i,] <- (temp$P / WTU_specs[,4]) * 1000
  print(paste("Values year ",i," aggregated...",sep=""))
}
print("Done")

#calculate interannual variability per unit
print("Calculate interannual variability per unit, add to table and write as grids...")
df_var_inter2 <- as.data.frame(matrix(ncol=3,nrow=nrow(df)))
for(i in 1:nrow(df))
{
  df_var_inter2[i,1] <- 1-((max(df_var_inter[1,,i])-min(df_var_inter[1,,i]))/max(df_var_inter[1,,i]))
  df_var_inter2[i,2] <- 1-((max(df_var_inter[2,,i])-min(df_var_inter[2,,i]))/max(df_var_inter[2,,i]))
  df_var_inter2[i,3] <- 1-((max(df_var_inter[3,,i])-min(df_var_inter[3,,i]))/max(df_var_inter[3,,i]))
}

#add to table
df[,8] <- df_var_inter2[,1]
df[,9] <- df_var_inter2[,2]
df[,10] <- df_var_inter2[,3]

#write as grids
Pannual_var_WTU_grid <- subs(WTU,df,by=1,which=8)
writeRaster(Pannual_var_WTU_grid,paste(outdir,"P_var_interannual_WT.tif",sep=""),overwrite=T)

Pannual_var_downstream_grid <- subs(downstream,df,by=1,which=9)
writeRaster(Pannual_var_downstream_grid,paste(outdir,"P_var_interannual_DS.tif",sep=""),overwrite=T)

Pannual_var_basin_grid <- subs(basins,df,by=1,which=10)
writeRaster(Pannual_var_basin_grid,paste(outdir,"P_var_interannual_basin.tif",sep=""),overwrite=T)
print("Done")

##CALCULATE INTRA-ANNUAL variability
print("Calculating intra-annual variability...")
#shift coordinates and convert from m to mm
Pymonmean_shift <- shift(Pymonmean,x=0.125)
P1 <- crop(Pymonmean_shift, c(xmn=0, xmx=180,ymn=-90,ymx=90))
P2 <- crop(Pymonmean_shift, c(xmn=180, xmx=360,ymn=-90,ymx=90))
P2_s <- shift(P2,x=-360)
Pmonthly2 <- merge(P1,P2_s)*1000

#write global grid
print("Calculating and writing global grid of intra-annual variability...")
Pmonthly_var <- 1-((max(Pmonthly2)-min(Pmonthly2))/max(Pmonthly2))
#Pmonthly_var_res <- resample(Pmonthly_var,template,method="ngb")
writeRaster(Pmonthly_var,paste(outdir,"P_var_intraannual.tif",sep=""),overwrite=T)
print("Done")

#create 3-dimensional array, dim1: upstream,downstream,basin. dim2:months. dim3:WTUs
print("Calculating aggregated average P per month...")
df_var_intra <- array(NA,dim=c(3,nlayers(Pmonthly2),nrow(df)))
#loop over years, calculate WTU aggregated P per year
for(i in 1:nlayers(Pmonthly2))
{
  P1 <- subset(Pmonthly2,i)
  P2 <- resample(P1,template,method="bilinear")
  #weighting for area
  P3 <- P2*area_km2
  #zonal statistics over upstream
  P_WTU <- zonal(P3,WTU,fun='sum')
  df_var_intra[1,i,] <- (P_WTU[,2] / WTU_specs[,2])
  #zonal statistics over downstream
  P_downstream <- zonal(P3,downstream,fun='sum')
  colnames(P_downstream)<-c("WTU_ID","P")
  temp <- merge(df, P_downstream, by = "WTU_ID",all.x=T)
  df_var_intra[2,i,] <- (temp$P / WTU_specs[,3])
  #zonal statistics over basin
  P_basin <- zonal(P3,basins,fun='sum')
  colnames(P_basin)<-c("WTU_ID","P")
  temp <- merge(df, P_basin, by = "WTU_ID",all.x=T)
  df_var_intra[3,i,] <- (temp$P / WTU_specs[,4])
  print(paste("Values month ",i," aggregated...",sep=""))
}
print("Done")

#calculate intraannual variability per unit
print("Calculate intraannual variability per unit, add to table and write as grids...")
df_var_intra2 <- as.data.frame(matrix(ncol=3,nrow=nrow(df)))
for(i in 1:nrow(df))
{
  df_var_intra2[i,1] <- 1-((max(df_var_intra[1,,i])-min(df_var_intra[1,,i]))/max(df_var_intra[1,,i]))
  df_var_intra2[i,2] <- 1-((max(df_var_intra[2,,i])-min(df_var_intra[2,,i]))/max(df_var_intra[2,,i]))
  df_var_intra2[i,3] <- 1-((max(df_var_intra[3,,i])-min(df_var_intra[3,,i]))/max(df_var_intra[3,,i]))
}

#add to table
df[,11] <- df_var_intra2[,1]
df[,12] <- df_var_intra2[,2]
df[,13] <- df_var_intra2[,3]

#write as grids
Pmonthly_var_WTU_grid <- subs(WTU,df,by=1,which=11)
writeRaster(Pmonthly_var_WTU_grid,paste(outdir,"P_var_intraannual_WT.tif",sep=""),overwrite=T)

Pmonthly_var_downstream_grid <- subs(downstream,df,by=1,which=12)
writeRaster(Pmonthly_var_downstream_grid,paste(outdir,"P_var_intraannual_DS.tif",sep=""),overwrite=T)

Pmonthly_var_basin_grid <- subs(basins,df,by=1,which=13)
writeRaster(Pmonthly_var_basin_grid,paste(outdir,"P_var_intraannual_basin.tif",sep=""),overwrite=T)
print("Done")

#write table
print("Writing CSV table with P indicators...")
write.csv(df,paste(outdir,"WTU_P_indicators.csv",sep=""))
print("Done")
