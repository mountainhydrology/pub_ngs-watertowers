##arthurlutz 20190129
##Script to preprocess grids for snow cover indicators of water tower units for NGS-Water Tower Index project
rm(list=ls(all=TRUE))
library(raster)
library(maptools)
library(rgdal)

##SETTINGS
base <- "e:\\Dropbox (FutureWater)\\Team\\Projects\\Active\\2019001_NGS_WaterTowers\\data\\"
WTU <- raster(paste(base,"index\\units\\WTU.tif",sep=""))
WTU_specs <- read.csv(paste(base,"index\\units\\WTU_specs.csv",sep=""),colClasses=c("NULL",NA,NA,NA,NA))
outdir <- paste(base,"Snow\\processed\\",sep="")
resolution <- 0.05

#snow input
Sannual <- brick(paste(base,"Snow\\MOD10CM006_yearmean_2001-2017.tif",sep=""))
Symonmean <- brick(paste(base,"Snow\\MOD10CM006_ymonmean_2001-2017.tif",sep=""))

##SETTINGS END

##increase speed of raster operations
rasterOptions(maxmemory = 1e+09)

#create raster template at project resolution
template <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,resolution=c(resolution,resolution))

#create raster with area per cell (km2)
area_km2 <- area(template)

#turn off scientific notation
options(scipen=999)

###SNOW INDICATORS
#initiate dataframe to store snow values per WTU
no <- length(unique(WTU))
df<-as.data.frame(matrix(ncol=4,nrow=no))
df[,1]<-1:no
colnames(df)<-c("WTU_ID","Mean_annual_snow_persistence_WT","Var_interannual_snow_persistence_WT","Var_intraannual_snow_persistence_WT")

##CALCULATE AVERAGE ANNUAL SNOW PERSISTENCE
print("Calculate and write global grid of average annual snow persistence...")
S_avgannual <- mean(Sannual)/100
#write grid
writeRaster(S_avgannual,paste(outdir,"Snow_persistence_avg_annual.tif",sep=""),overwrite=T)
print("Done")

#zonal statistics over WTs
print("Aggregating for upstream WTU, adding to table and writing as grid...")
S_avgannual_areacor <- S_avgannual * area_km2
S_avgannual_WTU <- zonal(S_avgannual_areacor,WTU,fun='sum')
df[,2] <- S_avgannual_WTU[,2] / WTU_specs[,2]
S_WTU_grid <- subs(WTU,df,by=1,which=2)
writeRaster(S_WTU_grid,paste(outdir,"Snow_persistence_avg_annual_WT.tif",sep=""),overwrite=T)
rm(S_avgannual_areacor)
rm(S_WTU_grid)
rm(S_avgannual)
print("Done")

##CALCULATE INTERANNUAL VARIABILITY
print("Calculating interannual variability...")
print("Calculating and writing global grid of interannual variability...")
Sannual_var <- 1-((max(Sannual/100)-min(Sannual/100))/max(Sannual/100))
writeRaster(Sannual_var,paste(outdir,"Snow_persistence_var_interannual.tif",sep=""),overwrite=T)
print("Done")

#create dataframe, dim1:years. dim2:WTUs
print("Calculating aggregated annual snow persistence per year...")
df_var_inter <- as.data.frame(matrix(ncol=nlayers(Sannual),nrow=nrow(df)))
#loop over years, calculate WTU aggregated snow persistence per year
for(i in 1:nlayers(Sannual))
{
  S1 <- subset(Sannual,i)/100
  #weighting for area
  S2 <- S1*area_km2
  #zonal statistics over WTU
  S_WTU <- zonal(S2,WTU,fun='sum')
  df_var_inter[,i] <- (S_WTU[,2] / WTU_specs[,2])
  print(paste("Values year ",i," aggregated...",sep=""))
}
print("Done")

#calculate interannual variability per unit
print("Calculate interannual variability per unit, add to table and write as grids...")
df_var_inter2 <- as.data.frame(matrix(ncol=1,nrow=nrow(df)))
for(i in 1:nrow(df))
{
  df_var_inter2[i,1] <- 1-((max(df_var_inter[i,])-min(df_var_inter[i,]))/max(df_var_inter[i,]))
}

#add to table
df[,3] <- df_var_inter2[,1]

#write as grids
Sannual_var_WTU_grid <- subs(WTU,df,by=1,which=3)
writeRaster(Sannual_var_WTU_grid,paste(outdir,"Snow_persistence_var_interannual_WT.tif",sep=""),overwrite=T)
print("Done")

##CALCULATE INTRA-ANNUAL VARIABILITY
print("Calculating intraannual variability...")
print("Calculating and writing global grid of intraannual variability...")
Sintra_var <- 1-((max(Symonmean/100)-min(Symonmean/100))/max(Symonmean/100))
writeRaster(Sintra_var,paste(outdir,"Snow_persistence_var_intraannual.tif",sep=""),overwrite=T)
print("Done")

#create dataframe, dim1:months. dim2:WTUs
print("Calculating aggregated annual snow persistence per month...")
df_var_intra <- as.data.frame(matrix(ncol=nlayers(Symonmean),nrow=nrow(df)))
#loop over years, calculate WTU aggregated snow persistence per month
for(i in 1:nlayers(Symonmean))
{
  S1 <- subset(Symonmean,i)/100
  #weighting for area
  S2 <- S1*area_km2
  #zonal statistics over WTU
  S_WTU <- zonal(S2,WTU,fun='sum')
  df_var_intra[,i] <- (S_WTU[,2] / WTU_specs[,2])
  print(paste("Values month ",i," aggregated...",sep=""))
}
print("Done")

#calculate intraannual variability per unit
print("Calculate intraannual variability per unit, add to table and write as grids...")
df_var_intra2 <- as.data.frame(matrix(ncol=1,nrow=nrow(df)))
for(i in 1:nrow(df))
{
  df_var_intra2[i,1] <- 1-((max(df_var_intra[i,])-min(df_var_intra[i,]))/max(df_var_intra[i,]))
}

#add to table
df[,4] <- df_var_intra2[,1]

#write as grids
Sintra_var_WTU_grid <- subs(WTU,df,by=1,which=4)
writeRaster(Sintra_var_WTU_grid,paste(outdir,"Snow_persistence_var_intraannual_WT.tif",sep=""),overwrite=T)
print("Done")

#write table
print("Writing CSV table with Snow indicators...")
write.csv(df,paste(outdir,"WTU_Snow_indicators.csv",sep=""))
print("Done")
