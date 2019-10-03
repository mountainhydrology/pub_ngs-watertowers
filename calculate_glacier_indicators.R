##arthurlutz 20190129 / 20190724
##Script to preprocess grids for glacier indicators of water tower units for NGS-Water Tower Index project
rm(list=ls(all=TRUE))
library(raster)
library(maptools)
library(rgdal)
library(sf)

##SETTINGS
base <- "e:\\Dropbox (FutureWater)\\Team\\Projects\\Active\\2019001_NGS_WaterTowers\\data\\"
WTU <- raster(paste(base,"index\\units\\WTU.tif",sep=""))
WTU_specs <- read.csv(paste(base,"index\\units\\WTU_specs.csv",sep=""))
outdir <- paste(base,"Glaciers\\processed\\",sep="")
resolution <- 0.05

#glaciers input
glacier_area <- raster(paste(base,"Glaciers\\Farinotti_NGEO_2019\\p05_degree_glacier_area_km2.tif",sep=""))
glacier_volume <- raster(paste(base,"Glaciers\\Farinotti_NGEO_2019\\p05_degree_glacier_volume_km3.tif",sep=""))
#glacier mass balance (calculated for WTUs manually by Walter from WGMS data)
wtu_mb <- read_sf(paste(base,"Glaciers\\WGMS\\WTU_MB.shp",sep=""))

#precipitation input
P <- raster(paste(base,"ERA5\\processed\\P_avg_annual_mm.tif",sep=""))

##SETTINGS END

##increase speed of raster operations
rasterOptions(maxmemory = 1e+09)

#create raster template at project resolution
template <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,resolution=c(resolution,resolution))

###GLACIER INDICATORS
#initiate dataframe to store glacier values per WTU
no <- length(unique(WTU))
df<-as.data.frame(matrix(ncol=7,nrow=no))
df[,1]<-1:no
colnames(df)<-c("WTU_ID","Ice_Volume_WT_km3","Ice_Area_WT_km2","P_over_glacier_WT_mmyr-1","Glacier_MB_WT_mmyr-1","Meltwater_Yield_Glaciers_mmyr-1","Meltwater_Yield_WTU_km3yr-1")

#sum glacier volume per WTU
G_volume_WTU <- zonal(glacier_volume,WTU,fun='sum')
df[,2] <- G_volume_WTU[,2]
G_volume_WTU_grid <- subs(WTU,df,by=1,which=2)
writeRaster(G_volume_WTU_grid,paste(outdir,"Glac_volume_WT_km3.tif",sep=""),overwrite=T)

#sum glacier area per WTU
G_area_WTU <- zonal(glacier_area,WTU,fun='sum')
df[,3] <- G_area_WTU[,2]
G_area_WTU_grid <- subs(WTU,df,by=1,which=3)
writeRaster(G_area_WTU_grid,paste(outdir,"Glac_area_WT_km2.tif",sep=""),overwrite=T)

#mean P over glacier area (mm/yr)
df[,4] <- (zonal(P * glacier_area,WTU,fun="sum") / zonal(glacier_area,WTU,fun="sum"))[,2]

#add glacier MB
MB <- cbind(wtu_mb$WTU_ID,wtu_mb$MB)
colnames(MB)<-c("WTU_ID","MB")
temp <- merge(df, MB, by = "WTU_ID",all.x=T)
df[,5] <- temp$MB

#calculate annual meltwater yield
#for glacier area
df[,6] <- df$`P_over_glacier_WT_mmyr-1`-df$`Glacier_MB_WT_mmyr-1`
df[,6][is.na(df[,6])]<-0
#for WTu
df[,7] <- df$`Meltwater_Yield_Glaciers_mmyr-1` * zonal(glacier_area,WTU,fun="sum")[,2] * 0.000001

#write table
write.csv(df,paste(outdir,"WTU_Glacier_indicators.csv",sep=""))
