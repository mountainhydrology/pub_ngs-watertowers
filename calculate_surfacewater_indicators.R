##arthurlutz 20190129 / 20190725
##Script to preprocess grids for surface water indicators of water tower units for NGS-Water Tower Index project
rm(list=ls(all=TRUE))
library(raster)
library(maptools)
library(rgdal)
library(sf)
library(tibble)
library(rgeos)
library(lwgeom)

##SETTINGS
base <- "e:\\Dropbox (FutureWater)\\Team\\Projects\\Active\\2019001_NGS_WaterTowers\\data\\"
WTU <- raster(paste,base,"index\\units\\WTU.tif",sep=""))
WTU_shp <- read_sf(paste,base,"index\\units\\WTU_vector.shp",sep=""))
WTU_specs <- read.csv(paste(base,"index\\units\\wTU_specs.csv",sep=""))
outdir <- paste(base,"HydroLakes\\processed\\",sep="")
resolution <- 0.05

#hydrolakes input
hydrolakes_shp <- read_sf(paste(base,"HydroLakes\\HydroLAKES_polys_v10_shp\\HydroLAKES_polys_v10.shp",sep=""))

##SETTINGS END

###SURFACE WATER INDICATORS
#initiate dataframe to store surface values per WTU
no <- length(unique(WTU))
df<-as.data.frame(matrix(ncol=2,nrow=no))
df[,1]<-1:no
colnames(df)<-c("WTU_ID","Lake_Storage_Volume_WT_km3")

##make polygon geometry valid
print("Preprocessing polygon files for valid geometry...")
WTU_shp1 <- st_make_valid(WTU_shp)

hydrolakes_shp2 <- hydrolakes_shp
hydrolakes_shp1 <- st_make_valid(hydrolakes_shp2)

rm(hydrolakes_shp)
rm(hydrolakes_shp2)

print("Done")

#loop over WTUs
print("Looping over WTUs to extract lake volumes...")

for(i in 1:no)
{
  #select WTU feature
  WTU_shp2 <- WTU_shp1[WTU_shp1$WTU_ID == i,]
  
  #select contained features from hydrolakes
  lakes_WTU <- st_join(WTU_shp2,hydrolakes_shp1,join=st_contains)
  lakes_WTU_df <- as.data.frame(lakes_WTU)
    
  #sum volumes of lakes and convert to km3 (source values are million m3)
  df[i,2] <- sum(lakes_WTU_df[,13])*0.001
  print(paste("WTU ",i," of ",no," done"))
}

#write table and grid
print("Write table and grid...")
surface_water_grid <- subs(WTU,df,by=1,which=2)
writeRaster(surface_water_grid,paste(outdir,"WTU_surface_water_storage_km3.tif",sep=""),overwrite=T)
write.csv(df,paste(outdir,"WTU_lake_storage_volume.csv",sep=""))
print("Done")
