##Walter Immerzeel / Arthur Lutz 20190718
##Script to delineate water tower units for NGS-Water Tower Index project
##Assigns ice volume, area and snow persistence to GMBA polygons and subsets based on thresholds
##Updated for revision to delineate only downstream basin parts connected to WTU
##Result is saved as shapefile
rm(list=ls(all=TRUE))
library(sf)
library(raster)
library(maptools)
library(rgdal)
library(rgeos)
library(fasterize)
library(tmaptools)
library(stars)
library(spex)
library(tidyverse)
library(lwgeom)

##SETTINGS
base = "d:\\Dropbox (FutureWater)\\Team\\Projects\\Active\\2019001_NGS_WaterTowers\\data\\"
outdir <- "d:\\Dropbox (FutureWater)\\Team\\Projects\\Active\\2019001_NGS_WaterTowers\\data\\index\\units_v3\\"
hydrobasins <- st_make_valid(read_sf(paste(base,"Basins\\FAO\\major_hydrobasins.shp",sep="")))
continents <- c("europe","austpacific","centralam","neareast","northam","asia","southam")
resolution <- 0.05
##SETTINGS END

# Read the GMBA shapefile
gmbafn = paste(base,"Delineation/GMBA_MountainRanges/GMBA Mountain Inventory_v1.2-World.shp",sep="") 
gmba_ori <- read_sf(gmbafn)

# Filenames glacier ice volume and area and snow persistence rasters
glacvfn = paste(base,"Glaciers/Farinotti_NGEO_2019/p05_degree_glacier_volume_km3.tif",sep="")
glacafn = paste(base,"Glaciers/Farinotti_NGEO_2019/p05_degree_glacier_area_km2.tif",sep="")
snowfn =  paste(base,"snow/processed/Snow_persistence_avg_annual.tif",sep="")

# Summarize zones for glacier volume, glacier area and snow persistence
# Create a subset based on thresholds and save as new polygon

gmbafn2 = paste(base,"index/units/gmba_all.shp",sep="") 
gmba <- ZonalPipe(gmbafn, glacvfn, stat="sum")
names(gmba)[4]<-"vol_km3"
writeOGR(gmba, gmbafn2,"gmba_new","ESRI Shapefile",overwrite_layer=TRUE)

# Glacier area
gmba <- ZonalPipe(gmbafn2, glacafn, stat="sum")
names(gmba)[5]<-"area_km2"
writeOGR(gmba, gmbafn2,"gmba_all","ESRI Shapefile",overwrite_layer=TRUE)

# Snow persistence
gmba <- ZonalPipe(gmbafn2, snowfn, stat="mean")
names(gmba)[6]<-"snow_p"
writeOGR(gmba, gmbafn2,"gmba_all","ESRI Shapefile",overwrite_layer=TRUE)

# Thesholds for ice volume, glacier area and snow persistence
gv = 0.1
ga = 0.1
sp = 0.1

gmba2 <- subset(gmba,((vol_km3 > gv)|(snow_p > sp))&(area_km2 > ga))
gmba2
gmbafn3 = paste(base,"index/units/gmba_ss.shp",sep="") 
writeOGR(gmba2, gmbafn3,"gmba_ss","ESRI Shapefile",overwrite_layer=TRUE)

# intersect hydrobasins and filtered GMBA mountain ranges, then dissolve boundaries to get WTUs
mntns <- read_sf(gmbafn3)
WTU1 <- st_intersection(hydrobasins,mntns)
WTU2 <- st_make_valid(WTU1)
WTU3 <- WTU2 %>% group_by(WTU2$MAJ_BAS) %>% summarize_all(mean)

write_sf(WTU3,paste(outdir,"wtu_temp.shp",sep=""))

#read WTU shp
WTU <- st_make_valid(read_sf(paste(outdir,"wtu_temp.shp",sep="")))

#extract downstream areas of WTUs
#loop over continents using continent-based FAO hydrobasins shapefiles with subbasins
for (con in continents)
{
  #get subbasins for continent
  subbasins <- st_make_valid(read_sf(paste(base,"Basins\\FAO\\hydrobasins_",con,"\\hydrobasins_",con,".shp",sep="")))
  
  #subset of all subbasins that overlap with WTUs, and store their IDs
  subbasins_overlap_wtu <- st_intersection(WTU,subbasins)
  df_subbas <- unique(subbasins_overlap_wtu$SUB_BAS)
  
  #get IDs of downstream subbasins, downstream of subset
  downstream <- unique(subbasins_overlap_wtu$TO_BAS)
  
  #subset the subbasins with the downstream IDs, and store their IDs
  subbasins_x <- subset(subbasins,subbasins$SUB_BAS %in% downstream)
  df_subbas <- c(df_subbas,unique(subbasins_x$SUB_BAS))
  
  #repeat the above until no more downstream subbasins remain
  i <- 1
  repeat {
    downstream <- unique(subbasins_x$TO_BAS)
    subbasins_y <- subset(subbasins,subbasins$SUB_BAS %in% downstream)
    df_subbas <- c(df_subbas,unique(subbasins_x$SUB_BAS))
    subbasins_x <- subbasins_y
    print(paste(con,i,sep=" "))
    i <- i+1
  
    if(nrow(subbasins_x)<=1)
    {
      break
    }
    ##in neareast there is a sequence of 6 subbasins draining into each other resulting in endless loop. made exception below to exit the loop earlier
    if(con == "neareast" & nrow(subbasins_x)<=7)
    {
      break
    }
  }
  #extract polygons with WTU basins per continent
  basins_final <- unique(df_subbas)
  wtu_connected <- subset(subbasins,subbasins$SUB_BAS %in% basins_final)
  WTU_ds <- wtu_connected %>% group_by(wtu_connected$MAJ_BAS) %>% summarize_all(mean)
  
  #write as shapefile
  write_sf(WTU_ds,paste(outdir,"WTU_basins_",con,"_temp.shp",sep=""))
}

#create raster template at project resolution
template <- raster(xmn=-180,xmx=180,ymn=-90,ymx=90,resolution=c(resolution,resolution))

#rasterize continent-wise basins and merge rasters
hydrobasins_ras <- raster(template)
for (con in continents)
{
  grid <- fasterize(read_sf(paste(outdir,"WTU_basins_",con,"_temp.shp",sep="")), template, field="MAJ_BAS")
  hydrobasins_ras <- merge(hydrobasins_ras,grid)
}

#rasterize WTUs
WTU_grid <- fasterize(WTU, template, field="MAJ_BAS")

#get unique IDs of water towers and remove basins without water towers
ids <- unique(WTU_grid)
hydrobasins_wtu <- hydrobasins_ras
hydrobasins_wtu[!hydrobasins_wtu %in% ids] <- NA

#grid with downstream areas
downstream <- hydrobasins_wtu
downstream[WTU_grid>0]<-NA

#generate new IDs for water tower systems and get WTU name
no <- length(unique(WTU_grid))
df<-as.data.frame(matrix(ncol=3,nrow=no))
colnames(df)<-c("MAJ_BAS","WTU_ID","MAJ_NAME")
df[,1]<-unique(WTU_grid)
df[,2]<-1:no
hydrobasins_df <- as.data.frame(hydrobasins[,1:2])
hydrobasins_df2 <- hydrobasins_df[,1:2]
hydrobasins_df3 <- unique(hydrobasins_df2)
temp <- merge(df,hydrobasins_df3, by="MAJ_BAS",all.x=T)
df[,3] <- temp[,4]
WTU_names <- df

WTU_reclass <- subs(WTU_grid,df,by=1,which=2)
hydrobasins_wtu_reclass <- subs(hydrobasins_wtu,df,by=1,which=2)
downstream_reclass <- subs(downstream,df,by=1,which=2)

#write grids
writeRaster(as.integer(hydrobasins_wtu_reclass),paste(outdir,"basins.tif",sep=""),overwrite=T)
writeRaster(as.integer(WTU_reclass),paste(outdir,"WTU.tif",sep=""),overwrite=T)
writeRaster(as.integer(downstream_reclass),paste(outdir, "basins_downstream.tif",sep=""),overwrite=T)

##read grids if already written before
#hydrobasins_wtu_reclass <- raster(paste(outdir,"basins.tif",sep=""))
#WTU_reclass <- raster(paste(outdir,"WTU.tif",sep=""))
#downstream_reclass <- raster(paste(outdir, "basins_downstream.tif",sep=""))

#store WTU properties in table
no <- length(unique(WTU_reclass))
df<-as.data.frame(matrix(ncol=7,nrow=no))
df[,1]<-1:no
colnames(df)<-c("WTU_ID","Area_WT_km2","Area_downstream_km2","Area_basin_km2","WTU_name","WTU_continent","WTU_name_short")
area_km2 <- area(WTU_reclass)
area_WT <- zonal(area_km2,WTU_reclass,fun='sum')
df[,2] <- area_WT[,2]
area_downstream <- zonal(area_km2,downstream_reclass,fun='sum')
colnames(area_downstream)<-c("WTU_ID","A")
temp <- merge(df, area_downstream, by = "WTU_ID",all.x=T)
df[,3] <- temp$A
area_basin <- zonal(area_km2,hydrobasins_wtu_reclass,fun='sum')
colnames(area_basin)<-c("WTU_ID","A")
temp <- merge(df, area_basin, by = "WTU_ID",all.x=T)
df[,4] <- temp$A
df[,5] <- WTU_names[,3]
basin_maj_code <- zonal(hydrobasins_ras,WTU_grid,fun='mean')
df[,6] <- substr(basin_maj_code[,2],0,1)
#add manual exceptions and change to categorical
df[,6][df[,6]==8]<-5
df[,6][df[,6]==6]<-5
df[,6][df[,6]==2]<-1
df[,6][df[,6]==1]<-"A"
df[,6][df[,6]==3]<-"B"
df[,6][df[,6]==4]<-"C"
df[,6][df[,6]==5]<-"D"

##add abbreviated WTU names
df[,7] <- df[,5]
for (i in 1:nrow(df))
{
  df[i,7] <- str_replace(df[i,7], "Northwestern ", "NW-")
  df[i,7] <- str_replace(df[i,7], "Northwest ", "NW-")
  df[i,7] <- str_replace(df[i,7], "North ", "N-")
  df[i,7] <- str_replace(df[i,7], " - ", "-")
  df[i,7] <- str_replace(df[i,7], "South ", "S-")
  df[i,7] <- str_replace(df[i,7], "East ", "E-")
  df[i,7] <- str_replace(df[i,7], "West ", "W-")
  df[i,7] <- str_replace(df[i,7], "United States", "US")
  df[i,7] <- str_replace(df[i,7], "Pacific", "Pac.")
  df[i,7] <- str_replace(df[i,7], "Atlantic", "Atl.")
  df[i,7] <- str_replace(df[i,7], "Arctic", "Arc.")
  df[i,7] <- str_replace(df[i,7], "Arctic", "Arc.")
  df[i,7] <- str_replace(df[i,7], "S-Argentina, South Atl. Coast", "S-Arg., S-Atl. Coast")
  df[i,7] <- str_replace(df[i,7], "Colombia-Ecuador, Pac. Coast", "Col.-Ecuad., Pac. Coast")
  df[i,7] <- str_replace(df[i,7], "Adriatic Sea-Greece - Black Sea Coast", "Adr. Sea-Bl. Sea Coast")
  df[i,7] <- str_replace(df[i,7], "Interior", "Int.")
  df[i,7] <- str_replace(df[i,7], "Plateau of Tibet Int.", "Tibetan Plateau")
  df[i,7] <- str_replace(df[i,7], "Huang He", "Yellow River")
  df[i,5] <- str_replace(df[i,5], "Plateau of Tibet Interior", "Tibetan Plateau")
  df[i,5] <- str_replace(df[i,5], "Huang He", "Yellow River")
  df[i,5] <- str_replace(df[i,5], "Northwest Territories", "Northwest Territories and Nunavut")
  df[i,7] <- str_replace(df[i,7], "NW-Territories", "NW-Territ. and Nunavut")
}

write.csv(df,paste(outdir,"WTU_specs.csv",sep=""))

##create vector version of WTU grid and write as shapefile with WTU_ID and WTU_NAME attributes
WTU_vector <- rasterToPolygons(WTU_reclass,na.rm=T,dissolve=T)
WTU_names2 <- WTU_names[,2:3]
colnames(WTU_names2)<-c("WTU_ID","WTU_NAME")
WTU_vector2 <- append_data(WTU_vector, WTU_names2, key.shp = "WTU_ID", key.data = "WTU_ID")
rgdal::writeOGR(WTU_vector2,dsn=paste(outdir,"WTU_vector.shp",sep=""),layer="WTU_vector", driver="ESRI Shapefile", overwrite_layer=T)

##create vector version of basins grid and write as shapefile with WTU_ID and WTU_NAME attributes
basins_vector <- spex::polygonize(hydrobasins_wtu_reclass)
#dissolve polygons
basins_vector2 <- basins_vector %>% group_by(basins_vector$WTU_ID) %>% summarize_all(mean)
#write as shapefile
st_write(basins_vector2, paste(outdir,"basins_vector.shp",sep=""),delete_dsn=T)

##create vector version of downstream basins grid and write as shapefile with WTU_ID and WTU_NAME attributes
downstream_vector <- spex::polygonize(downstream_reclass)
#dissolve polygons
downstream_vector2 <- downstream_vector %>% group_by(downstream_vector$WTU_ID) %>% summarize_all(mean)

#write as shapefile
st_write(downstream_vector2, paste(outdir,"downstream_vector.shp",sep=""),delete_dsn=T)
