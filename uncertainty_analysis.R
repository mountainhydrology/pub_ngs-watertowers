##arthur lutz 20190711
##Script to include input uncertainty in WT Index calculation 
rm(list=ls())
library(sf)

##SETTINGS
# Specific base folder
base = "d:/Dropbox (FutureWater)/Team/Projects/Active/2019001_ngs_watertowers/data/"

#local folder with water gap runs output
base_local <- "d:/workdir/"

#file of final indicators used to extract ranking as in manuscript
indicators = read.csv(paste(base, "Indicators/indicators.csv",sep=""))

# specify indicator files of subindicators as input for uncertainty analysis
pindf = paste(base, "ERA5/processed/WTU_P_indicators.csv",sep="")
sindf = paste(base,"Snow/processed/WTU_Snow_indicators.csv",sep="")
gindf = paste(base,"Glaciers/processed/WTU_Glacier_indicators.csv",sep="")
swindf = paste(base,"HydroLakes/processed/WTU_lake_storage_volume.csv",sep="")

# outputfiles
outfile <- paste(base,"Indicators/uncertainty/uncertainty_analysis_ranking.csv",sep="")
outfile_scoring_sd <- paste(base,"Indicators/uncertainty/uncertainty_analysis_scoring.csv",sep="")

# specify file with WTU specs
wtuf = paste(base,"index/units/WTU_specs.csv",sep="")

##uncertainty ranges input
p_wtu_u <- read.csv(paste(base,"uncertainty/precipitation/P_uncertainty_per_WTU.csv",sep=""))
s_u <- 0.05 #sd for snow persistence multiplication (limit resulting snow persistence at 0 and 1)
i_v_u <- read.csv(paste(base,"uncertainty/ice_volume/WTU_IceVol_uncertainty.csv",sep=""))
i_mb_u <- as.data.frame(read_sf(paste(base,"Glaciers/WGMS/WTU_MB.shp",sep="")))
l_u <- c(0.1) #sd for lake and reservoir volume multiplication

#number of repetitions for WTI calculation
n <- 1000
##SETTINGS END

##get WTU ranking and scoring (normalized and non-normalized) of equal weighted WTI and add as first column to dataframe
WTI_ranked <-indicators[order(-indicators$fin_sd),]
rank <- as.data.frame(c(1:nrow(WTI_ranked)))
colnames(rank)<-"rank"
WTI_ranked <- cbind(WTI_ranked,rank)
WTI_ranked <- WTI_ranked[order(WTI_ranked$WTU_ID),]
df_ranking <- as.data.frame(matrix(nrow = nrow(WTI_ranked),ncol=1+n))
df_ranking[,1] <- WTI_ranked$rank
orig_rank <- WTI_ranked$rank

df_scoring <- as.data.frame(matrix(nrow = nrow(WTI_ranked),ncol=1+n))
df_scoring[,1]<-indicators$fin_sd_nor

df_scoring_nn <- as.data.frame(matrix(nrow = nrow(WTI_ranked),ncol=1+n))
df_scoring_nn[,1]<-indicators$fin_sd

## Read CSV with base indicators
pind <- read.csv(file=pindf, header=TRUE, sep=",")
str(pind)
sind <- read.csv(file=sindf, header=TRUE, sep=",")
str(sind)
gind <- read.csv(file=gindf, header=TRUE, sep=",")
str(gind)
swind <- read.csv(file=swindf, header=TRUE, sep=",")
str(swind)

#read CSV with WTU specs

wtu <- read.csv(file=wtuf, header=TRUE, sep=",")
str(wtu)

#loop over number of repetitions
for (i in 1:n)
{
  #create indicator df
  indicator <- subset(wtu,select = c(WTU_ID,WTU_name))
  
  # precipitation
  #read P downstream multiplication factor
  p_ds_mult <- read.csv(paste(base_local,"uncertainty\\water_gap_runs\\",sprintf("%05d",i),"_P_downstream_multiplication_factor.csv",sep=""))
  p_ds_mult2 <- merge(indicator,p_ds_mult,by="WTU_ID",all.x=T)
  p_ds_mult2 <- replace(p_ds_mult2,is.na(p_ds_mult2),1)
  p_wtu_mult <- rnorm(nrow(wtu),1,min(1,p_wtu_u$sd_normalized))
  p_wtu <- pind$Pannual_WT_km3 * p_wtu_mult
  pind <- replace(pind,is.na(pind),0)
  p_ds <- pind$Pannual_basin_downstream_km3 * p_ds_mult2$P_downstream_mult
  p_basin <- p_wtu + p_ds
  indicator$pup_tot <- p_wtu / p_basin
  indicator$pinter <- pind$Var_inter_upstream
  indicator$pintra <- pind$Var_intra_upstream
  # final indicator
  # first average pinter and pintra, then multiply with pup_tot
  indicator$ptot <- (0.5*(indicator$pintra+indicator$pinter))*indicator$pup_tot
  
  # snow
  indicator$snow_p <- sapply(sind$Mean_annual_snow_persistence_WT * rnorm(nrow(wtu),1,s_u),function(y) min(max(y,0),1))
  indicator$snow_intra <- sind$Var_intraannual_snow_persistence_WT
  indicator$snow_inter <- sind$Var_interannual_snow_persistence_WT
  # final indicator
  # first average snow_inter and snow_intra, then multiply with snow_p
  indicator$stot <- (0.5*(indicator$snow_intra+indicator$snow_inter))*indicator$snow_p
  
  # glaciers
  ice_vol <- gind$Ice_Volume_WT_km3 * rnorm(nrow(wtu),1,i_v_u$IceVol.SD.normalized)
  indicator$glac_v <- ice_vol / (p_wtu + ice_vol)
 
  #melt water yield indicator
  Pglac <- gind$P_over_glacier_WT_mmyr.1 * p_wtu_mult
  MB <- cbind(i_mb_u$WTU_ID,i_mb_u$MB,i_mb_u$MB_UNC)
  colnames(MB)<-c("WTU_ID","MB","MB_UNC")
  temp <- merge(wtu, MB, by = "WTU_ID",all.x=T)
  
  mb_sd_normalized <- abs(temp$MB_UNC/temp$MB)
  glac_mb <- rnorm(nrow(wtu),gind$Glacier_MB_WT_mmyr.1,i_mb_u$MB_UNC)
  glac_meltyield <- (Pglac-glac_mb) * gind$Ice_Area_WT_km2*0.000001
  
  indicator$glac_m <- glac_meltyield/(glac_meltyield+p_wtu)
  #set negative and NA to zero
  indicator$glac_m[indicator$glac_m<0]<-0
  indicator$glac_m[is.na(indicator$glac_m)]<-0
  
  #final indicator
  indicator$gtot <- (indicator$glac_v+indicator$glac_m)/2
  
  # surface water
  sw_mult <- rnorm(nrow(wtu),1,l_u)
  indicator$swtot <- (swind$Lake_Storage_Volume_WT_km3 * sw_mult) / (p_wtu  + (swind$Lake_Storage_Volume_WT_km3 * sw_mult))
  indicator$swtot[is.na(indicator$swtot)] <- 0
  
  # demand indicators
  dind <- read.csv(paste(base_local,"uncertainty/water_gap_runs/",sprintf("%05d",i),"_WTU_Demand_indicators.csv",sep=""))
  dind_orig <- read.csv(paste(base,"WaterDemandsWRI\\processed\\WTU_Demand_indicators.csv",sep=""))
  dirrind <- read.csv(paste(base_local,"uncertainty/water_gap_runs/",sprintf("%05d",i),"_WTU_Irrigation_Water_Gap_monthly.csv",sep=""))
  dindind <- read.csv(paste(base_local,"uncertainty/water_gap_runs/",sprintf("%05d",i),"_WTU_Industrial_Water_Gap_monthly.csv",sep=""))
  ddomind <- read.csv(paste(base_local,"uncertainty/water_gap_runs/",sprintf("%05d",i),"_WTU_Domestic_Water_Gap_monthly.csv",sep=""))
  dnatind <- read.csv(paste(base_local,"uncertainty/water_gap_runs/",sprintf("%05d",i),"_WTU_Natural_Water_Gap_monthly.csv",sep=""))
  irr_wg <- dirrind$Water_gap_total
  indicator$dem_irr <- irr_wg / dind$Irr_demand_basin_km3yr.1
  indicator$dem_irr[is.na(indicator$dem_irr)] <- 0
  
  ind_wg <- dindind$Water_gap_total
  indicator$dem_ind <-  ind_wg / dind$Ind_demand_basin_km3yr.1
  indicator$dem_ind[is.na(indicator$dem_ind)] <- 0
  
  dom_wg <- ddomind$Water_gap_total
  indicator$dem_dom <- dom_wg / dind$Dom_demand_basin_km3yr.1
  indicator$dem_dom[is.na(indicator$dem_dom)] <- 0
  
  nat_wg <- dnatind$Water_gap_total
  indicator$dem_nat <- nat_wg / dind$Nat_demand_basin_km3yr.1
  indicator$dem_nat[is.na(indicator$dem_nat)] <- 0
  
  indicator$demtot <- (indicator$dem_irr + indicator$dem_ind + indicator$dem_dom + indicator$dem_nat)/4
  indicator$suptot <-  (indicator$ptot + indicator$stot+indicator$gtot + indicator$swtot)/4
  indicator$fin_sd <- indicator$demtot * indicator$suptot
  indicator$fin_sd_nor <- (indicator$fin_sd-min(indicator$fin_sd))/(max(indicator$fin_sd)-min(indicator$fin_sd))
  
  #add fin_sd_norm to scoring dfs
  df_scoring[,i+1] <- indicator$fin_sd_nor
  df_scoring_nn[,i+1] <- indicator$fin_sd
  
  #rank by WTI and add ranking to dataframe
  WTI_ranked <-indicator[order(-indicator$fin_sd),]
  rank <- as.data.frame(c(1:nrow(WTI_ranked)))
  colnames(rank)<-"rank"
  WTI_ranked <- cbind(WTI_ranked,rank)
  WTI_ranked <- WTI_ranked[order(WTI_ranked$WTU_ID),]
  df_ranking[,1+i] <- WTI_ranked$rank
}
  
  #calculate position difference in ranking
  WTI_deviation <- df_ranking
  for (i in 1:n)
  {
    WTI_deviation[,(1+i)] <- abs(WTI_deviation[,1]-df_ranking[,(1+i)])
  }
  
  #count number of deviations per class
  WTI_deviation_classes <- as.data.frame(matrix(nrow=nrow(WTI_deviation),ncol=6))
  for(i in 1:nrow(WTI_deviation))
  {
    a <- table(as.vector(as.numeric(WTI_deviation[i,2:(n+1)])))
    WTI_deviation_classes[i,1] <- if(length(a[names(a)==0])>0){a[names(a)==0]} else 0
    WTI_deviation_classes[i,2] <- if(length(a[names(a)==1])>0){a[names(a)==1]} else 0
    WTI_deviation_classes[i,3] <- if(length(a[names(a)==2])>0){a[names(a)==2]} else 0
    WTI_deviation_classes[i,4] <- if(length(a[names(a)==3])>0){a[names(a)==3]} else 0
    WTI_deviation_classes[i,5] <- if(length(a[names(a)==4])>0){a[names(a)==4]} else 0
    WTI_deviation_classes[i,6] <- if(length(a[names(a)==5])>0){sum(a[names(a)==5])} else 0
    WTI_deviation_classes[i,7] <- if(length(a[names(a)==6])>0){sum(a[names(a)==6])} else 0
    WTI_deviation_classes[i,8] <- if(length(a[names(a)==7])>0){sum(a[names(a)==7])} else 0
    WTI_deviation_classes[i,9] <- if(length(a[names(a)==8])>0){sum(a[names(a)==8])} else 0
    WTI_deviation_classes[i,10] <- if(length(a[names(a)==9])>0){sum(a[names(a)==9])} else 0
    WTI_deviation_classes[i,11] <- n - sum(WTI_deviation_classes[i,1:10])
    
  }
  
  
  #express deviation as percentage and replace NA by zeros
  WTI_deviation_classes_perc <- (WTI_deviation_classes / n)*100
  WTI_deviation_classes_perc[is.na(WTI_deviation_classes_perc)] <- 0
  colnames(WTI_deviation_classes_perc)<- c(0,1,2,3,4,5,6,7,8,9,'>=9')
  WTI_deviation_classes[,12] <- sum(WTI_deviation_classes[,1:11],na.rm = T) 
  
  #append WTU_ID, WTI_rank, and WTU name, and rank deviations
  df <- cbind(indicators$WTU_ID,indicators$WTU_name,orig_rank,WTI_deviation_classes_perc)
  df2 <-df[order(df$orig_rank),]
  write.csv(df2,outfile)
  
  #calculate SD of scoring and uncertainty
  df <- as.data.frame(matrix(nrow=nrow(indicators),ncol=6))
  df[,1]<-indicators$WTU_ID
  df[,2]<- indicators$fin_sd
  df[,4]<- indicators$fin_sd_nor
  for (y in 1:nrow(df))
  {
    df[y,3] <- sd(df_scoring_nn[y,2:n+1])
    df[y,5] <- sd(df_scoring[y,2:n+1])
    df[y,6] <- paste(format(round(df[y,4],digits=2),nsmall=2)," ± ",format(round(df[y,5],digits=2),nsmall=2),sep="")
  }
  colnames(df)<-c("WTU_ID","WTI_unnorm_final","WTI_unnorm_sd","WTI_norm_final","WTI_norm_sd","merged")
  write.csv(df,outfile_scoring_sd)





