##arthurlutz 20190711
##Script to analyse the sensitivity of Water Tower Index ranking to weighting of indicators 
rm(list=ls(all=TRUE))
library(dplyr)
##SETTINGS
base = "d:/Dropbox (FutureWater)/Team/Projects/Active/2019001_ngs_watertowers/data/"
indicators = read.csv(paste(base, "Indicators/indicators.csv",sep=""))
n <- 10000 #number of index calculations
outfile <- paste(base,"Indicators/sensitivity_indicators_weights.csv",sep="")
min_weight <- 1
max_weight <- 3


##SETTINGS end

##get WTU ranking of equal weighted WTI and add as first column to dataframe
WTI_ranked <-indicators[order(-indicators$fin_sd),]
rank <- as.data.frame(c(1:nrow(WTI_ranked)))
colnames(rank)<-"rank"
WTI_ranked <- cbind(WTI_ranked,rank)
WTI_ranked <- WTI_ranked[order(WTI_ranked$WTU_ID),]
df_ranking <- as.data.frame(matrix(nrow = nrow(WTI_ranked),ncol=1+n))
df_ranking[,1] <- WTI_ranked$rank
orig_rank <- WTI_ranked$rank


df <- as.data.frame(matrix(ncol=5,nrow=n,1))
##weights of indicators up to max times larger, for 4 indicators
df[1:(n/8),1] <- runif(1:(n/8),min=min_weight,max=max_weight)
df[((n/8)+1):(2*n/8),2] <- runif(1:(n/8),min=min_weight,max=max_weight)
df[((2*n/8)+1):(3*n/8),3] <- runif(1:(n/8),min=min_weight,max=max_weight)
df[((3*n/8)+1):(4*n/8),4] <- runif(1:(n/8),min=min_weight,max=max_weight)
##weights of indicators up to 1/max times smaller, for 4 indicators
df[((4*n/8)+1):(5*n/8),1] <- runif(1:(n/8),min=1/max_weight,max=min_weight)
df[((5*n/8)+1):(6*n/8),2] <- runif(1:(n/8),min=1/max_weight,max=min_weight)
df[((6*n/8)+1):(7*n/8),3] <- runif(1:(n/8),min=1/max_weight,max=min_weight)
df[((7*n/8)+1):(8*n/8),4] <- runif(1:(n/8),min=1/max_weight,max=min_weight)
##sum of weights
df[,5]<- df[,1]+df[,2]+df[,3]+df[,4]

##randomize order of rows for supply indicators and demand indicators
weights_sup <- sample_n(df, size=nrow(df), replace = FALSE)
weights_dem <- sample_n(df, size=nrow(df), replace = FALSE)
  
for (i in 1:n)
{

#calculate WTI
  indicators$demtot <- ((indicators$dem_irr * weights_dem[i,1]) + (indicators$dem_ind * weights_dem[i,2]) + (indicators$dem_dom * weights_dem[i,3]) + (indicators$dem_nat * weights_dem[i,4]))/weights_dem[i,5]
  indicators$suptot <-  ((indicators$ptot * weights_sup[i,1]) + (indicators$stot * weights_sup[i,2]) + (indicators$gtot * weights_sup[i,3]) + (indicators$swtot * weights_sup[i,4]))/weights_sup[i,5]
  indicators$fin_sd <- indicators$demtot * indicators$suptot
  
  #rank by WTI and add ranking to dataframe
  WTI_ranked <-indicators[order(-indicators$fin_sd),]
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
WTI_deviation_classes <- as.data.frame(matrix(nrow=nrow(WTI_deviation),ncol=11))
for(i in 1:nrow(WTI_deviation))
{
  a <- table(as.vector(as.numeric(WTI_deviation[i,2:(n+1)])))
  WTI_deviation_classes[i,1] <- if(length(a[names(a)==0])>0){a[names(a)==0]} else 0
  WTI_deviation_classes[i,2] <- if(length(a[names(a)==1])>0){a[names(a)==1]} else 0
  WTI_deviation_classes[i,3] <- if(length(a[names(a)==2])>0){a[names(a)==2]} else 0
  WTI_deviation_classes[i,4] <- if(length(a[names(a)==3])>0){a[names(a)==3]} else 0
  WTI_deviation_classes[i,5] <- if(length(a[names(a)==4])>0){a[names(a)==4]} else 0
  WTI_deviation_classes[i,6] <- if(length(a[names(a)==5])>0){a[names(a)==5]} else 0
  WTI_deviation_classes[i,7] <- if(length(a[names(a)==6])>0){a[names(a)==6]} else 0
  WTI_deviation_classes[i,8] <- if(length(a[names(a)==7])>0){a[names(a)==7]} else 0
  WTI_deviation_classes[i,9] <- if(length(a[names(a)==8])>0){a[names(a)==8]} else 0
  WTI_deviation_classes[i,10] <- if(length(a[names(a)==9])>0){a[names(a)==9]} else 0
  WTI_deviation_classes[i,11] <- n - sum(WTI_deviation_classes[i,1:10])
}

#express deviation as percentage and replace NA by zeros
WTI_deviation_classes_perc <- (WTI_deviation_classes / n)*100
WTI_deviation_classes_perc[is.na(WTI_deviation_classes_perc)] <- 0
colnames(WTI_deviation_classes_perc)<- c(0,1,2,3,4,5,6,7,8,9,'>9')
WTI_deviation_classes[,12] <- sum(WTI_deviation_classes[,1:10],na.rm = T) 

#append WTU_ID, WTI_rank, and WTU name, and rank deviations
df <- cbind(indicators$WTU_ID,indicators$WTU_name,orig_rank,WTI_deviation_classes_perc)
df2 <-df[order(df$orig_rank),]
write.csv(df2,outfile)
