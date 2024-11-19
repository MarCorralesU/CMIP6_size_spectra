library(plyr)

# from table 3 (Corrales-Ugalde et al.)  mean size spectra coefficients from CMIP6 models and PSSdb (UVP and Scanner products)
source<- c('CESM', 'CMCC', 'CNRM', 'GFDL', 'GISS', 'IPSL', 'UKESM', 'PSSdb')
slope<- c(-1.03, -1.14, -1.18, -1.14, -1.61, -1.27, -1.22, -1.07)
intercept<- c(2.29e13, 4.43e13, 4.14e14, 1.29e15, 3.02e18, 4.22e14, 2.20e14, 3.88e13)

coeff_df<-data.frame(source, slope, intercept)
coeff_df$log10_intercept<- log10(coeff_df$intercept)

# get biovol (um3/m2), this dataset is after performing depth integration

d<-read.csv(path.expand('~/GIT/CMIP6_size_spectra/data/PSSdb_data_full.csv'), sep=',')


i=0

meanlogBiovol = mean(log10(d$total_biovolume), na.rm=TRUE)
sdlogBiovol =  sd(log10(d$total_biovolume), na.rm=TRUE)

set.seed(i)
randomLogBiovolume = rlnorm(n=nrow(d), meanlog=meanlogBiovol+i, sdlog=sdlogBiovol+i)
hist(randomLogBiovolume)
hist(log10(randomLogBiovolume))

hist(d$total_biovolume)
hist(log10(d$total_biovolume))
# biovolume_sum should be lognormal, roughly

# get NB using the model coefficients
source<-c()
biovols<-c()

for (s in coeff_df$source){
   set.seed(i)
   randomLogBiovolume <-rlnorm(n=nrow(d), meanlog=meanlogBiovol+i, sdlog=sdlogBiovol+i)
   logNB<-(coeff_df$slope[coeff_df$source==s]*log10(randomLogBiovolume))+coeff_df$log10_intercept[coeff_df$source==s] # use size spectra to get log normalized biovolume
   biovol<- (10^(logNB))*randomLogBiovolume # multiply it by the biovol 'size class' to get the total biovolume of each size
   source<-c(source, s)
   biovols<-c(biovols, sum(biovol)) # sum the biovolume across comparable size ranges 
}

coeff_df$modelled_sumBiovol<- biovols

write.csv(coeff_df, file.path(path.expand('~/GIT/CMIP6_size_spectra/data/Biovol_estimates_same_size_range.csv')), row.names=FALSE)


