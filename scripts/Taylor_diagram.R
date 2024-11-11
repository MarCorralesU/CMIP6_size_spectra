#library(plotrix)
#trace(taylor.diagram, edit=TRUE)
#taylor.diagram_ed <- edit(taylor.diagram) #or fix
#environment(taylor.diagram_ed) <- asNamespace('plotrix')
#assignInNamespace("taylor.diagram", taylor.diagram_ed, ns = "plotrix")
setwd(file.path(path.expand('~/GIT/CMIP6_size_spectra/scripts')))
source('Taylor_custom_function.R')
library(dplyr)
library(tidyverse)
library(doBy)
library(Metrics)


ref = read.csv(file.path(path.expand('~/GIT/CMIP6_size_spectra/scripts/PSSdb_data_Taylor.csv')), as.is=TRUE)
ref <-ref[order(ref$biomes,ref$month),]
ref <- ref[!(ref$biomes=='HCSS' & ( ref$month==12)),]#remove winter data for HCSS since it has very few data points
ref$log_biovolume<-log10(ref$total_biovolume)


d = read.csv(file.path(path.expand('~/GIT/CMIP6_size_spectra/scripts/CMIP6_data_Taylor.csv')), as.is=TRUE)
d <-d[order(d$source,d$experiment,d$biomes,d$month),]
d <- d[!(d$biomes=='HCSS' & (d$month==12)),]
d$log_biovolume<-log10(d$total_biovolume)
d$intercept<-log10(d$intercept)
d <- subset(d, d$experiment=='hist')

d_full = rbind(ref,d)


# section to get tables with correlation coefficients

# get corr coefficient
biomes<-c()
models<-c()
R2_slope<-c()
R2_intercept<-c()
R2_biovol<-c()

d_PSSdb<- d_full[d_full$source=='PSSdb',]


for (b in unique(d_full$biomes)){
  d_f<- d_full[d_full$biomes==b,]
  for (m in unique(d_f$source[d_f$experiment=='hist'])){
    biomes<-c(biomes, b)
    
    
    r_slope<-cor(d_f$slope[d_f$source=='PSSdb'], d_f$slope[d_f$source==m], use = "pairwise")
    R_slope<-c(R_slope, r_slope)
    
    
    r_intercept<-corr(d_f$intercept[d_f$source=='PSSdb'], d_f$intercept[d_f$source==m], use = "pairwise")
    R_intercept<-c(R_intercept, r_intercept)
    
    
    r_biovol<-corr(d_f$log_biovolume[d_f$source=='PSSdb'], d_f$log_biovolume[d_f$source==m], use = "pairwise")
    R_biovol<-c(R_biovol, r_biovol)
    
    models<-c(models, m)
  }
}
R_df<-data.frame(models, biomes, R_slope, R_intercept, R_biovol)
write.csv(R_df, file.path(path.expand('~/GIT/CMIP6_size_spectra/data/corr_CMIP6.csv')), row.names=FALSE)

# get standard deviation
biomes<-c()
models<-c()
SD_slope<-c()
SD_intercept<-c()
SD_biovol<-c()

d_PSSdb<- d_full[d_full$source=='PSSdb',]


for (b in unique(d_full$biomes)){
  d_f<- d_full[d_full$biomes==b,]
  for (m in unique(d_f$source)){
    biomes<-c(biomes, b)
    
    
    sd_slope<-sd(d_f$slope[d_f$source==m])
    SD_slope<-c(SD_slope, sd_slope)
    
    
    sd_intercept<-sd(d_f$intercept[d_f$source==m])
    SD_intercept<-c(SD_intercept, sd_intercept)
    
    
    sd_biovol<-sd(d_f$log_biovolume[d_f$source==m])
    SD_biovol<-c(SD_biovol, sd_biovol)
    
    models<-c(models, m)
  }
}
SD_df<-data.frame(models, biomes, SD_slope, SD_intercept, SD_biovol)
write.csv(SD_df, file.path(path.expand('~/GIT/CMIP6_size_spectra/data/stDev_CMIP6.csv')), row.names=FALSE)

# get RMSE
biomes<-c()
models<-c()
RMSE_slope<-c()
RMSE_intercept<-c()
RMSE_biovol<-c()

d_PSSdb<- d_full[d_full$source=='PSSdb',]


for (b in unique(d_full$biomes)){
  d_f<- d_full[d_full$biomes==b,]
  for (m in unique(d_f$source[d_f$experiment=='hist'])){
    biomes<-c(biomes, b)
  

    rmse_slope<-rmse(d_f$slope[d_f$source=='PSSdb'], d_f$slope[d_f$source==m])
    RMSE_slope<-c(RMSE_slope, rmse_slope)
  
   
    rmse_intercept<-rmse(d_f$intercept[d_f$source=='PSSdb'], d_f$intercept[d_f$source==m])
    RMSE_intercept<-c(RMSE_intercept, rmse_intercept)
  
    
    rmse_biovol<-rmse(d_f$log_biovolume[d_f$source=='PSSdb'], d_f$log_biovolume[d_f$source==m])
    RMSE_biovol<-c(RMSE_biovol, rmse_biovol)
  
    models<-c(models, m)
}
}
RMSE_df<-data.frame(models, biomes, RMSE_slope, RMSE_intercept, RMSE_biovol)
write.csv(RMSE_df, file.path(path.expand('~/GIT/CMIP6_size_spectra/data/RMSE_CMIP6.csv')), row.names=FALSE)


#Centered RMSE
centered_rmse <- function(x, y){
  xprime = x - mean(x, na.rm=TRUE)
  yprime = y - mean(y, na.rm=TRUE)
  return(sqrt(sum((xprime-yprime)^2)/length(complete.cases(x))))
}


biomes<-c()
models<-c()
CRMSE_slope<-c()
CRMSE_intercept<-c()
CRMSE_biovol<-c()



for (b in unique(d_full$biomes)){
  d_f<- d_full[d_full$biomes==b,]
  for (m in unique(d_f$source[d_f$experiment=='hist'])){
    biomes<-c(biomes, b)
    
    #crmse_slope<-sqrt(((sqrt(mean((d_f$slope[d_f$source=='PSSdb']-d_f$slope[d_f$source==m])^2)))^2) - (((mean(d_f$slope[d_f$source=='PSSdb']))-(mean(d_f$slope[d_f$source==m])))^2))
    crmse_slope<-centered_rmse(d_f$slope[d_f$source=='PSSdb'], d_f$slope[d_f$source==m])
    CRMSE_slope<-c(CRMSE_slope, crmse_slope)
    
    #crmse_intercept<-sqrt(((sqrt(mean((d_f$intercept[d_f$source=='PSSdb']-d_f$intercept[d_f$source==m])^2)))^2) - (((mean(d_f$intercept[d_f$source=='PSSdb']))-(mean(d_f$intercept[d_f$source==m])))^2))
    crmse_intercept<-centered_rmse(d_f$intercept[d_f$source=='PSSdb'], d_f$intercept[d_f$source==m])
    CRMSE_intercept<-c(CRMSE_intercept, crmse_intercept)
    
    #crmse_biovol<-sqrt(((sqrt(mean((d_f$log_biovolume[d_f$source=='PSSdb']-d_f$log_biovolume[d_f$source==m])^2)))^2) - (((mean(d_f$log_biovolume[d_f$source=='PSSdb']))-(mean(d_f$log_biovolume[d_f$source==m])))^2))
    crmse_biovol<-centered_rmse(d_f$log_biovolume[d_f$source=='PSSdb'], d_f$log_biovolume[d_f$source==m])
    CRMSE_biovol<-c(CRMSE_biovol, crmse_biovol)
    
    models<-c(models, m)
  }
}
CRMSE_df<-data.frame(models, biomes, CRMSE_slope, CRMSE_intercept, CRMSE_biovol)
write.csv(CRMSE_df, file.path(path.expand('~/GIT/CMIP6_size_spectra/data/CenteredRMSE_CMIP6.csv')), row.names=FALSE)

# modified 11/11/2024, separate taylor diagram by biome and coefficient

###.       SLOPE
for (biom in c('LC', 'HCSS', 'HCPS')){
# for this, remember how you had to use fix(taylor.diagram) to change the maxsd parameter
  if (biom=='LC'){
    symbol<-17
  }
  else if(biom=='HCSS'){
    symbol<-16
  }
  else if(biom=='HCPS'){
    symbol<-15
  }
   pdf(file = paste("/Users/mc4214/GIT/CMIP6_size_spectra/figures/Taylor_slopes_mod_SCALE_", biom,".pdf", sep=""),   # The directory you want to save the file in
       width = 7, # The width of the plot in inches
       height = 7)
   taylor.diagram_mod(ref$slope[ref$biomes==biom], d$slope[d$source=='CESM' & d$biomes==biom & d$experiment=='hist'],col="#0077bb", pch=symbol,ref.sd=FALSE,
                        grad.corr.lines=c(0.2,0.4,0.6,0.8, 0.9), pos.cor = FALSE, 
                        #pcex=1.5,cex.axis=0.1,
                        mar=c(5,5,5,5), gamma.col=19, 
                        pcex=2, cex.axis=1, cex.lab=1, alpha=0.5,
                        lwd=10,font=3,lty=3, x_axis_range=1)
   taylor.diagram_mod(ref$slope[ref$biomes==biom], d$slope[d$source=='CMCC' & d$biomes==biom & d$experiment=='hist'], add =TRUE, col='#33BBEE', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
   taylor.diagram_mod(ref$slope[ref$biomes==biom], d$slope[d$source=='CNRM' & d$biomes==biom & d$experiment=='hist'], add =TRUE, col='#009988', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
   taylor.diagram_mod(ref$slope[ref$biomes==biom], d$slope[d$source=='GFDL' & d$biomes==biom & d$experiment=='hist'], add =TRUE,col='#EE7733', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
   taylor.diagram_mod(ref$slope[ref$biomes==biom], d$slope[d$source=='GISS' & d$biomes==biom & d$experiment=='hist'], add =TRUE, col='#CC3311', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
   taylor.diagram_mod(ref$slope[ref$biomes==biom], d$slope[d$source=='IPSL' & d$biomes==biom & d$experiment=='hist'], add =TRUE,col='darkgoldenrod1', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
   taylor.diagram_mod(ref$slope[ref$biomes==biom], d$slope[d$source=='UKESM' & d$biomes==biom & d$experiment=='hist'],add =TRUE, col='#BBBBBB', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)

   lpos<-sd(ref$slope)*0.6# remember this controls x position of legend
   legend(x=lpos*0.95, y=lpos*1.1, cex= 0.5, legend=c('LC', 'HCSS', 'HCPS'),pch=c(17, 16, 15))
   legend(x=lpos*0.7, y=lpos*1.1, legend=c("CESM","CMCC","CNRM","GFDL","GISS","IPSL","UKESM"),pch=16,cex = 0.5, col=c("#0077BB","#33BBEE","#009988","#EE7733","#CC3311",'darkgoldenrod1','#BBBBBB'))

   dev.off()

}

###.       INTERCEPT
for (biom in c('LC', 'HCSS', 'HCPS')){
  # for this, remember how you had to use fix(taylor.diagram) to change the maxsd parameter
  if (biom=='LC'){
    symbol<-17
  }
  else if(biom=='HCSS'){
    symbol<-16
  }
  else if(biom=='HCPS'){
    symbol<-15
  }
  pdf(file = paste("/Users/mc4214/GIT/CMIP6_size_spectra/figures/Taylor_intercept_mod_experiment_zoomedbitmore_", biom,".pdf", sep=""),   # The directory you want to save the file in
      width = 7, # The width of the plot in inches
      height = 7)
  taylor.diagram_mod(ref$intercept[ref$biomes==biom], d$intercept[d$source=='CESM' & d$biomes==biom & d$experiment=='hist'],col="#0077bb", pch=symbol,ref.sd=FALSE,
                     grad.corr.lines=c(0.2,0.4,0.6,0.8, 0.9), pos.cor = FALSE, 
                     #pcex=1.5,cex.axis=0.1,
                     mar=c(5,5,5,5), gamma.col=19, 
                     pcex=2, cex.axis=1, cex.lab=1, alpha=0.5,
                     lwd=10,font=3,lty=3, x_axis_range=0.25)
  taylor.diagram_mod(ref$intercept[ref$biomes==biom], d$intercept[d$source=='CMCC' & d$biomes==biom & d$experiment=='hist'], add =TRUE, col='#33BBEE', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$intercept[ref$biomes==biom], d$intercept[d$source=='CNRM' & d$biomes==biom & d$experiment=='hist'], add =TRUE, col='#009988', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$intercept[ref$biomes==biom], d$intercept[d$source=='GFDL' & d$biomes==biom & d$experiment=='hist'], add =TRUE,col='#EE7733', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$intercept[ref$biomes==biom], d$intercept[d$source=='GISS' & d$biomes==biom & d$experiment=='hist'], add =TRUE, col='#CC3311', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$intercept[ref$biomes==biom], d$intercept[d$source=='IPSL' & d$biomes==biom & d$experiment=='hist'], add =TRUE,col='darkgoldenrod1', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$intercept[ref$biomes==biom], d$intercept[d$source=='UKESM' & d$biomes==biom & d$experiment=='hist'],add =TRUE, col='#BBBBBB', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  
  lpos<-sd(ref$intercept)*2# remember this controls x position of legend
  #legend(x=lpos*0.95, y=lpos*1.1, cex= 0.5, legend=c('LC', 'HCSS', 'HCPS'),pch=c(17, 16, 15))
  legend(x=lpos*0.7, y=lpos*1.1, legend=c("CESM","CMCC","CNRM","GFDL","GISS","IPSL","UKESM"),pch=16,cex = 0.5, col=c("#0077BB","#33BBEE","#009988","#EE7733","#CC3311",'darkgoldenrod1','#BBBBBB'))
  
  dev.off()
  
}

###.       biovolume
for (biom in c('LC', 'HCSS', 'HCPS')){
  # for this, remember how you had to use fix(taylor.diagram) to change the maxsd parameter
  if (biom=='LC'){
    symbol<-17
  }
  else if(biom=='HCSS'){
    symbol<-16
  }
  else if(biom=='HCPS'){
    symbol<-15
  }
  pdf(file = paste("/Users/mc4214/GIT/CMIP6_size_spectra/figures/Taylor_biovolume_mod_experiment_superzoomed", biom,".pdf", sep=""),   # The directory you want to save the file in
      width = 7, # The width of the plot in inches
      height = 7)
  taylor.diagram_mod(ref$log_biovolume[ref$biomes==biom], d$log_biovolume[d$source=='CESM' & d$biomes==biom & d$experiment=='hist'],col="#0077bb", pch=symbol,ref.sd=FALSE,
                     grad.corr.lines=c(0.2,0.4,0.6,0.8, 0.9), pos.cor = FALSE, 
                     #pcex=1.5,cex.axis=0.1,
                     mar=c(5,5,5,5), gamma.col=19, 
                     pcex=2, cex.axis=1, cex.lab=1, alpha=0.5,
                     lwd=10,font=3,lty=3, x_axis_range=0.3)
  taylor.diagram_mod(ref$log_biovolume[ref$biomes==biom], d$log_biovolume[d$source=='CMCC' & d$biomes==biom & d$experiment=='hist'], add =TRUE, col='#33BBEE', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$log_biovolume[ref$biomes==biom], d$log_biovolume[d$source=='CNRM' & d$biomes==biom & d$experiment=='hist'], add =TRUE, col='#009988', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$log_biovolume[ref$biomes==biom], d$log_biovolume[d$source=='GFDL' & d$biomes==biom & d$experiment=='hist'], add =TRUE,col='#EE7733', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$log_biovolume[ref$biomes==biom], d$log_biovolume[d$source=='GISS' & d$biomes==biom & d$experiment=='hist'], add =TRUE, col='#CC3311', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$log_biovolume[ref$biomes==biom], d$log_biovolume[d$source=='IPSL' & d$biomes==biom & d$experiment=='hist'], add =TRUE,col='darkgoldenrod1', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  taylor.diagram_mod(ref$log_biovolume[ref$biomes==biom], d$log_biovolume[d$source=='UKESM' & d$biomes==biom & d$experiment=='hist'],add =TRUE, col='#BBBBBB', pcex=2,alpha=0.5,pch=symbol, x_axis_range=0.5)
  
  lpos<-sd(ref$log_biovolume)*2# remember this controls x position of legend
  #legend(x=lpos*0.95, y=lpos*1.1, cex= 0.5, legend=c('LC', 'HCSS', 'HCPS'),pch=c(17, 16, 15))
  legend(x=lpos*0.7, y=lpos*1.1, legend=c("CESM","CMCC","CNRM","GFDL","GISS","IPSL","UKESM"),pch=16,cex = 0.5, col=c("#0077BB","#33BBEE","#009988","#EE7733","#CC3311",'darkgoldenrod1','#BBBBBB'))
  
  dev.off()
  
}


                                  ## explore pssdb file climatology, to determine whether data removal in HCSS is valid
library(tidyverse)
library (dplyr)
df_UVP = read.csv('/Users/mc4214/GIT/PSSdb/raw/NBSS_data/NBSS_ver_10_2023/UVP_1b_Size-spectra-fit_v2023-10.csv', as.is=TRUE)
df_Scanner = read.csv('/Users/mc4214/GIT/PSSdb/raw/NBSS_data/NBSS_ver_10_2023/Scanner_1b_Size-spectra-fit_v2023-10.csv', as.is=TRUE)
df_PSSdb <-rbind(df_UVP, df_Scanner)
df_basins <- df_PSSdb[,c('latitude', 'longitude', 'ocean')] %>% rename(lat = latitude, lon=longitude)
df_PSSdb_biome <- read.csv('/Users/mc4214/Documents/CMIP6_PSS_paper/data/PSSdb_data_full.csv', as.is=TRUE)
df_PSSdb_biome = na.omit(merge(df_PSSdb_biome, df_basins, by=c('lat', 'lon'), all=TRUE)) %>% distinct()


plot_slope <- ggplot(df_PSSdb_biome, aes(x=month, y=slope))+   # fill colours the ribbons and other shapes
                    geom_point(aes(colour=ocean))+
                    stat_summary(geom='ribbon', fun.data = mean_sdl, fun.args=list(mult=1),alpha = 0.5)+   # added transparency
                    stat_summary(geom='line', fun.y = mean, size=1)+
                    stat_summary(fun.data = "mean_cl_normal", colour = "#009988", size = 0.5, linewidth=0)+
                    #ylim(-1.4,-0.8)+
                    scale_x_continuous(breaks = c(2,4,6,8,10,12), labels = seq(2, 12, by=2))+
  
                    xlab('Month')+
                    ylab(expression(slope~(L^-1~mu~m^-3)))+
                    theme(axis.text=element_text(size=12),
                          axis.title=element_text(size=14,face="bold"),
                            # Hide panel borders and remove grid lines
                            panel.border = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            # Remove panel background
                            panel.background = element_blank(),
                            # Change axis line
                            axis.line = element_line(colour = "#009988"))+
  facet_wrap(~factor(biomes, levels=c('LC', 'HCSS', 'HCPS')))

   
plot_slope

plot_intercept <- ggplot(df_PSSdb_biome, aes(x=month, y=intercept))+   # fill colours the ribbons and other shapes
  geom_point(aes(colour=ocean))+
  stat_summary(geom='ribbon', fun.data = mean_sdl, fun.args=list(mult=1),alpha = 0.5)+   # added transparency
  stat_summary(geom='line', fun.y = mean, size=1)+
  stat_summary(fun.data = "mean_cl_normal", colour = "#009988", size = 0.5, linewidth=0)+
  #ylim(-1.4,-0.8)+
  scale_x_continuous(breaks = c(2,4,6,8,10,12), labels = seq(2, 12, by=2))+
  
  xlab('Month')+
  ylab(expression(intercept~(mu~m^-3~L^-1~mu~m^-3)))+ #Intercept ($\mu$m$^{-3}$ m$^{-2}$ $\mu$m$^{-3}$)
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "#009988"))+
  facet_wrap(~factor(biomes, levels=c('LC', 'HCSS', 'HCPS')))


plot_intercept

                                ## explore pssdb file climatology, to determine whether data removal in HCSS is valid
df_PSSdb_subset  <- df_PSSdb_biome[df_PSSdb_biome$biomes =='LC', ]
df_PSSdb_subset_2  <- df_PSSdb_subset[df_PSSdb_subset$ocean %in% c('South Pacific Ocean', 'South Atlantic Ocean', 'North Pacific Ocean', 'North Atlantic Ocean'), ]

plot_slope <- ggplot(df_PSSdb_subset_2, aes(x=month, y=slope))+   # fill colours the ribbons and other shapes
  geom_point()+
  stat_summary(geom='ribbon', fun.data = mean_sdl, fun.args=list(mult=1),alpha = 0.5)+   # added transparency
  stat_summary(geom='line', fun.y = mean, size=1)+
  stat_summary(fun.data = "mean_cl_normal", colour = "#009988", size = 0.5, linewidth=0)+
  ylim(-2, -0.7)+
  scale_x_continuous(breaks = c(2,4,6,8,10,12), labels = seq(2, 12, by=2))+
  
  xlab('Month')+
  ylab(expression(slope~(L^-1~mu~m^-3)))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "#009988"),
        strip.text = element_text(size = 14))+
  facet_wrap(~factor(ocean, levels=c('North Pacific Ocean', 'North Atlantic Ocean','South Pacific Ocean', 'South Atlantic Ocean')), ncol=2)


plot_slope

plot_slope_all_basins <- ggplot(df_PSSdb_subset, aes(x=month, y=slope))+   # fill colours the ribbons and other shapes
  geom_point()+
  stat_summary(geom='ribbon', fun.data = mean_sdl, fun.args=list(mult=1),alpha = 0.5)+   # added transparency
  stat_summary(geom='line', fun.y = mean, size=1)+
  stat_summary(fun.data = "mean_cl_normal", colour = "#009988", size = 0.5, linewidth=0)+
  ylim(-2, -0.7)+
  scale_x_continuous(breaks = c(2,4,6,8,10,12), labels = seq(2, 12, by=2))+
  
  xlab('Month')+
  ylab(expression(slope~(L^-1~mu~m^-3)))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "#009988"),
        strip.text = element_text(size = 14))+
  facet_wrap(~factor(ocean))


plot_slope_all_basins


plot_slope_LC<- ggplot(df_PSSdb_subset, aes(x=month, y=slope))+   # fill colours the ribbons and other shapes
  geom_point()+
  stat_summary(geom='ribbon', fun.data = mean_sdl, fun.args=list(mult=1),alpha = 0.5)+   # added transparency
  stat_summary(geom='line', fun.y = mean, size=1)+
  stat_summary(fun.data = "mean_cl_normal", colour = "#009988", size = 0.5, linewidth=0)+
  ylim(-2, -0.7)+
  scale_x_continuous(breaks = c(2,4,6,8,10,12), labels = seq(2, 12, by=2))+
  
  xlab('Month')+
  ylab(expression(slope~(L^-1~mu~m^-3)))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "#009988"),
        strip.text = element_text(size = 14))
plot_slope_LC
