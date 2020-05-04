closeAllConnections() 
rm(list=ls())
dev.off()

library(dplyr)
library(rotl)
library(brms)
library(MCMCglmm)
library(ggplot2)
library(gridExtra)
library(ggstance)
library(viridis)

# Set working directory
#setwd(~/...)

dat<-read.csv('Marine_invertebrate_data.csv', sep=',')

tree<-read.tree('MarineInvertebrateTree.tre')


# Get the matrix for the random effects
inv.matrix<-MCMCglmm::inverseA(tree, nodes='TIPS', scale=T)
Matrix<-solve(inv.matrix$Ainv)

# Match genera from the matrix to genera in the data
changenames<-as.vector(rownames(inv.matrix$Ainv))

names1<-c(); otts<-c()
for (i in 1:length(changenames)){
  names1[i]<-unlist(strsplit(as.character(changenames[i]), '_'))[1]
  n<-length(unlist(strsplit(as.character(changenames[i]), '_')))
  otts[i]<-unlist(strsplit(as.character(changenames[i]), '_'))[n]
}

rownames(Matrix)<-names1
colnames(Matrix)<-names1



# Select only the data required
datPD<-dat[,c(7,8,11:14)]
datPD<-datPD[complete.cases(datPD),]

# Exclude genera that is not in the matrix
common_genera<-intersect(factor(datPD$Genus), names1)
Matrix<-Matrix[which(names1 %in%  common_genera),which(names1 %in%  common_genera)]
datPD<-datPD[datPD$Genus %in% common_genera,]
datPD<-subset(datPD, !datPD$PlanktonicTime==0)
datPD$LogPlank<-log10(datPD$PlanktonicTime)




prior1<-get_prior(LogPlank~Hemisphere:absLatitude+absLatitude*DevelopmentalMode+
                    log10(EggSize)*DevelopmentalMode+(1|Genus),
                  family=gaussian(),data=datPD)

model<-brm(LogPlank~Hemisphere:absLatitude+absLatitude*DevelopmentalMode+
             log10(EggSize)*DevelopmentalMode+(1|Genus),
           data = datPD,
           family = gaussian(), cov_ranef = list(Matrix),
           control = list(max_treedepth = 13, adapt_delta=0.99), chains=3, prior=prior1,thin = 5, iter=30000,
           warmup=5000, file='PlankDuration')


datPD<-dat[,c(7,8,11:14,16)]
datPD<-datPD[complete.cases(datPD),]
mod_temp<-brm(absLatitude~s(temperature_annual), data=datPD, control=list(adapt_delta=0.99), chains=3,
              thin = 5, iter=20000, warmup=10000, file='Temperature')


datn<-with(datPD[datPD$DevelopmentalMode=='L',], expand.grid('EggSizeLog'=seq((min(log10(EggSize))),(max(log10(EggSize))), l=100),
                                                    'temperature_annual'=seq(min(temperature_annual), max(temperature_annual), by=1), 
                                                    'DevelopmentalMode'=c('L')))
datn$EggSize<-10^datn$EggSizeLog
datn$absLatitude<-as.data.frame(fitted(mod_temp, datn))$Estimate
datn$Hemisphere<-'N'

pr<-as.data.frame(fitted(model, datn, re_formula=NA))
datn$fit<-10^pr$Estimate




library(viridis)

L<-ggplot()+
  geom_tile(data=datn[datn$Hemisphere=='N',], aes(x=log10(EggSize), y=temperature_annual, fill=fit))+
  scale_fill_viridis(name='PD (days)')+
  geom_point(data=datPD[datPD$DevelopmentalMode=='L',], aes(x=log10(EggSize), y=temperature_annual, size=PlanktonicTime), 
             col='#CCCCCC', alpha=0.4)+
  theme_classic()+
  scale_x_continuous(breaks=c(log10(100), log10(200),log10(300), log10(400),log10(500), log10(600), log10(700), log10(800), log10(900),
                              log10(1000)),
                     labels= c(expression(10^2),'','','','','','','','',expression(10^3)),
                     limits=c(log10(45), log10(1200)))+
  scale_y_continuous(breaks=c(5,10,15,20,25,30), labels=c('','10','','20', '', '30'))+
  ylab('Temperature (annual mean),°C')+
  xlab('Egg size, \u03BCm')+
  theme(text = element_text(size=20), plot.title=element_text(hjust=0.5))+
  theme( axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'))+
  ggtitle('b- Non-feeding larvae')+labs(size='Observed')+
  guides( fill = guide_legend(order = 1),count = guide_colorbar(order = 0))

L



datn<-with(datPD[datPD$DevelopmentalMode=='P',], expand.grid('EggSizeLog'=seq((min(log10(EggSize))),(max(log10(EggSize))), l=100),
                                                             'temperature_annual'=seq(min(temperature_annual), max(temperature_annual), by=1), 
                                                             'DevelopmentalMode'=c('P')))
datn$EggSize<-10^datn$EggSizeLog
datn$absLatitude<-as.data.frame(fitted(mod_temp, datn))$Estimate
datn$Hemisphere<-'N'

pr<-as.data.frame(fitted(model, datn, re_formula=NA))
datn$fit<-10^pr$Estimate







P<-ggplot()+
  geom_tile(data=datn[datn$Hemisphere=='N',], aes(x=log10(EggSize), y=temperature_annual, fill=fit))+
  scale_fill_viridis(name='PD (days)')+
  geom_point(data=datPD[datPD$DevelopmentalMode=='P',], aes(x=log10(EggSize), y=temperature_annual, size=PlanktonicTime), 
             col='#CCCCCC', alpha=0.4)+
  theme_classic()+
  scale_x_continuous(breaks=c(log10(100), log10(200),log10(300), log10(400),log10(500), log10(600), log10(700), log10(800), log10(900),
                              log10(1000)),
                     labels= c(expression(10^2),'','','','','','','','',expression(10^3)),
                     limits=c(log10(45), log10(1200)))+
  scale_y_continuous(breaks=c(5,10,15,20,25,30), labels=c('','10','','20', '', '30'))+
  ylab('Temperature (annual mean), °C')+
  xlab('Egg size, \u03BCm')+
  theme(text = element_text(size=20), plot.title=element_text(hjust=0.5))+
  theme( axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'))+
  ggtitle('a- Feeding larvae')+labs(size='Observed')+
  guides( fill = guide_legend(order = 1),count = guide_colorbar(order = 0))

P


gA<-ggplotGrob(P)
gB<-ggplotGrob(L+ylab(''))


#Fig 2
#library(Cairo)

plot(gtable_cbind(gA,gB))
#ggsave("Fig2.pdf",
#       plot=plot(gtable_cbind(gA,gB)),
#        device=cairo_pdf,  dpi=300)


