#setwd('~/...')

# Load dispersal predictions
dat<-read.csv('mercator_predictions.csv', sep=',')


# Load correlations at different spatial scales
d1<-read.csv('correlation_10km.csv', sep=',')
d1$km<-log10(10)
d2<-read.csv('correlation_100km.csv', sep=',')
d2$km<-log10(100)
d3<-read.csv('correlation_1000km.csv', sep=',')
d3$km<-log10(1000)

d<-rbind(d1,d2,d3)


library(dplyr)
library(ggplot2)
library(ggstance)


# Get median and 95% quantiles
sum_d<-as.data.frame(d %>% group_by(km) %>% summarise('median'=median(x), 'lwr'=quantile(x, 0.025),
                                                           'upr'=quantile(x, 0.975)))



library(ggplot2)

FigS6<-ggplot()+
  geom_jitter(data=d, aes(x=km, y=x), height=0, pch=21, fill='grey', alpha=0.4)+
  geom_point(data=sum_d, aes(x=km, y=median, ymin=lwr, ymax=upr), size=4)+
  geom_errorbar(data=sum_d, aes(x=km, y=median, ymin=lwr, ymax=upr), width=0.1)+
  ylab("Pearson's correlation coefficient")+
  xlab("Maximum distance (km)")+
  theme_classic()+
  scale_x_continuous(breaks=c(log10(0.1),log10(1),log10(10), log10(100), log10(1000)),
                     labels=c(0.1,1,10,100,1000))+
  theme(text=element_text(size=18))+
  scale_y_continuous(breaks=c(-0.4,-0.2,0,0.2,0.4,0.6,0.8,1))+
  geom_rect(aes(xmin=min(dat$dispersal,na.rm=T), xmax=max(dat$dispersal,na.rm=T), ymin=1.05, ymax=1.15), fill='#CCCCCC', col='white')+
  geom_boxploth(data=dat, aes(x=dispersal, y=1.1), size=0.5,width=0.05,alpha=0.5,outlier.size=1, outlier.alpha=0.01, fill='white',
                col='black')

FigS6
