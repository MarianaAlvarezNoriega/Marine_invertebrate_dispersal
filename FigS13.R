
rm(list=ls())
dev.off()

library(dplyr)
library(rotl)
library(brms)
library(MCMCglmm)
library(ggplot2)
library(gridExtra)
library(ggstance)
library(broom)

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
datPD<-dat[,c(7,8,11:12,14)]
datPD<-datPD[complete.cases(datPD),]

# Exclude genera that are not in the matrix
common_genera<-intersect(factor(datPD$Genus), names1)
Matrix<-Matrix[which(names1 %in%  common_genera),which(names1 %in%  common_genera)]
datPD<-datPD[datPD$Genus %in% common_genera,]



# Fit model
prior1<-get_prior(log10(EggSize)~DevelopmentalMode:absLatitude:Hemisphere+absLatitude*DevelopmentalMode,
                  family=gaussian(),data=datPD)

model<-brm(log10(EggSize)~DevelopmentalMode:absLatitude:Hemisphere+absLatitude*DevelopmentalMode,
          data = datPD,seed=17,
          family = gaussian(), 
          control = list(max_treedepth = 13, adapt_delta=0.99), chains=3, prior=prior1,thin = 5, iter=40000,
          warmup=5000, file='EggSize_noPhylogeny')

# Get a table with summary statistics
tidy(model, conf_level=0.95, chain=F)[1:6,]






Egg_post<-posterior_samples(model)[1:6]
#write.table(Egg_post, 'egg_posterior_no phylo.csv', sep=',', row.n=F)

datPD<-dat[,c(7:12,14)]
datPD<-datPD[complete.cases(datPD),]
new<-with(datPD[datPD$DevelopmentalMode=='P',], expand.grid('Latitude'=seq(min(Latitude), max(Latitude), by=1), 
                                                  'Longitude'=seq(-180, 180, by=5), 'DevelopmentalMode'='P'))
new$absLatitude<-abs(new$Latitude)
new$Hemisphere<-ifelse(new$Latitude<0,'S','N')
pr<-as.data.frame(fitted(model, new, re_formula=NA))
new$fit<-pr$Estimate


world<-map_data(map='world')
world<-map_data("world") %>% 
  filter(region != "Antarctica")


P<-ggplot()+geom_raster(data=new, aes(x=Longitude, y=Latitude, fill=fit))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group))+theme_classic()+
  scale_fill_gradient( high='#990000', low='#FFCCCC',
                       breaks=c(log10(120),log10(110),log10(100),log10(90)),
                       labels=c(120,110,100,90))+
  geom_count(data=datPD[datPD$DevelopmentalMode=='P',], aes(x=Longitude,y=Latitude),  col='white', alpha=0.4, show.legend = F)+
  guides(fill=guide_legend(title='Egg size \n(\u03BCm)'))+
  theme(axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'),
        text=element_text(size=20), plot.title=element_text(hjust=0.5))+
  ggtitle('a- Planktotrophic')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  xlab('')+ylab('')
P



new<-with(datPD[datPD$DevelopmentalMode=='L',], expand.grid('Latitude'=seq(min(Latitude), max(Latitude), by=1), 
                                                  'Longitude'=seq(-180, 180, by=5),'DevelopmentalMode'='L'))
new$absLatitude<-abs(new$Latitude)
new$Hemisphere<-ifelse(new$Latitude<0,'S','N')
pr<-as.data.frame(fitted(model, new, re_formula=NA))
new$fit<-pr$Estimate




L<-ggplot()+geom_raster(data=new, aes(x=Longitude, y=Latitude, fill=fit))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group))+theme_classic()+
  scale_fill_gradient( high='#990000', low='#FFCCCC',breaks=c(log10(280),log10(260),log10(240), log10(220),log10(200)),
                       labels=c(280,260,240,220,200))+
  geom_count(data=datPD[datPD$DevelopmentalMode=='L',], aes(x=Longitude,y=Latitude),  col='white', alpha=0.4, show.legend = F)+
  guides(fill=guide_legend(title='Egg size \n(\u03BCm)'))+
  theme(axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'),
        text=element_text(size=20), plot.title=element_text(hjust=0.5))+
  ggtitle('b- Lecithotrophic')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  xlab('')+ylab('')
L

grid.arrange(P,L, ncol=2)
gA<-ggplotGrob(P)
gB<-ggplotGrob(L)

FigS13<-gtable_cbind(gA,gB)

plot(FigS13)

