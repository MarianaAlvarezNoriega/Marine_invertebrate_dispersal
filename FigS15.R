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
datPD<-dat[,c(7,8,11:14)]
datPD<-datPD[complete.cases(datPD),]

# Exclude genera that are not in the matrix
common_genera<-intersect(factor(datPD$Genus), names1)
Matrix<-Matrix[which(names1 %in%  common_genera),which(names1 %in%  common_genera)]
datPD<-datPD[datPD$Genus %in% common_genera,]
datPD<-subset(datPD, !datPD$PlanktonicTime==0)
datPD$LogPlank<-log10(datPD$PlanktonicTime)



# Get prior
prior1<-get_prior(LogPlank~Hemisphere:absLatitude+absLatitude*DevelopmentalMode+
                    log10(EggSize)*DevelopmentalMode,
                  family=gaussian(),data=datPD)

# Fit model
model<-brm(LogPlank~Hemisphere:absLatitude+absLatitude*DevelopmentalMode+
             log10(EggSize)*DevelopmentalMode,
          data = datPD,
          family = gaussian(), 
          control = list(max_treedepth = 13, adapt_delta=0.99), chains=3, prior=prior1,thin = 5, iter=30000,
          warmup=5000, file='PlankDuration_no_phylo')

# Get a table with summary statistics
tidy(model, conf_level=0.95, chain=F)[1:7,]



# Save the posterior distributions of the fixed effects
Plank_post<-posterior_samples(model)[1:7]
#write.table(Plank_post, 'PLD_posterior_no_phylo.csv', sep=',', row.n=F)


# Plot predictions

datPD<-dat[,c(7:14)]
datPD<-datPD[complete.cases(datPD),]
datPD<-datPD[datPD$Genus %in% common_genera,]
datPD<-subset(datPD, !datPD$PlanktonicTime==0)
datPD$LogPlank<-log10(datPD$PlanktonicTime)


new<-with(datPD[datPD$DevelopmentalMode=='P',], expand.grid('Latitude'=seq(min(Latitude), max(Latitude), by=1), 
                                                  'Longitude'=seq(-180, 180, by=1),'EggSize'=mean(EggSize),
                                                  'DevelopmentalMode'='P'))
new$absLatitude<-abs(new$Latitude)
new$Hemisphere<-ifelse(new$Latitude<0,'S','N')
pr<-as.data.frame(fitted(model, new))
new$fit<-pr$Estimate


world<-map_data(map='world')
world<-map_data("world") %>% 
  filter(region != "Antarctica")



Q2a<-ggplot()+geom_raster(data=new, aes(x=Longitude, y=Latitude, fill=fit))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group))+theme_classic()+
  scale_fill_gradient( high='#990000', low='#FFCCCC',breaks=c(log10(30),log10(25),log10(20),log10(15),log10(10)),
                       labels=c(30,25,20,15,10))+
  geom_count(data=datPD[datPD$DevelopmentalMode=='P',], aes(x=Longitude,y=Latitude),  col='white', alpha=0.4, show.legend = F)+
  guides(fill=guide_legend(title='PD \n(days)'))+
  theme(axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'),
        text=element_text(size=18), plot.title=element_text(hjust=0.5))+
  ggtitle('a- Planktotrophic larvae')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  xlab('')+ylab('')
Q2a


new<-with(datPD[datPD$DevelopmentalMode=='L',], expand.grid('Latitude'=seq(min(Latitude), max(Latitude), by=1), 
                                                  'Longitude'=seq(-180, 180, by=1),'EggSize'=mean(EggSize), 'DevelopmentalMode'='L'))
new$absLatitude<-abs(new$Latitude)
new$Hemisphere<-ifelse(new$Latitude<0,'S','N')
pr<-as.data.frame(fitted(model, new, re_formula=NA))
new$fit<-pr$Estimate



Q2b<-ggplot()+geom_raster(data=new, aes(x=Longitude, y=Latitude, fill=fit))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group))+theme_classic()+
  scale_fill_gradient(low='#FFCCCC', high='#990000', breaks=c(log10(10),log10(8), log10(6),log10(4), log10(2)),
                      labels=c('10','8', '6', '4','2'))+
  geom_count(data=datPD[datPD$DevelopmentalMode=='L',], aes(x=Longitude,y=Latitude),  col='white', alpha=0.4, show.legend = F)+
  guides(fill=guide_legend(title='PD \n(days)'))+
  theme(axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'),
        text=element_text(size=18), plot.title=element_text(hjust=0.5))+
  ggtitle('b- Lecithotrophic larvae')+
  scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                     labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))+
  scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                     labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
  xlab('')+ylab('')
Q2b



FigS15<-grid.arrange(Q2a, Q2b, ncol=2)

FigS15



