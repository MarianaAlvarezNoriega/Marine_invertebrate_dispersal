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

# Create a column with values of 1 for planktotrophic larvae and 0 for 
# lecithotrophic larvae
dat$Planktotrophic<-ifelse(dat$DevelopmentalMode=='P',1,0)

# Select only the data required
datPD<-dat[,c(7,8,11,14,19)]
datPD<-datPD[complete.cases(datPD),]

# Exclude genera that are not in the matrix
common_genera<-intersect(factor(datPD$Genus), names1)
Matrix<-Matrix[which(names1 %in%  common_genera),which(names1 %in%  common_genera)]
datPD<-datPD[datPD$Genus %in% common_genera,]


# Fit model

model<-brm(Planktotrophic~absLatitude+absLatitude:Hemisphere+(1|Genus),
          family=binomial(),
          data = datPD,seed=17,
          cov_ranef = list(Matrix),
          control = list(max_treedepth = 17, adapt_delta=0.999), chains=3,thin = 5, iter=30000, 
          warmup=5000, file='Planktotrophy')

# Compute R2 of the model
bayes_R2(model)

# Get a table with summary statistics
tidy(model, conf_level=0.95, chain=F)[1:3,]


# Get posterior distributions for fixed effects
Plank_post<-posterior_samples(model)[1:3]

#write.table(Plank_post, 'feeding_posterior.csv', sep=',', row.n=F)


datPD<-dat[,c(7,8:11,14,19)]
new<-with(datPD, expand.grid('Latitude'=seq(min(Latitude), max(Latitude), by=1), 
                             'Longitude'=seq(-180, 180, by=1)))
new$absLatitude<-abs(new$Latitude)
new$Hemisphere<-ifelse(new$Latitude<0,'S','N')
pr<-as.data.frame(fitted(model, new, re_formula=NA))
new$fit<-pr$Estimate


world<-map_data(map='world')
world<-map_data("world") %>% 
  filter(region != "Antarctica")



FigS8<-ggplot()+geom_raster(data=new, aes(x=Longitude, y=Latitude, fill=fit))+
  geom_polygon(data=world, aes(x=long, y=lat, group=group))+theme_classic()+
  scale_fill_gradient(high='#990000', low='#FFCCCC', breaks=c(0.71,0.70,0.69,0.68,67))+
  geom_count(data=datPD, aes(x=Longitude,y=Latitude, shape=DevelopmentalMode),  col='grey', alpha=0.4, show.legend = T)+
  guides(fill=guide_legend(title=''))+scale_size_continuous(name='')+
    ggtitle('')+scale_shape_discrete(name='')+
  theme(text = element_text(size=15), plot.title=element_text(hjust=1))+
  theme(axis.line = element_line(size = 1), axis.ticks.length = unit(0.25, 'cm'),
        text=element_text(size=15))+
    scale_y_continuous(breaks=c(-80, -60, -40,-20,0,20,40,60,80), 
                       labels=c('80°S','60°S',  '40°S', '20°S','0°', '20°N', '40°N',  '60°N','80°N'))+
    scale_x_continuous(breaks=c(-180,  -120, -60,0,60,120,180), 
                       labels=c('180°W', '120°W', '60°W', '0°', '60°E', '120°E', '180°E'))
FigS8




