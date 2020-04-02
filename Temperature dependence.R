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



# Select only the data required
datPD<-dat[,c(7,8,11:14)]
datPD<-datPD[complete.cases(datPD),]

# Exclude genera that is not in the matrix
common_genera<-intersect(factor(datPD$Genus), names1)
Matrix<-Matrix[which(names1 %in%  common_genera),which(names1 %in%  common_genera)]
datPD<-datPD[datPD$Genus %in% common_genera,]
datPD<-subset(datPD, !datPD$PlanktonicTime==0)
datPD$LogPlank<-log10(datPD$PlanktonicTime)



# Get prior
prior1<-get_prior(LogPlank~Hemisphere:absLatitude+absLatitude*DevelopmentalMode+
                    log10(EggSize)*DevelopmentalMode+(1|Genus),
                  family=gaussian(),data=datPD)

# Fit model
model<-brm(LogPlank~Hemisphere:absLatitude+absLatitude*DevelopmentalMode+
             log10(EggSize)*DevelopmentalMode+(1|Genus),
           data = datPD,
           family = gaussian(), cov_ranef = list(Matrix),
           control = list(max_treedepth = 13, adapt_delta=0.99), chains=3, prior=prior1,thin = 5, iter=30000,
           warmup=5000, file='PlankDuration')



# Estimate temperature dependence following O'Connor et al. 2007

datPD<-dat[,c(7,8,11:14, 16)]
datPD<-datPD[complete.cases(datPD),]
datPD<-subset(datPD, !datPD$PlanktonicTime==0)
datPD$LogPlank<-log10(datPD$PlanktonicTime)

# Boltzmann constant
Boltz<-8.617333262145*10^(-5)

# Scaled inverse temperature
datPD$Arren<-1/(Boltz*(datPD$temperature_annual+273))


# Get prior
prior1<-get_prior(LogPlank~Arren+(1|Genus),
                  family=gaussian(),data=datPD[datPD$DevelopmentalMode=='P',])

model1<-brm(LogPlank~Arren+(1|Genus),
          data = datPD[datPD$DevelopmentalMode=='P',],
          family = gaussian(), cov_ranef = list(Matrix),
          control = list(max_treedepth = 15, adapt_delta=0.99), chains=3,
          prior=prior1,thin = 5, iter=30000, warmup=5000,
           file='TemperatureDependenceP')

prior2<-get_prior(LogPlank~Arren+(1|Genus),
                  family=gaussian(),data=datPD[datPD$DevelopmentalMode=='L',])

model2<-brm(LogPlank~Arren+(1|Genus),
          data = datPD[datPD$DevelopmentalMode=='L',],
          family = gaussian(), cov_ranef = list(Matrix),
          control = list(max_treedepth = 15, adapt_delta=0.99), chains=3,
          prior=prior1,thin = 5, iter=30000, warmup=5000,
          file='TemperatureDependenceL')



