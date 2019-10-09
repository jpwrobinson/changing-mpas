

library(rethinking)
library(tidyverse)
source('scripts/scaling_function.R')

## Fit bayesian linear regression to each site, for each metric, and extract posterior median of slope
## From raw data biom, focal is scaled data fitted to model m1, all provided in results/biomass_model.Rdata
## biom is the mean biomass at each reef, per year, with regime (recovering/shifted) and protection (fished/no-take)

## scale and center response covariates
focal<-scaler(biom, ID = c('diff','Year', 'Location', 'base', 'site.year'))
## rescale recovery years from 1998 bleaching
focal$Year <- focal$Year - 2004

# fit model with reef-specific temporal trends, and regime and fishing fixed effects
m1 <- map2stan(
	alist(
	    diff ~ dnorm( mu , sigma ) ,
	    mu <- a + ar[Location] + 
	    				(bA)*Shifted.Recovering.dummy + 
	    				(bB+br)*Year + 
	    				# (bB2+br)*Year^2 + 
	    				(bC)*Shifted.Recovering.dummy*Year +
	    				(bD)*Fished.Protected.dummy +
	    				(bE)*Fished.Protected.dummy*Year,
	    a ~ dnorm(15, 10),
	    c(ar, br)[Location] ~  dmvnorm2(0,sigmar, Rho), 
	    c(bA, bB, bC, bD, bE) ~ dnorm(0, 5),
	    c(sigma, sigmar) ~ dcauchy( 0 , 2 ),
	    Rho ~ dlkjcorr(2)
), data=focal, warmup =1500, iter=7000, chains=3)

save(m1, biom, focal, file='results/biomass_model.Rdata')

## make posterior predictions over 2005-2014 for each regime and protection/fished 

pred.dat<-expand.grid(Year = unique(focal$Year), Location = factor('Cousin Carbonate'), 
		Shifted_Recovering_dummy = unique(focal$Shifted.Recovering.dummy),
		Fished_Protected_dummy = unique(focal$Fished.Protected.dummy))
pred.dat$state<-ifelse(pred.dat$Shifted_Recovering_dummy == 0, 'Shifted', 'Recovering')
pred.dat$Management<-ifelse(pred.dat$Fished_Protected_dummy == 0, 'Fished', 'Protected')
pred.dat$group <- with(pred.dat, paste(state, Management))

# replace varying intercept samples with zeros
# 1000 samples by 21 sites
a_location_zeros <- matrix(0,1000,21)
mu<-link(m1, data=pred.dat, n=1000, replace=list(ar=a_location_zeros))
mu.mean <- apply( mu , 2 , mean )
mu.ci <- apply( mu , 2 , PI )
pred.dat$pred<-mu.mean
pred.dat$pred.lower<-mu.ci[1,]
pred.dat$pred.upper<-mu.ci[2,]


## now posterior distributions for each fixed effect - year, management + state
postyear <-as.data.frame(extract.samples(m1, 1000)$bB); colnames(postyear)<-c('shifted.fished')
postyear$shifted.protected <-as.data.frame(extract.samples(m1, 1000)$bB)[,1] + as.data.frame(extract.samples(m1, 1000)$bE)[,1]
postyear$recovering.fished <-as.data.frame(extract.samples(m1, 1000)$bB)[,1] + as.data.frame(extract.samples(m1, 1000)$bC)[,1]
postyear$recovering.protected <-as.data.frame(extract.samples(m1, 1000)$bB)[,1] + as.data.frame(extract.samples(m1, 1000)$bC)[,1] + as.data.frame(extract.samples(m1, 1000)$bE)[,1]
postyear<- postyear %>% gather(param, mu)
