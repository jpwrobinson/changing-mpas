

library(tidyverse)
library(rethinking)
source('scripts/scaling_function.R')

## From raw data fg2.biom, focal and focal94 are scaled data fitted to models m1 and m2, respectively, all provided in results/FG_models.Rdata

## Carn vs. Herb biomass
## Fit bayesian linear regression, with regime, year and protection as fixed effects, estimated separately for each FG (herbivore or carnivore)
focal<-scaler(fg2.biom[fg2.biom$Year != 1994,], ID = c('biom', 'Location','site.year', 'base'))

m1 <- map2stan(
	alist(
	    biom ~ dgamma2( mu , scale ) ,
	    log(mu) <- A +
	    				bA*Shifted.Recovering.dummy + 
	    				bB*Year + 
	    				bC*Shifted.Recovering.dummy*Year +
	    				bD*Fished.Protected.dummy +
	    				bE*Fished.Protected.dummy*Year,
		A <- a + a_loc[Location] + a_fg[FG2],
		bA <- ba + ba_fg[FG2],
		bB <- bb + bb_fg[FG2],
		bC <- bc + bc_fg[FG2],
		bD <- bd + bd_fg[FG2],
		bE <- be + be_fg[FG2],
	    
	    ## adaptive priors
	    c(a_loc)[Location] ~  dnorm(0,sigmar),
	    c(a_fg, ba_fg, bb_fg, bc_fg, bd_fg, be_fg)[FG2] ~  dmvnorm2(0,sigmar_fg, Rho),

	    ## fixed priors
	    a ~ dnorm(5.46, 10),
	    c(ba, bb, bc, bd, be) ~ dnorm(0, 3),
	    c(scale, sigmar, sigmar_fg) ~ dcauchy( 0 , 2 ),
	    Rho ~ dlkjcorr(4)
), data=focal, warmup=1500, iter=7000, chains=3)


## now fit to 1994 data only, without regime fixed effect
focal94<-scaler(fg2.biom[fg2.biom$Year == 1994,], ID = c('biom', 'Location','site.year', 'base'))

m2 <- map2stan(
	alist(
	    biom ~ dgamma2( mu , scale ) ,
	    log(mu) <- A + bD*Fished.Protected.dummy,
		A <- a + a_loc[Location] + a_fg[FG2],
		bD <- bd + bd_fg[FG2],
	    ## adaptive priors
	    c(a_loc)[Location] ~  dnorm(0,sigmar),
	    c(a_fg, bd_fg)[FG2] ~   dmvnorm2(0,sigmar_fg, Rho),

	    ## fixed priors
	    a ~ dnorm(5.26, 10),
	    bd ~ dnorm(0, 3),
	    c(scale, sigmar, sigmar_fg) ~ dcauchy( 0 , 2 ),
	    Rho ~ dlkjcorr(4)
), data=focal94, warmup=1500, iter=7000, chains=3)


save(m1,m2, fg2.biom, focal, focal94, file='results/FG_models')

### For posteriors of FG biomass by each state and management, ignoring year
pred.dat<-expand.grid(Year = 0,
		Location = factor('Cousin Carbonate'), 
		Shifted_Recovering_dummy = unique(focal$Shifted.Recovering.dummy),
		Fished_Protected_dummy = unique(focal$Fished.Protected.dummy),
		FG2=unique(focal$FG2))

pred.dat$state<-ifelse(pred.dat$Shifted_Recovering_dummy == 0, 'Shifted', 'Recovering')
pred.dat$Management<-ifelse(pred.dat$Fished_Protected_dummy == 0, 'Fished', 'Protected')
pred.dat$group <- with(pred.dat, paste(state, Management, FG2))

# replace varying intercept samples with zeros
# 1000 samples by 21 sites
a_location_zeros <- matrix(0,1000,21)
mu<-link(m1, data=pred.dat, n=1000, replace=list(ar=a_location_zeros))
mu.dat<-data.frame(group = rep(pred.dat$group, each = 1000))
mu.dat$state<-pred.dat$state[match(mu.dat$group, pred.dat$group)]
mu.dat$Management<-pred.dat$Management[match(mu.dat$group, pred.dat$group)]
mu.dat$FG2<-pred.dat$FG2[match(mu.dat$group, pred.dat$group)]
mu.dat$mu<-gather(data.frame(mu$mu))$value


## Repeat for 1994
pred.dat<-expand.grid(Location = factor('Cousin Carbonate'),
			Fished_Protected_dummy = unique(focal$Fished.Protected.dummy), 
			FG2=unique(focal$FG2))
pred.dat$Management<-ifelse(pred.dat$Fished_Protected_dummy == 0, 'Fished', 'Protected')
pred.dat$group <- with(pred.dat, paste(Management, FG2))

# replace varying intercept samples with zeros
# 1000 samples by 21 sites
a_location_zeros <- matrix(0,1000,21)
mu<-link(m2, data=pred.dat, n=1000, replace=list(ar=a_location_zeros))

mu94<-data.frame(group = rep(pred.dat$group, each = 1000))
mu94$Management<-pred.dat$Management[match(mu94$group, pred.dat$group)]
mu94$FG2<-pred.dat$FG2[match(mu94$group, pred.dat$group)]
mu94$mu<-gather(data.frame(mu$mu))$value


save(mu94, mu.dat, file = 'results/posterior_FG_Fig3.Rdata')
