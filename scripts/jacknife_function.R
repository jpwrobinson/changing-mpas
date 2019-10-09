


## Jacknife function to estimate reserve response ratios uncertainty limits
## function refits ratio values, dropping one site at a time, looping through all sites

jacklims<-function(dataset, variable, ID1, ID2=NA, ID3 = NA, site.var) {

	library(tidyverse)
	# enquotes variables to be used with dplyr-functions
	variable <- enquo(variable)
	group.var1 <- enquo(ID1) ## this is management variable

	if(!missing(ID2)){
	group.var2 <- c(enquo(ID2), group.var1)} ## additional grouping variables, eg. FG

	if(!missing(ID3)){
	group.var2 <- c(enquo(ID3), enquo(ID2), group.var1)} ## additional grouping variables, eg. FG
	

	sites<-data.frame(unique(dataset[,site.var]))[,1]
	nsites <- length(unique(sites))

	ratio<-numeric()

	for(i in 1:nsites){

	dataset<-data.frame(dataset)
	dat<-dataset[!(dataset[,site.var] %in% sites[i]),]

	if(!missing(ID2)){
	jack <- dat %>% 
		group_by(!!!group.var2) %>%
		summarise(b = mean(!!variable)) %>% ungroup() %>%
		spread(!!group.var1, b) %>% 
		mutate(ratio = log(Protected / Fished), jack = i)}

	if(missing(ID2)){
	jack <- dat %>% 
		group_by(!!group.var1) %>%
		summarise(b = mean(!!variable)) %>% ungroup() %>%
		spread(!!group.var1, b) %>% 
		mutate(ratio = log(Protected / Fished), jack = i)}

	ratio<-rbind(ratio, jack)

	}
	ratio
}


## same function as jacklims, but using observed data (e.g. biomass or richness) not response ratios

jacklims_raw<-function(dataset, variable, ID1, ID2=NA, ID3 = NA, site.var) {

	library(tidyverse)
	# enquotes variables to be used with dplyr-functions
	variable <- enquo(variable)
	group.var1 <- enquo(ID1) ## this is management variable

	if(!missing(ID2)){
	group.var2 <- c(enquo(ID2), group.var1)} ## additional grouping variables, eg. FG

	if(!missing(ID3)){
	group.var2 <- c(enquo(ID3), enquo(ID2), group.var1)} ## additional grouping variables, eg. FG
	

	sites<-data.frame(unique(dataset[,site.var]))[,1]
	nsites <- length(unique(sites))

	net<-numeric()

	for(i in 1:nsites){

	dataset<-data.frame(dataset)
	dat<-dataset[!(dataset[,site.var] %in% sites[i]),]

	if(!missing(ID2)){
	jack <- dat %>% 
		group_by(!!!group.var2) %>%
		summarise(b = mean(!!variable)) %>% ungroup() %>%
		spread(!!group.var1, b) %>% 
		mutate(net = Protected - Fished, jack = i)}

	if(missing(ID2)){
	jack <- dat %>% 
		group_by(!!group.var1) %>%
		summarise(b = mean(!!variable)) %>% ungroup() %>%
		spread(!!group.var1, b) %>% 
		mutate(net = Protected - Fished, jack = i)}

	net<-rbind(net, jack)

	}
	net
}