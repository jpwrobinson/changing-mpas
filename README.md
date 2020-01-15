# changing-mpas
Data, statistical models, and R scripts accompanying Graham et al. Changing Role of Coral Reef Marine Reserves in a Warming Climate (in revision).

[data](/data) contains .csv files of data underlying Figs 1-4 and S2-S4.

[scripts](/scripts) contains .R scripts with Bayesian model structures, code for estimating the posterior distributions plotted in main text figures, and functions for estimating jacknife confidence intervals and for standardizing datasets before model fitting.

[results](/results) contains .Rdata files with analysed datasets and fitted Bayesian models.

Data viz and analyses were run in R 3.6.0, using packages tidyverse, rethinking, Rstan and cowplot.