#DFTD phylodynamics code 

#Written by Nick Fountain-Jones (nfountainjones@gmail.com)
library(treeio)
library(ggtree)
library(tidyverse)
library(treestructure) 
library(skygrowth)
library(coda)
library(ggmcmc)
library(phylodyn)

#------------------------------------------------------------------------
##################Visualize MCC Tree##################
#------------------------------------------------------------------------
#read a .tre file
beast <- read.beast('WP_tumours_StrctC_GMRF_Test8.tre')
str(beast)
#can read nexus too


#plot the tree

ph <- ggtree(beast) + geom_tiplab(size=2)+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) 
ph  
get_taxa_name(tree_view = NULL, node = NULL) #tree has to be in the viewer


#doesn't work as yet
beastPoly <- as.polytomy(beast@phylo, feature='posterior', fun=function(x) as.numeric(x) < 50)
plot(beastPoly)

#convert into a newick
new <- read.nexus("WP_tumours_StrctC_GMRF_Test8.tre")
write.tree(new, file = "DFTD.newick")

DFTD <- read.tree("DFTD.newick")
  #------------------------------------------------------------------------
##################look for non-random tree structure using treestructure (Volz et al)##################
#------------------------------------------------------------------------
#not used in the manuscript 

#try min 10 for clade size 
treeSt <-  trestruct(DFTD, minCladeSize = 10, minOverlap = -Inf, nsim = 10000,
                     level = 0.05, ncpu = 1, verbosity = 1)
treeSt_df <- as.data.frame(treeSt)

plot(treeSt, use_ggtree = TRUE)

#------------------------------------------------------------------------
##################Look at NE through time for various linneages (Karcher et al)##################
#------------------------------------------------------------------------  

# skygrowth BP
fit <- skygrowth.map(DFTD ,  
                     , res = 2*15  # Ne changes every 6 months over a 14 year period
                     , tau0 = .1    # Smoothing parameter. If prior is not specified, this will also set the scale of the prior
)
plot(fit)
growth.plot(fit)

#fit with mcmc - I get a different result here than with skygrid map - not sure why yet
mcmcfit <- skygrowth.mcmc(DFTD, res = 2*15, tau0=.1 ) #not sure how to tune this
plot( mcmcfit )  + scale_y_log10(limits=c(.01, 1e5))
growth.plot( mcmcfit )


#add covariates
cov <- read.csv("DFTF_cov.csv", h=T)

#fix dates
library(lubridate)
cov$date <- dmy(cov$date) 

library(DataExplorer)
plot_missing(cov) ## Are there missing values, and what is the missing data profile? None in this case
cov_noNA<- drop_na(cov)
plot_missing(cov_noNA)


str(cov)

#plot covariates

FOIplot <- ggplot(data=cov_noNA, aes(x=date, y=forceMean)) +
  geom_line(colour="red",size=1) + geom_point()+theme_bw()+scale_x_date(breaks = '2 years')
FOIplot

Prevplot <- ggplot(data=cov_noNA, aes(x=date, y=prevalenceMean)) +
  geom_line(colour="black",size=1) + geom_point()+theme_bw()+scale_x_date(breaks = '2 years')
Prevplot

PopSIzeplot <- ggplot(data=cov_noNA, aes(x=date, y=Nestimate)) +
  geom_line(colour="blue",size=1) + geom_point()+theme_bw()+scale_x_date(breaks = '2 years')
PopSIzeplot

#convert date to a decimel
cov_noNA$date <- decimal_date(cov_noNA$date)
 
cov2017 <- dplyr::filter(cov_noNA, time < 2018) #only genomic data from DFTD until 2017
str(cov2017)


#this reduces it to year rather than every 6 months
covS <- cov2017%>% group_by(time) %>% summarise_all( list(mean))
str(covS)

#plot averages
FOIplotAvg <- ggplot(data=covS, aes(x=time, y=forceMean)) +
  geom_line(colour="red",size=1) + geom_point()+theme_bw()
FOIplotAvg
PrevplotAvg <- ggplot(data=covS, aes(x=time, y=prevalenceMean)) +
  geom_line(colour="black",size=1) + geom_point()+theme_bw()
PrevplotAvg
PopSizeplotAvg <- ggplot(data=covS, aes(x=time, y=Nestimate)) +
  geom_line(colour="black",size=1) + geom_point()+theme_bw()
PopSizeplotAvg


#------------------------------------------------------------------------
##################Covariate analysis##################
#------------------------------------------------------------------------  


#scale them

cov2017$prevalenceMean <- scale(cov2017$prevalenceMean)
cov2017$forceMean <- scale(cov2017$forceMean)
cov2017$Nestimate <- scale(cov2017$Nestimate)


#last sampling time
DLS <- 2017.3

set.seed(123)
covS$date <- NULL
str(covS)


# due to lack of signal the analyes did not zonverge so were not included in the manuscript.

covFitNeAvg <- skygrowth.mcmc.covar(DFTD,~Nestimate ,covS,maxSampleTime=DLS, res=22, iter=80000000, quiet=T) #iter0=thinging (10 by default)
summary(covFitNeAvg$beta)
neAvgMCMC <- as.mcmc(as.data.frame(covFitNeAvg$beta))
effectiveSize(neAvgMCMC)#time series length N, square  root{var x)}/n
#make a ggmcmc object

ggBNe <- ggs(neAvgMCMC)
HPDinterval(neAvgMCMC)
ggs_traceplot(ggBNe)

covFitFOI <- skygrowth.mcmc.covar(DFTD,~forceMean ,covS,maxSampleTime=DLS, res=45, iter=80000000, quiet=T) #iter0=thinging (10 by default)
summary(covFitFOI$beta)
save(covFitNeAvg, file="covFitNeAvg.Rdata")
load("covFitPrev.Rdata")
growth.plot(covFitFOI)

b <- as.mcmc(as.data.frame(covFitFOI $beta))
effectiveSize(b)#time series length N, square  root{var x)}/n

#make a ggmcmc object
ggB <- ggs(b)
HPDinterval(b)
ggs_density(ggB)
ggs_traceplot(ggB)


covFitPrevAvg <- skygrowth.mcmc.covar(DFTD,~prevalenceMean,covS,maxSampleTime=DLS, res=22, iter=80000000, quiet=T) 
save(covFitNeAvg, file="covFitNeAvg.Rdata")
load("covFitPrev.Rdata")
growth.plot(covFitFOI)
bprev <- as.mcmc(as.data.frame(covFitPrevAvg$beta))
effectiveSize(bprev)#time series length N, square  root{var x)}/n

#make a ggmcmc object
ggB <- ggs(bprev )
HPDinterval(bprev )
ggs_density(ggB)
ggs_traceplot(ggB)


#for final model
growth.plot(covFitM)
plot(covFitM)
summary(covFitM$beta)

plot(ts(covFitM$beta))


#---------------------------------------------------------------------------------------------
## Vector Autoregression Analysis
#---------------------------------------------------------------------------------------------
#Covariate data from May and November only. Im not sure if there is enough df s here considering 23 data points and 16 lagged variables

targetMonths<- c(5,11)
targetYears <- c(2006:2017)
cov_NovMay <- cov %>% filter(month %in% targetMonths) %>% filter(time %in% targetYears)

#we will use August covariate data in 2006 as well - as it is a close match to the DFTD genetic diversity data
cov_NovMay_extraMonth <- rbind(cov[1,],cov_NovMay )

res2006_2017 <- res[18:40,]

#no DFTD genetic diversity data from November 2017

data6month <-cbind(DFTD_diversity = res2006_2017[,2], cov_NovMay_extraMonth[1:23,])

#add var packages

library(tseries)
library(vars)
library(forecast)
library(TSstudio)


divTS <- ts(data6month$DFTD_diversity, start = c(2006,1), end= c(2017,1), frequency = 2)
ts_plot(divTS ) 
pp.test(divTS ) #no stationary - not that it matterns for VAR

neTS <- ts(data6month$Nestimate, start = c(2006,1),end= c(2017,1), frequency = 2)
ts_plot(neTS ) 
pp.test(neTS )

prevTS <- ts(data6month$prevalenceMean, start = c(2006,1),end= c(2017,1), frequency = 2)
ts_plot(prevTS ) 
pp.test(prevTS )

forceTS <- ts(data6month$forceMean, start = c(2006,1),end= c(2017, 1), frequency = 2)
ts_plot(forceTS)
pp.test(forceTS ) #stationary as p > 0.05 and we can reject the null hypothesis

#construct simplified ts data frame

data_v1 <- cbind(neTS, divTS , prevTS, forceTS)
colnames(data_v1) <- cbind("devilPopSize","DFTDgeneticDiv", "DFTDprev", "DFTDFOI")


#cross correlogram

ccf_v1 <- forecast::ggCcf(divTS , neTS)+theme_bw()
#positive values of k indicate that future values of y present y values on future y values. 
#negative values infdicates historic values of x shape current values of y

ccf_v2 <- forecast::ggCcf(forceTS , neTS)+theme_bw()

ccf_v3 <- forecast::ggCcf(prevTS , neTS)+theme_bw()

ccf_v4 <- forecast::ggCcf(forceTS , divTS)+theme_bw()

#correlation heatmap

corr <- round(cor(data_v1), 1)

library(ggcorrplot)

ggcorrplot(corr, hc.order = TRUE, type = "lower",
           outline.col = "white",  lab = TRUE)
