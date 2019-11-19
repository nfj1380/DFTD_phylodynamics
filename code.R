#DFTD phylodynamics code 

#Written by Nick Fountain-Jones (nfountainjones@gmail.com)
library(treeio)
library(ggtree)
library(tidyverse)
library(treestructure) 
library(skygrowth)
library(coda)
library(ggmcmc)
library(treeio)
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

#try min 10 for clade size 
treeSt <-  trestruct(DFTD, minCladeSize = 10, minOverlap = -Inf, nsim = 10000,
                     level = 0.05, ncpu = 1, verbosity = 1)
treeSt_df <- as.data.frame(treeSt)

plot(treeSt, use_ggtree = TRUE)

#------------------------------------------------------------------------
##################Look at NE through time for various linneages (Karcher et al)##################
#------------------------------------------------------------------------  
#add covariates
cov <- read.csv("DFTF_cov.csv", h=T)
#fix dates
library(lubridate)
cov$date<- mdy(cov$date) 
cov$date <- decimal_date(cov$date) 

library(DataExplorer)
plot_missing(cov) ## Are there missing values, and what is the missing data profile? None in this case
cov_noNA<- drop_na(cov)
plot_missing(cov_noNA)

cov2017 <- dplyr::filter(cov_noNA, time < 2018) #only genomic data from DFTD until 2017

#this reduces it to year
covS <- cov2017%>% group_by(time) %>% summarise_all( list(mean))

str(cov)

#scale them

covS[,2:19] <- scale(covS[,2:19])


# skygrowth BP
fit <- skygrowth.map(DFTD ,  
                     , res = 2*14  # Ne changes every 6 months over a 14 year period
                     , tau0 = .1    # Smoothing parameter. If prior is not specified, this will also set the scale of the prior
)
plot(fit)
growth.plot(fit)

#fit with mcmc - I get a different result here than with skygrid map - not sure why yet
mcmcfit <- skygrowth.mcmc(DFTD, res = 2*14, tau0=.1 ) #not sure how to tune this
plot( mcmcfit )  #+ scale_y_log10(limits=c(.01, 1e5))
growth.plot( mcmcfit )


#compare to phylodynn

b0 <- BNPR(DFTD)
plot_BNPR( b0 )

#last sampling time
DLS <- 2017.3

set.seed(123)
covS$date <- NULL
str(covS)
covFitNe <- skygrowth.mcmc.covar(DFTD,~Nestimate,cov2017,maxSampleTime=DLS, res=50, iter=40000000) #iter0=thinging (10 by default)
summary(covFitFOI$beta)
save(covFitNe, file="covFitNe.Rdata")
load("covFitPrev.Rdata")
growth.plot(covFitFOI)

b <- as.mcmc(as.data.frame(covFitNe$beta))
effectiveSize(b)#time series length N, square  root{var x)}/n

#make a ggmcmc object
ggB <- ggs(b)
HPDinterval(b)
ggs_density(ggB)
ggs_traceplot(ggB)

#males
set.seed(123)
covFitM <- skygrowth.mcmc.covar(bFIVnoFR,~m,cov,maxSampleTime=DLS,res=50,iter=4000000,iter0=20, tau0 =10, quiet=T) #res= number of time points
# could alter tau log prior - number of time events
a <- as.mcmc(as.data.frame(covFitM$beta))
effectiveSize(a)#time series length N, square  root{var x)}/n

#make a ggmcmc object
ggA <- ggs(a)
HPDinterval(a)
ggs_density(ggA)
ggs_traceplot(ggA)

#for final model
growth.plot(covFitM)
plot(covFitM)
summary(covFitM$beta)

plot(ts(covFitM$beta))



str(b)