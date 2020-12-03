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

########################################################
#-----------Epidemilogical covariates------------------
########################################################

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

#------------------------------------------------------------------------
##################Look at NE through time for various linneages (Karcher et al)##################
#------------------------------------------------------------------------  

b0 <- BNPR(DFTD)
plot_BNPR( b0 )


########################################################
#---Correlations with epidemiological features---------
########################################################

res <- cbind(rev(2018-b0$x),rev(b0$effpop)) #rev - reverse element of effective population size

#gr=res[,2];gr=diff(gr)/gr[-length(gr)]#calculate the growth rate
#res=cbind((res[-1,1]+res[-nrow(res),1])/2,gr)

#interpolate effective population size for years with epi characteristics
xs <- seq(2006,2017)
DFTD_EffectivePopSize<- approx(res[,1],res[,2],xs,rule=2)$y #interpolating estimate for each year

CovDataCombined <- cbind(covS, DFTD_EffectivePopSize)

#correlation between devil population size and DFTD population size
cor(CovDataCombined$Nestimate,CovDataCombined$DFTD_EffectivePopSize)
cor.test(CovDataCombined$Nestimate,CovDataCombined$DFTD_EffectivePopSize)
ccf(CovDataCombined$Nestimate,CovDataCombined$DFTD_EffectivePopSize)

#calculate and plot correlations
ReducedCovar <- cbind( devilPopSize=CovDataCombined$Nestimate, Prev=CovDataCombined$prevalenceMean, FOI=CovDataCombined$forceMean, DFTDEffPopSize=CovDataCombined$DFTD_EffectivePopSize)
cormat <- round(cor(ReducedCovar),2)
library(reshape2)
melted_cormat <- melt(cormat)

library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Final Heatmap

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

library(ggpubr)
cor1 <- ggscatter(as.data.frame(ReducedCovar) , x = "devilPopSize", y = "DFTDEffPopSize",
                  add = "reg.line",  # Add regression line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
cor1 + stat_cor(method = "pearson", label.x = 1, label.y = 2)

cor2 <- ggscatter(as.data.frame(ReducedCovar) , x = "devilPopSize", y =  "FOI",
                     add = "reg.line",  # Add regressin line
                     add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                     conf.int = TRUE # Add confidence interval
                     
)
# Add correlation coefficient
cor2 + stat_cor(method = "pearson", label.x = 0.8, label.y = 2)


########################################################
#-------------------Skygrowth model---------------------
########################################################


globalgrowth <- skygrowth.mcmc(DFTD, res = 2*14, tau0=0.1,tau_logprior = function (x) dexp(x,0.1,T), mhsteps= 1e+06, control=list(thin=1e3) ) 

#check convergence
globalMCMC <- as.mcmc(cbind(globalgrowth$growthrate[,1:(ncol(globalgrowth$growthrate)-1)],globalgrowth$ne,globalgrowth$tau))
effectiveSize(globalMCMC)

#plots

growth.plot(globalgrowth)+theme_bw()
neplot(globalgrowth)+theme_bw()


#last sampling time
DLS <- 2017.3

set.seed(123)
covS$date <- NULL
str(covS)


growth.plot(covFitFOI)

b <- as.mcmc(as.data.frame(covFitNe$beta))
effectiveSize(b)#time series length N, square  root{var x)}/n

#make a ggmcmc object
ggB <- ggs(b)
HPDinterval(b)
ggs_density(ggB)
ggs_traceplot(ggB)

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