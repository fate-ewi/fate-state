set.seed(123)
library(bayesdfa)
library(ewidata)
library(reshape2)
mcmc_iter=6000
max_trends=4
mcmc_chains=3
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

ich_data = ewidata[grep("ICH.",ewidata$code),]
# reshape data
melted = melt(ich_data[,c("code","year","value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y)=="code")])
# fit a model with equal variances
dfa_mod = fit_dfa(y=Y, num_trends=2, iter=mcmc_iter, 
  varIndx=rep(1,nrow(Y)), chains=mcmc_chains)

rotated = rotate_trends(dfa_mod)
#plot_trends(rotated, years=as.numeric(colnames(Y)))

p1 = plot_loadings(rotated, facet=FALSE, names=names, threshold=NULL)
# reduces to ~ 32-33
p2 = plot_loadings(rotated, facet=FALSE, names=names, threshold=0.8)
# reduces to ~ 
p3 = plot_loadings(rotated, facet=FALSE, names=names, threshold=0.85)
