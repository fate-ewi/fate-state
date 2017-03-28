library(ewidata)
library(knitr)
library(reshape2)
library(rstan)

devtools::install_github("fate-ewi/bayesdfa")
library(bayesdfa)

mcmc_iter = 4000
max_trends = 4
mcmc_chains = 3
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

regions = unique(ewidata$system)

for (i in 1:length(regions)) {
  # if the file doesn't exist, run the models
  if (!file.exists(paste0("Ecosystem_", regions[i], ".rds"))) {
    sub_data = ewidata[ewidata$system %in% regions[i],]
    # reshape data
    melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
    Y <- dcast(melted, code ~ year)
    names = Y$code
    Y = as.matrix(Y[,-which(names(Y) == "code")])
    
    # do the trend search, and save the table of model selection, along with the best model. By default, this isn't comparing student-t
    # versus normal models, but just estimating the student-t df parameter
    dfa_summary = find_dfa_trends(
      y = Y,
      kmax = min(max_trends, nrow(Y)),
      iter = mcmc_iter,
      compare_normal = FALSE,
      variance = c("unequal", "equal"),
      chains = mcmc_chains
    )
    saveRDS(dfa_summary, file = paste0("Ecosystem_", regions[i], ".rds"))
    
    # Make default plots (currently work in progress)
    pdf(paste0("Ecosystem_", regions[i], "_plots.pdf"))
    rotated = rotate_trends(dfa_summary$best_model)
    # trends
    plot_trends(rotated_modelfit, years = colnames(Y))
    # loadings
    plot_loadings(rotated_modelfit)
    # predicted values with data
    plot_fitted(dfa_summary$best_model)
    dev.off()
    
  }
}


