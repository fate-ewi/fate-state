library(ewidata)
library(knitr)
library(reshape2)
library(rstan)

devtools::install_github("fate-ewi/bayesdfa")
library(bayesdfa)

mcmc_iter = 4000
max_trends = 3
mcmc_chains = 2
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#regions = unique(ewidata$system)
regions = "NCC"

for (i in 1:length(regions)) {
  # process data
  sub_data = ewidata[ewidata$system %in% regions[i],]
  # reshape data
  melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
  Y <- dcast(melted, code ~ year)
  names = Y$code
  Y = as.matrix(Y[,-which(names(Y) == "code")])
  
  # if the file doesn't exist, run the models
  if (!file.exists(paste0("Ecosystem_", regions[i], "3.rds"))) {
    # do the trend search, and save the table of model selection, along with the best model. By default, this isn't comparing student-t
    # versus normal models, but just estimating the student-t df parameter
    dfa_summary = find_dfa_trends(
      y = Y,
      kmax = min(max_trends, nrow(Y)),
      control = list(adapt_delta = 0.999, max_treedepth = 20),
      iter = mcmc_iter,
      compare_normal = FALSE,
      variance = c("unequal", "equal"),
      chains = mcmc_chains
    )
    saveRDS(dfa_summary, file = paste0("Ecosystem_", regions[i], "3.rds"))
  }
  
  if(file.exists(paste0("Ecosystem_", regions[i], "3.rds"))) dfa_summary = readRDS(file = paste0("Ecosystem_", regions[i], "3.rds"))
  
  # Make default plots (currently work in progress)
  pdf(paste0("Ecosystem_", regions[i], "_plots3.pdf"))
  rotated = rotate_trends(dfa_summary$best_model)
  # trends
  print(plot_trends(rotated, years = as.numeric(colnames(Y))))
  # loadings
  print(plot_loadings(rotated, names = names))
  if(ncol(rotated$Z_rot_mean)==2) {
    plot(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], col="white", 
      xlab="Loading 1", ylab = "Loading 2")
    text(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], names, cex=0.3)
    lines(c(-10,10),c(0,0))
    lines(c(0,0), c(-10,10))
  }
  # predicted values with data
  print(plot_fitted(dfa_summary$best_model, names=names) + 
    theme(strip.text.x = element_text(size = 6)))
  summary_table<-dfa_summary$summary
  capture.output(summary_table, file = paste0("Ecosystem_", regions[i], "_summary3.txt"))
  dev.off()  
  
}
