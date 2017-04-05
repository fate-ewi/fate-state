library(Rcpp)
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

dfa_data = data.frame(
  grep = c(".LAT", "BLKI.", "CHI.", "PAV.", ".CO$", ".PI$", "ICY.", 
    "ICH.", "SEWARDLINE_", "AKCLIM_EBS_", "AKCLIM_GOA_", "calcofiMEAN.", "NLCOPE.","SBRD"),
  names = c(
    "Latitude",
    "BLKI",
    "CHI",
    "PAV",
    "Coho",
    "Pink",
    "ICY",
    "AK_Ichthyoplankton",
    "Seward Line",
    "AK Climate - GOA",
    "AK Climate - EBS",
    "CALCOFI",
    "NCC copepods",
    "CC_SBRD"
  )
)

for (i in 1:nrow(dfa_data)) {
  # if the file doesn't exist, run the models
  if (!file.exists(paste0(dfa_data$names[i], ".rds"))) {
    sub_data = ewidata[grep(dfa_data$grep[i], ewidata$code),]
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
    saveRDS(dfa_summary, file = paste0(dfa_data$names[i], ".rds"))
  
    # Make default plots (currently work in progress)
    pdf(paste0(dfa_data$names[i], "_plots.pdf"))
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
    print(plot_fitted(dfa_summary$best_model,names=names) + 
            theme(strip.text.x = element_text(size = 6)))
    dev.off()
}

}
