
library(ewidata)
library(knitr)
library(reshape2)
mcmc_iter = 4000
max_trends = 4
mcmc_chains = 3
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dfa_data = data.frame(
  grep = c(".LAT", "BLKI.", "CHI.", "PAV.", ".CO$", ".PI$", "ICY.", "ICH."),
  names = c(
    "Latitude",
    "BLKI",
    "CHI",
    "PAV",
    "Coho",
    "Pink",
    "ICY",
    "Ichthyoplankton"
  )
)

for (i in 1:nrow(dfa_data)) {
  # if the file doesn't exist, run the models
  if (!file.exists(paste0(dfa_data$names[i], ".rds"))) {
    sub_data = ewidata[grep(dfa_data$grep[i], ewidata$code), ]
    # reshape data
    melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
    Y <- dcast(melted, code ~ year)
    names = Y$code
    Y = as.matrix(Y[, -which(names(Y) == "code")])
    
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
  }
}