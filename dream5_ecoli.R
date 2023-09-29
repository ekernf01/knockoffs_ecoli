# These lines are not needed when using the wrapper script run_on_aws, but they are useful for interacting with it in rstudio.
# dir.create("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v34/")
# setwd("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v34/")
source("../dream5_ecoli_gold_standards.R")
source("../dream5_ecoli_setup.R")
dir.create("logs")

successful_jobs = read.csv("successful_jobs.txt")[[1]]
# The full run.
parallel::mclapply(
  (1:nrow(conditions)) %>% setdiff(successful_jobs),
  do_one_safe, 
  reuse_results = F, 
  mc.cores = parallel::detectCores()-1
)

# For the appendix: Gaussian mixture models don't fit these data well. Better BIC for a single Gaussian.
compare_mixture_models(X = ecoli_tf_expression)

source("../dream5_ecoli_plots.R")
