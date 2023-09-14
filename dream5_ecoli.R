dir.create("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v33/")
setwd("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v33/")
source("../dream5_ecoli_gold_standards.R")
source("../dream5_ecoli_setup.R")
dir.create("logs")

# This is for internal use (interactive debugging).
# do_one(1,  reuse_results = F, test_mode = T)
do_one_safe(1,  reuse_results = F, test_mode = F)
do_one_safe(3,  reuse_results = F, test_mode = F)
do_one_safe(2,  reuse_results = F, test_mode = F)
do_one_safe(4,  reuse_results = F, test_mode = F)

# The full run.
parallel::mclapply(
  1:nrow(conditions),
  do_one_safe, 
  reuse_results = F, 
  mc.cores = parallel::detectCores()-1
)

# For the appendix: Gaussian mixture models don't fit these data well. Better BIC for a single Gaussian.
compare_mixture_models(X = ecoli_tf_expression)

source("../dream5_ecoli_plots.R")
