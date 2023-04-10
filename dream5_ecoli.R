# dir.create("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v32/")
# setwd("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v32/")
source("../dream5_ecoli_setup.R")
dir.create("logs")

# This is for internal use to develop new experiments case by case.
# do_one(nrow(conditions),  omit_knockout_evaluation_samples = F, reuse_results = F)

# The full run.
parallel::mclapply(nrow(conditions):1, do_one, omit_knockout_evaluation_samples = F, reuse_results = F, mc.cores = parallel::detectCores()-1)
parallel::mclapply(nrow(conditions):1, do_one, omit_knockout_evaluation_samples = T, reuse_results = F, mc.cores = parallel::detectCores()-1)


# For the appendix: Gaussian mixture models don't fit these data well. Better BIC for a single Gaussian.
compare_mixture_models(X = ecoli_tf_expression)
# Compare also GeneNets, a linear + Gaussian competitor
source("../dream5_ecoli_genets.R")
source("../dream5_ecoli_plots.R")
