# setwd("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v28/")
source("../dream5_ecoli_setup.R")

dir.create("logs")

# This is for internal use to develop new experiments case by case.
do_one(18,  omit_knockout_evaluation_samples = F, reuse_results = F)
parallel::mclapply(which(conditions$proportion_removed>0), 
                   do_one, omit_knockout_evaluation_samples = F, reuse_results = F, mc.cores = parallel::detectCores())

# This is for internal use to get fast results when only updating the gold standards and not the knockoff procedure, 
parallel::mclapply(nrow(conditions):1, do_one, omit_knockout_evaluation_samples = F, reuse_results = T, mc.cores = parallel::detectCores())
parallel::mclapply(nrow(conditions):1, do_one, omit_knockout_evaluation_samples = T, reuse_results = T, mc.cores = parallel::detectCores())



# # The full run.
# parallel::mclapply(nrow(conditions):1, do_one, omit_knockout_evaluation_samples = F, reuse_results = F, mc.cores = parallel::detectCores())
# parallel::mclapply(nrow(conditions):1, do_one, omit_knockout_evaluation_samples = T, reuse_results = F, mc.cores = parallel::detectCores())

