


# Test our best setting using the DREAM5 new results (53 new predictions explicitly tested)
explicitly_tested = 
  read.csv("omit_knockout_evaluation_samples=FALSE/seed=1/shrinkage_param=0.001/condition_on=pert_labels_plus_pca50/address_genetic_perturbations=TRUE/knockoff_type=glasso/chip_augmented/results_with_evaluation.csv.gz") %>% 
  merge(subset(ecoli_networks$dream5_new_test, is_positive_control=="no"),
        by = c("Gene1_name", "Gene2_name"),
        all.x = F,
        suffixes = c("chip_augmented", "DREAM5_new_tests"))
explicitly_tested %>% 
  subset(rlookc::knockoffQvals(knockoff_stat)<0.5) %>% 
  extract2("is_confirmedDREAM5_new_tests") %>% 
  table

# Record the top false positives in the most successful setting
# For later manual inspection
read.csv("omit_knockout_evaluation_samples=FALSE/seed=1/shrinkage_param=0.001/condition_on=pert_labels_plus_pca50/address_genetic_perturbations=TRUE/knockoff_type=glasso/chip_augmented/results_with_evaluation.csv.gz") %>%
  subset( !is_confirmed) %>%
  dplyr::arrange(q) %>% 
  write.csv("top_false_positives.csv")

# Enable direct comparison to dream5 consensus method.
# How many findings at q<0.5 and how many are supported by the benchmark network used in dream5?
read.csv("omit_knockout_evaluation_samples=FALSE/seed=1/shrinkage_param=0.001/condition_on=pert_labels_plus_pca50/address_genetic_perturbations=TRUE/knockoff_type=glasso/dream5/results_with_evaluation.csv.gz") %>%
  subset(q<0.5) %>% 
  extract2("is_confirmed") %>%
  fillna("unknown") %>%
  table %>%
  t %>%
  t %>%
  add_totals("total", colSums)

