#setwd("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v24")
source("../dream5_ecoli_setup.R")
# Collect results to make certain key summary plots
cat("Loading experimental results and making UMAPs...")
conditions_with_summaries = conditions = read.csv("experiments_to_run.csv", row.names = 1)
exchangeability = umaps = knockoffs = calibration_sim = calibration_gs = list()
i = 0
for(condition_index in nrow(conditions):1){
  try({
    write.table(conditions[condition_index,], sep = "\t", quote = F)
    attach(conditions[condition_index,], warn.conflicts = F)
    working_directory = ""
    for(j in seq_along(conditions[condition_index,])){
      prefix = paste0( colnames(conditions)[[j]], "=", conditions[condition_index, j] )
      working_directory %<>% file.path(prefix, .)
    }
    # Check out the knockoffs via UMAP and KNN exchangeability test
    withr::with_dir(working_directory, {
      knockoffs[[condition_index]] = read.csv("knockoffs.csv", row.names = 1)
      exchangeability[[condition_index]] = rlookc::KNNTest(X = ecoli_tf_expression, X_k = knockoffs[[condition_index]])
      conditions_with_summaries[condition_index, "KNN_exchangeability_p"]          =  exchangeability[[condition_index]] %>% extract2("p_value")
      conditions_with_summaries[condition_index, "KNN_exchangeability_proportion"] =  exchangeability[[condition_index]] %>% extract2("prop_not_swapped")
      umaps[[condition_index]] = umap::umap(knockoffs[[condition_index]], random.state = 0) %>% extract2("layout")
      calibration_sim[[condition_index]] = read.csv("average_case_calibration.csv", row.names = 1) %>%
        colMeans %>%
        (function(x) data.frame(nominal_fdr = gsub("^X", "", names(x)) %>% as.numeric, empirical_fdr = x))
    })
    umaps[[condition_index]]           %<>% merge(conditions[condition_index,])
    calibration_sim[[condition_index]] %<>% merge(conditions[condition_index,])
    # check vs various gold standards
    for( gold_standard_name in AVAILABLE_GOLD_STANDARDS ){
      i = i + 1
      withr::with_dir(file.path(working_directory, gold_standard_name), {
        calibration_gs[[i]] = read.csv("ecoli_calibration.csv")
        calibration_gs[[i]] %<>% merge(conditions[condition_index,])
        calibration_gs[[i]][["gold_standard_name"]] = gold_standard_name
      })
    }
  })
}
umaps %<>% data.table::rbindlist()
calibration_gs %<>% data.table::rbindlist()
calibration_sim %<>% data.table::rbindlist()
conditions_with_summaries %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )
umaps %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )
calibration_gs %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )
calibration_sim %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )

# Vary confounder handling
calibration_gs %>%
  subset( T &
            gold_standard_name=="chip_augmented" &
            seed==1 &
            knockoff_method %in% c( "glasso_0.001")
  ) %>%
  ggplot() +
  geom_point(aes(x = nominal_fdr, y = empirical_fdr, color = condition_on, shape = condition_on)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  ggtitle("Calibration and counfounders") +
  facet_grid(knockoff_method~ifelse(address_genetic_perturbations, "Knockouts addressed", "Not addressed")) +
  scale_y_continuous(limits = 0:1)
ggsave("fig_confounders.pdf", width = 6, height = 3)

# same, on confounding, but make a simpler plot with error bars
fillna =function(x, filler){ x[is.na(x)] = filler; x}
calibration_gs %>%
  subset( T &
            gold_standard_name=="chip_augmented" &
            address_genetic_perturbations == T &
            seed == 1 &
            knockoff_method %in% c("glasso_0.001")
  ) %>%
  dplyr::mutate(moe_95 = fillna(moe_95, 1)) %>%
  ggplot(
    mapping = aes(x = nominal_fdr,
                  y = empirical_fdr,
                  ymin = pmax(empirical_fdr - moe_95, 0),
                  ymax = pmin(empirical_fdr + moe_95, 1))
  ) +
  geom_pointrange() +
  geom_abline(aes(slope = 1, intercept = 0)) +
  ggtitle("Calibration and counfounders") +
  facet_grid(knockoff_method~condition_on) +
  scale_y_continuous(limits = 0:1)
ggsave("fig_confounders_errorbars.pdf", width = 12, height = 3)

# Plot on confounding, but include power
fillna =function(x, filler){ x[is.na(x)] = filler; x}
calibration_gs %>%
  subset( T &
            gold_standard_name=="chip_augmented" &
            address_genetic_perturbations == T &
            seed == 1 &
            knockoff_method %in% c("sample", "glasso_0.001")
  ) %>%
  dplyr::mutate(moe_95 = fillna(moe_95, 1)) %>%
  ggplot(
    mapping = aes(x = nominal_fdr,
                  y = num_discoveries,
                  shape = condition_on, color = condition_on)
  ) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  ggtitle("Calibration and counfounders") +
  facet_grid(~knockoff_method)
ggsave("fig_confounders_power.pdf", width = 6, height = 3)

# Vary knockoff generation method (KNN exchangeability diagnostic)
conditions_with_summaries %>%
  subset(
      condition_on == "none" &
      seed==1 &
      !address_genetic_perturbations
  ) %>%
  head(8) %>%
  ggplot(mapping = aes(
    x = KNN_exchangeability_proportion,
    color = knockoff_method,
    label = knockoff_method,
    y = KNN_exchangeability_p
  )) +
  geom_vline(aes(xintercept = 0.5)) +
  geom_point() +
  ggrepel::geom_label_repel()
ggsave("fig_knockoff_type_knn.pdf", width = 6, height = 5)

# Vary knockoff generation method (real Y, calibration)
calibration_gs %>%
  subset(
    condition_on == "none" &
      seed==1 &
      !address_genetic_perturbations &
      gold_standard_name=="chip_augmented"
  ) %>%
  ggplot() +
  geom_point(aes(x = nominal_fdr, y = empirical_fdr, color = knockoff_method, shape = knockoff_method)) +
  geom_line(aes(x = nominal_fdr, y = empirical_fdr, color = knockoff_method, shape = knockoff_method)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  ggtitle("Calibration by knockoff method") +
  facet_grid(~address_genetic_perturbations) +
  scale_y_continuous(limits = 0:1)
ggsave("fig_knockoff_type.pdf", width = 6, height = 5)

# Vary knockoff generation method (real Y, AUPR)
calibration_gs %>%
  subset(
    condition_on == "none" &
      seed==1 &
      !address_genetic_perturbations
  ) %>%
  ggplot() +
  geom_point(aes(x = num_discoveries,
                 y = 1-empirical_fdr,
                 color = knockoff_method,
                 shape = knockoff_method)) +
  geom_line(aes(x = num_discoveries,
                 y = 1-empirical_fdr,
                 color = knockoff_method,
                 shape = knockoff_method)) +
  ggtitle("Enrichment vs various gold standards") +
  scale_y_continuous(limits = 0:1) +
  scale_x_log10() +
  xlab("Recall") +
  ylab("Precision") +
  facet_wrap(~gold_standard_name)
ggsave("fig_knockoff_type_aupr.pdf", width = 6, height = 5)

# Vary knockoff method (simulated Y)
calibration_sim %>%
  subset(condition_on == "none" & !address_genetic_perturbations & seed==1) %>%
  ggplot() +
  geom_line(aes(x = nominal_fdr, y = empirical_fdr, color = knockoff_method, shape = knockoff_method)) +
  geom_point(aes(x = nominal_fdr, y = empirical_fdr, color = knockoff_method, shape = knockoff_method)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  ggtitle("Calibration by method", subtitle = "Simulated target genes") +
  scale_y_continuous(limits = 0:1)
ggsave("fig_knockoff_type_sim.pdf", width = 6, height = 5)

# Fig number uncertain: various gold standards
calibration_gs %>%
  subset(condition_on == "none" &
           !address_genetic_perturbations &
           knockoff_type == "sample" ) %>%
  ggplot() +
  geom_line(aes(x = num_discoveries,
                y = 1-empirical_fdr,
                color = gold_standard_name,
                shape = gold_standard_name)) +
  ggtitle("Enrichment vs various gold standards") +
  scale_y_continuous(limits = 0:1) +
  scale_x_log10() +
  xlab("Recall") +
  ylab("Precision") +
  facet_wrap(~seed, ncol =2)
ggsave("fig_gs.pdf", width = 6, height = 5)

# Fig number uncertain: replicates with different random seed
calibration_gs %>%
  subset(condition_on == "none" &
           gold_standard_name == "chip_augmented" &
           !address_genetic_perturbations &
           knockoff_type == "sample" ) %>%
  ggplot() +
  geom_line(aes(x = nominal_fdr, y = empirical_fdr, group = seed)) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  ggtitle("Multiple knockoff realizations") +
  scale_y_continuous(limits = 0:1) +
  xlab("Recall") +
  ylab("Precision")
ggsave("fig_randomness.pdf", width = 6, height = 5)

# Visualize data and knockoffs: does the important structure seem to be captured?
umap_actual    = umap::umap(ecoli_tf_expression, random.state = 0) %>% extract2("layout")
umap_actual %<>% merge(data.frame("knockoff_type" = "real data",
                                  "knockoff_method" = "real data",
                                  "shrinkage_param" = NA,
                                  "address_genetic_perturbations" = F,
                                  "condition_on"  = "none"   ,
                                  "seed"    = 1))
umaps %>%
  rbind(umap_actual) %>%
  mutate(cluster = clusters$cluster %>% rep(times = dplyr::n()/805 )) %>%
  subset(condition_on == "none" & !address_genetic_perturbations & seed==1) %>%
  ggplot() +
  geom_point(aes(x = V1,  y = V2,color = as.character(cluster))) +
  coord_fixed() +
  xlab("UMAP1") + ylab("UMAP2") +
  facet_wrap(~knockoff_method, nrow = 2) +
  ggtitle("UMAP and K-means clusters from E. coli TF expression")
ggsave(("fig_umaps.pdf"), width = 10, height =6)

