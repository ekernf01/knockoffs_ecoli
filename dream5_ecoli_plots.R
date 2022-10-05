#setwd("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v28")
source("../dream5_ecoli_setup.R")
fillna =function(x, filler){ x[is.na(x)] = filler; x}
wrapFacetLabels = function(x){
  sapply(x, function(s){
    s %>%
      gsub( "_", " ", . ) %>%
      strwrap( width = 20 ) %>%
      paste0( collapse = "\n" )
  })
}
ggplot2::theme_update(text = element_text(family = "ArialMT"))
# Collect results to make certain key summary plots
# This will allow us to show calibration on any combination of gold standard, knockoff construction method,
# number of principal components, et cetera by subsetting a giant tidy dataframe. 
cat("Loading experimental results and making UMAPs...")
conditions_with_summaries = conditions = read.csv("experiments_to_run.csv", row.names = 1)
exchangeability = umaps = knockoffs = calibration_sim = calibration_gs = list()
i = 0
for(omit_knockout_evaluation_samples in c(T, F)){
  for(condition_index in nrow(conditions):1){
    try({
      conditions$omit_knockout_evaluation_samples = omit_knockout_evaluation_samples
      conditions_with_summaries$omit_knockout_evaluation_samples = omit_knockout_evaluation_samples
      offset = omit_knockout_evaluation_samples*nrow(conditions)
      write.table(conditions[condition_index,], sep = "\t", quote = F)
      attach(conditions[condition_index,], warn.conflicts = F)
      working_directory = ""
      for(j in seq_along(conditions[condition_index,])){
        prefix = paste0( colnames(conditions)[[j]], "=", conditions[condition_index, j] )
        working_directory %<>% file.path(prefix, .)
      }
      # Check out the knockoffs via UMAP and KNN exchangeability test
      # This is more complex if done with knockout samples included & excluded; we just do included
      if(!omit_knockout_evaluation_samples){ 
        withr::with_dir(working_directory, {
          knockoffs[[condition_index + offset]] = read.csv("knockoffs.csv", row.names = 1)
          exchangeability[[condition_index]] = rlookc::KNNTest(X = ecoli_tf_expression, X_k = knockoffs[[condition_index]])
          conditions_with_summaries[condition_index, "KNN_exchangeability_p"]          =  exchangeability[[condition_index]] %>% extract2("p_value")
          conditions_with_summaries[condition_index, "KNN_exchangeability_proportion"] =  exchangeability[[condition_index]] %>% extract2("prop_not_swapped")
          umaps[[condition_index]] = umap::umap(knockoffs[[condition_index]], random.state = 0) %>% extract2("layout")
          umaps[[condition_index]]           %<>% merge(conditions[condition_index,])
          calibration_sim[[condition_index]] = read.csv("average_case_calibration.csv", row.names = 1) %>%
            colMeans %>%
            (function(x) data.frame(nominal_fdr = gsub("^X", "", names(x)) %>% as.numeric, empirical_fdr = x))
          calibration_sim[[condition_index]] %<>% merge(conditions[condition_index,])
        })
      }
      # check vs various gold standards
      for( gold_standard_name in names(ecoli_networks) ){
        try(silent = T,
          {
            i = i + 1
            withr::with_dir(file.path(working_directory, gold_standard_name), {
              calibration_gs[[i]] = read.csv("ecoli_calibration.csv.gz")
              calibration_gs[[i]] %<>% merge(conditions[condition_index,])
              calibration_gs[[i]][["gold_standard_name"]] = gold_standard_name
            })
          })
      }
    })
  }
}
umaps %<>% data.table::rbindlist()
calibration_gs %<>% data.table::rbindlist()
calibration_sim %<>% data.table::rbindlist()
conditions_with_summaries %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )
umaps %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )
calibration_gs %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )
calibration_sim %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )

# Coarsely categorize gold standards for certain plots
calibration_gs %<>% 
  dplyr::mutate(gold_standard_category = "") %>%
  dplyr::mutate(gold_standard_category = ifelse(grepl("chip", gold_standard_name), "RegulonDB chip", gold_standard_category)) %>%
  dplyr::mutate(gold_standard_category = ifelse(grepl("dream5", gold_standard_name), "Exact DREAM5\nreference networks", gold_standard_category)) %>%
  dplyr::mutate(gold_standard_category = ifelse(grepl("10_9", gold_standard_name), "RegulonDB curated", gold_standard_category)) %>%
  dplyr::mutate(gold_standard_category = ifelse(grepl("M3D", gold_standard_name), "M3D knockouts", gold_standard_category)) %>%
  dplyr::mutate(gold_standard_category = ifelse(grepl("regulonDB_knockout", ignore.case = T, gold_standard_name), "RegulonDB knockouts", gold_standard_category)) %>%
  dplyr::mutate(gold_standard_category = ifelse(grepl("negative", ignore.case = T, gold_standard_name), "Biased negative", gold_standard_category)) %>%
  dplyr::mutate(gold_standard_category = ifelse(grepl("positive", ignore.case = T, gold_standard_name), "Biased positive", gold_standard_category)) 
  
dir.create("figures", recursive = T, showWarnings = F)

# Fig <ecoli> A) KNN exchangeability
{
  conditions_with_summaries %>%
    subset(
      !omit_knockout_evaluation_samples &
        condition_on == "none" &
        seed==1 &
        !address_genetic_perturbations
    ) %>%
    dplyr::mutate(proportion_removed = as.character(proportion_removed)) %>%
    dplyr::mutate(knockoff_method = reorder(knockoff_method, KNN_exchangeability_proportion)) %>%
    tidyr::pivot_longer(cols = c("KNN_exchangeability_proportion", "KNN_exchangeability_p"), 
                        names_prefix = "KNN_exchangeability_") %>%
    ggplot() +
    geom_bar(aes(
      x = knockoff_method,
      y = value, 
      fill = proportion_removed
    ), stat = "identity", position = "dodge") +
    facet_wrap(~name, scales = "free_y") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    geom_hline(data = data.frame(name = "proportion", value = 0.5), aes(yintercept = value), color = "red") +
    ggtitle("KNN exchangeability test")
  ggsave("figures/exchangeability.pdf", width = 4, height = 3)
  ggsave("figures/exchangeability.svg", width = 4, height = 3)
}

# Fig <ecoli> B) Simulated target genes
{
  calibration_sim %>%   
    dplyr::mutate(proportion_removed = as.character(proportion_removed)) %>%
    subset(condition_on == "none" & !address_genetic_perturbations & seed==1) %>%
    ggplot() +
    geom_line(aes(x = nominal_fdr, y = empirical_fdr, color = knockoff_method, linetype = proportion_removed)) +
    geom_point(aes(x = nominal_fdr, y = empirical_fdr, color = knockoff_method, shape = proportion_removed)) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    ggtitle("Calibration by method", subtitle = "Simulated target genes") +
    scale_x_continuous(breaks = ((0:2)/2) %>% setNames(c("0", "0.5", "1")), limits = 0:1) +  
    scale_y_continuous(breaks = (0:2)/2, limits = 0:1) +
    coord_fixed() 
  ggsave("figures/simulated_targets.pdf", width = 4, height = 3)
  ggsave("figures/simulated_targets.svg", width = 4, height = 3)
}

# Fig <ecoli> C) Real target genes (various gold standards)
{
  # Vary knockoff generation method (real Y, calibration)
  most_important_gold_standards = c("dream5_new_test", "regulondb10_9", "dream5", "chip_tu_augmented", 
                                    "chip_intersection_M3Dknockout", "chip_intersection_RegulonDB_knockout")
  calibration_gs %>%
    subset(
      !omit_knockout_evaluation_samples &
        condition_on == "none" &
        seed==1 &
        !address_genetic_perturbations &
        nominal_fdr %in% c(0.2, 0.5) & 
        gold_standard_name %in% most_important_gold_standards
      ) %>%
    dplyr::mutate(nominal_fdr = as.character(nominal_fdr)) %>%
    ggplot() +
    geom_pointrange(aes(y = knockoff_method, 
                        xmin = empirical_fdr - moe_95, 
                        xmax = empirical_fdr + moe_95, 
                        x = empirical_fdr,
                        color = nominal_fdr, 
                        shape = proportion_removed), 
                    position = position_dodge(0.5)) +    
    facet_wrap(~wrapFacetLabels( gold_standard_name ) ) +
    geom_vline(aes(xintercept = as.numeric(nominal_fdr), color = nominal_fdr)) +
    ggtitle("Calibration across various gold standards", subtitle = "Faceted by gold standard:") +
    scale_x_continuous(breaks = (0:2)/2, limits = 0:1) 
  ggsave("figures/real_target_genes.pdf", width = 7.5, height = 4)
  ggsave("figures/real_target_genes.svg", width = 7.5, height = 4)
}

# Fig <ecoli> D) confounder adjustment
{
  calibration_gs %>%
    subset(
      !omit_knockout_evaluation_samples &
        knockoff_method == "glasso_0.001" &
        seed==1 &
        !address_genetic_perturbations &
        nominal_fdr %in% c(0.2, 0.5)& 
        gold_standard_name %in% most_important_gold_standards
    ) %>%
    dplyr::mutate(nominal_fdr = as.character(nominal_fdr)) %>%
    ggplot() +
    geom_pointrange(aes(y = condition_on, 
                        xmin = empirical_fdr - moe_95, 
                        xmax = empirical_fdr + moe_95, 
                        x = empirical_fdr,
                        color = nominal_fdr, 
                        shape = do_omit_possible_spouses), 
                    position = position_dodge(0.5)) +    
    facet_wrap(~wrapFacetLabels( gold_standard_name ) )+
    geom_vline(aes(xintercept = as.numeric(nominal_fdr), color = nominal_fdr)) +
    ggtitle("Calibration when correcting for possible confounders") +
    scale_x_continuous(breaks = (0:2)/2, limits = 0:1) 
  ggsave("figures/confounders.svg", width = 7.5, height = 4)
  ggsave("figures/confounders.pdf", width = 7.5, height = 4)
}

# Fig <ecolisupp> A) UMAPs
{
  # Visualize data and knockoffs: does the important structure seem to be captured?
  umap_actual    = umap::umap(ecoli_tf_expression, random.state = 0) %>% extract2("layout")
  umap_actual %<>% merge(data.frame("knockoff_type" = "real data",
                                    "knockoff_method" = "real data",
                                    "shrinkage_param" = NA,
                                    "address_genetic_perturbations" = F,
                                    "omit_knockout_evaluation_samples" = F,
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
    scale_color_discrete(name = "Cluster") +
    ggtitle("UMAP and K-means clusters from E. coli TF expression")
  ggsave(("figures/umaps.pdf"), width = 10, height =6)
  ggsave(("figures/umaps.svg"), width = 10, height =6)
}

# Fig <ecolisupp> B) Gold standard power
{
  # Check how well high-throughput data capture known interactions from regulonDB
  check_chip_power = function( regulon_network, 
                               regulon_network_name = deparse(substitute(regulon_network)),
                               chip_network, 
                               chip_network_name = deparse(substitute(chip_network)) ) {
    regulon_network %<>% subset(is_confirmed) # This function wants everything present to be a positive.
    chip_network %<>% subset(is_confirmed) # This function wants everything present to be a positive.
    regulon_network %>%
      dplyr::mutate(n_expected = Gene1_name %in% chip_network$Gene1_name) %>%
      dplyr::mutate(n_present =
                      Gene1_name %in% chip_network$Gene1_name & 
                      Gene2_name %in% chip_network$Gene2_name) %>% 
      subset(n_expected) %>% 
      with(table(Gene1_name, factor(as.character(n_present), levels = c("FALSE", "TRUE")))) %>% 
      as.data.frame.matrix %>% 
      add_totals %>% 
      set_colnames(c("n_missing", "n_present")) %>%
      dplyr::add_rownames(var = "regulator") %>%
      dplyr::mutate(n_total =   n_missing + n_present ) %>%
      dplyr::mutate(prop_present = n_present/n_total ) %>%
      dplyr::mutate(source_network = regulon_network_name ) %>%
      dplyr::mutate(check_network = chip_network_name) %>%
      dplyr::arrange(n_total) %>% 
      dplyr::mutate(regulator = factor(regulator, levels = regulator))
  }
  chip_power = expand.grid(name1 = names(ecoli_networks), name2 = names(ecoli_networks)) %>%
    subset(name1!=name2) %>%
    with( data = ., 
          mapply(function(name1, name2) {
            check_chip_power(regulon_network = ecoli_networks[[name1]],
                             regulon_network_name = name1,
                             chip_network = ecoli_networks[[name2]],
                             chip_network_name = name2)
          }, 
          name1 = name1, 
          name2 = name2,
          SIMPLIFY = F)
    ) %>%
    Reduce(rbind, .)
  
  chip_power %>%  
    subset(source_network %in%  c("dream5_new_test", "regulondb10_9", "chip_tu_augmented", "regulonDB_knockout_tu_augmented", "M3Dknockout_tu_augmented")) %>%
    subset( check_network %in%  c("chip_tu_augmented", "regulonDB_knockout_tu_augmented", "M3Dknockout_tu_augmented")) %>%
    ggplot() + 
    geom_bar(stat = "identity", aes(x = regulator, y = -n_missing), fill = "red") + 
    geom_bar(stat = "identity", aes(x = regulator, y =  n_present), fill ="black") + 
    coord_flip() + 
    ggtitle("Gold standard power to verify prior findings") + 
    facet_grid( paste0("Edges compared to:\n", check_network) ~ paste0("Edges taken from:\n", source_network), scales = "free") +
    ylab("Number of hypotheses that are \n(not confirmed | confirmed)") 
  ggsave("figures/gs_comparison.pdf", width = 10, height = 7)
  ggsave("figures/gs_comparison.svg", width = 10, height = 7)
}

# Fig. <ecolisupp> Study of bias in gold standards
{
  # Vary knockoff generation method (real Y, calibration)
  most_important_gold_standards = c("chip_intersection_M3Dknockout"                        ,
                                    "chip_intersection_RegulonDB_knockout"                ,
                                    "chip_intersection_M3Dknockout_biased_positive"        ,
                                    "chip_intersection_RegulonDB_knockout_biased_positive",
                                    "chip_intersection_M3Dknockout_biased_negative"       ,
                                    "chip_intersection_RegulonDB_knockout_biased_negative" )
  calibration_gs %>%
    subset(
      !omit_knockout_evaluation_samples &
        seed==1 &
        !address_genetic_perturbations &
        nominal_fdr %in% c(0.2, 0.5) & 
        knockoff_method == "glasso_0.001" &
        gold_standard_name %in% most_important_gold_standards
    ) %>%
    dplyr::mutate(nominal_fdr = as.character(nominal_fdr)) %>%
    ggplot() +
    geom_pointrange(aes(y = condition_on, 
                        xmin = empirical_fdr - moe_95, 
                        xmax = empirical_fdr + moe_95, 
                        x = empirical_fdr,
                        color = nominal_fdr, 
                        shape = do_omit_possible_spouses), 
                    position = position_dodge(0.5)) +
    facet_wrap(~wrapFacetLabels(gold_standard_name) ) +
    geom_vline(aes(xintercept = as.numeric(nominal_fdr), color = nominal_fdr)) +
    ggtitle("Effect of gold standard bias") +
    scale_x_continuous(breaks = (0:2)/2, limits = 0:1) 
  ggsave("figures/biased_gold_standards.pdf", width = 7.5, height = 4)
  ggsave("figures/biased_gold_standards.svg", width = 7.5, height = 4)
}


  
  
  
  
  