# setwd("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v34")
source("../dream5_ecoli_gold_standards.R")
source("../dream5_ecoli_setup.R")
fillna =function(x, filler){ x[is.na(x)] = filler; x}
wrapFacetLabels = function(x){
  for(i in seq_along(x)){
    x[[i]] = x[[i]] %>%
      gsub( "_", " ", . ) %>%
      strwrap( width = 20 ) %>%
      paste0( collapse = "\n" )
  }
  return(x)
}
ggplot2::theme_update(text = element_text(family = "ArialMT"))
# Collect results to make certain key summary plots
# This will allow us to show calibration on any combination of gold standard, knockoff construction method,
# number of principal components, et cetera by subsetting a giant tidy dataframe. 
cat("Loading experimental results...")
conditions_with_summaries = conditions = read.csv("experiments_to_run.csv", row.names = 1)
exchangeability = tsnes = knockoffs = calibration_sim = calibration_gs = list()
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
    # Check out the knockoffs via tsne and KNN exchangeability test
    # This is more complex if done with knockout samples included & excluded; we just do included
    if(!omit_knockout_evaluation_samples){ 
      withr::with_dir(working_directory, {
        if(conditions_with_summaries[condition_index, "fdr_control_method"]=="knockoffs") {
          conditions_with_summaries[condition_index, "index_in_list_of_knockoffs"] = paste0(condition_index, "_")
          knockoffs[[paste0(condition_index, "_")]] = read.csv("knockoffs.csv", row.names = 1) %>% set_colnames(NULL)
          exchangeability[[condition_index]] = rlookc::KNNTest(X = ecoli_tf_expression, X_k = knockoffs[[condition_index]])
          conditions_with_summaries[condition_index, "KNN_exchangeability_p"]          =  exchangeability[[condition_index]] %>% extract2("p_value")
          conditions_with_summaries[condition_index, "KNN_exchangeability_proportion"] =  exchangeability[[condition_index]] %>% extract2("prop_not_swapped")
        } else {
          conditions_with_summaries[condition_index, "KNN_exchangeability_p"]          =  NA
          conditions_with_summaries[condition_index, "KNN_exchangeability_proportion"] =  NA
        }
        # Wrapped in trycatch because sometimes we don't do these simulations and there's no data to slurp in. 
        calibration_sim[[condition_index]] = tryCatch(
          {
            print(getwd())
            read.csv("average_case_calibration.csv", row.names = 1) %>%
              colMeans %>%
              (function(x) data.frame(expected_fdr = gsub("^X", "", names(x)) %>% as.numeric, observed_fdr = x)) %>% 
              merge(conditions[condition_index,])
          },
          error = function(e) data.frame()
        )
      })
    }
    # check vs various gold standards
    for( gold_standard_name in names(ecoli_networks) ){
      try(
        silent = T,
        {
          i = i + 1
          withr::with_dir(file.path(working_directory, gold_standard_name), {
            calibration_gs[[i]] = read.csv("ecoli_calibration.csv.gz")
            calibration_gs[[i]] %<>% merge(conditions[condition_index,])
            calibration_gs[[i]][["gold_standard_name"]] = gold_standard_name
          }
          )
        }
      )
    }
  })
}

calibration_gs %<>% data.table::rbindlist()
calibration_sim %<>% data.table::rbindlist()
conditions_with_summaries %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )
calibration_gs %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )
calibration_sim %<>% mutate( knockoff_method = paste0(knockoff_type, "_", shrinkage_param) %>% gsub("_NA", "", .) )

# Rename certain things to customize plots for our audience
calibration_gs %<>% rename(expected_fdr = nominal_fdr, observed_fdr = empirical_fdr)
prettify_names = function(df){
  df %<>% mutate(knockoff_type = gsub("naive", "permuted", knockoff_type))
  try({
    df %<>% mutate(knockoff_method = gsub("naive", "permuted", knockoff_method))
  }, silent = T)
  df %<>% mutate(condition_on = gsub("pert_labels", "Protocol annotations", condition_on))
  df %<>% mutate(condition_on = gsub("plus_pca", "and top PC's (", condition_on))
  df %<>% mutate(condition_on = gsub("_", " ", condition_on))
  df %<>% mutate(condition_on = gsub("0", "0)", condition_on))
  df %<>% mutate(condition_on = gsub("0))", "0)", condition_on))
  df %<>% dplyr::mutate(gold_standard_name = gsub("_intersection_", " and ", gold_standard_name)) 
  df %<>% dplyr::mutate(gold_standard_name = gsub("_tu_augmented|_", " ", gold_standard_name)) 
  df %<>% dplyr::mutate(gold_standard_name = gsub("new test", "validation", gold_standard_name)) 
  return(df)
}
conditions %<>% prettify_names
conditions_with_summaries %<>% prettify_names
calibration_sim %<>% prettify_names
calibration_gs %<>% prettify_names

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
    dplyr::mutate(knockoff_method = reorder(knockoff_method, KNN_exchangeability_proportion)) %>%
    tidyr::pivot_longer(cols = c("KNN_exchangeability_proportion", "KNN_exchangeability_p"), 
                        names_prefix = "KNN_exchangeability_") %>%
    ggplot() +
    geom_point(aes(
      x = knockoff_method,
      y = value, 
    )) +
    facet_wrap(~name) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    geom_hline(data = data.frame(name = "proportion", value = 0.5), aes(yintercept = value), color = "red") +
    ggtitle("KNN exchangeability test", "Real data")
  ggsave("figures/exchangeability.pdf", width = 5, height = 3)
  ggsave("figures/exchangeability.svg", width = 5, height = 3)
}

# Fig <ecoli> B) Simulated target genes
{
  calibration_sim %>%   
    subset(T & 
             condition_on == "none" & 
             !address_genetic_perturbations & 
             seed==1 
    ) %>%
    subset(do_simulate_y) %>%
    subset(!is.na(knockoff_method)) %>%
    ggplot(
      aes(
        x = expected_fdr, 
        y = observed_fdr, 
        color = paste(
          fdr_control_method,
          gsub("NA", "", knockoff_method)
        )
      )
    ) +
    labs(colour = "FDR control method") + 
    geom_line() +
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0)) +
    ggtitle("Calibration by method", subtitle = "Simulated target genes") +
    scale_x_continuous(breaks = ((0:2)/2) %>% setNames(c("0", "0.5", "1")), limits = 0:1) +  
    scale_y_continuous(breaks = (0:2)/2, limits = 0:1) +
    coord_fixed() 
  ggsave("figures/simulated_targets.pdf", width = 5, height = 4)
  ggsave("figures/simulated_targets.svg", width = 5, height = 4)
}

# Fig <ecoli> C) Real target genes (various gold standards)
{
  # Vary knockoff generation method (real Y, calibration)
  most_important_gold_standards = c(#"dream5_new_test", "regulondb10_9", "dream5", "chip_tu_augmented", 
    "chip and M3Dknockout", "chip and RegulonDB knockout")
  calibration_gs %>%
    subset(
      !omit_knockout_evaluation_samples &
        condition_on == "none" &
        do_simulate_y &
        seed==1 &
        !address_genetic_perturbations &
        gold_standard_name %in% most_important_gold_standards & 
        do_omit_possible_spouses 
    ) %>%
    ggplot(
      aes(x = expected_fdr, 
          y = observed_fdr, 
          color = paste(
            fdr_control_method,
            gsub("NA", "", knockoff_method)
          )
      )
    ) +
    labs(colour = "FDR control method") + 
    geom_point() + 
    geom_line() +
    geom_abline(intercept=0, slope=1) +
    facet_wrap(~wrapFacetLabels( gold_standard_name ), scales = "free_y" ) +
    ggtitle("Calibration in gold standards based on \nChIP and perturbation transcriptomics", "Real data") +
    scale_x_continuous(breaks = (0:2)/2, limits = 0:1)  
  ggsave("figures/real_target_genes.pdf", width = 6.5, height = 3.5)
  ggsave("figures/real_target_genes.svg", width = 6.5, height = 3.5)
}

# Fig <ecoli> D) confounder adjustment
{
  calibration_gs %>%
    subset(
      !omit_knockout_evaluation_samples &
        knockoff_method == "glasso_1e-04" &
        seed==1 &
        !address_genetic_perturbations &
        gold_standard_name %in% most_important_gold_standards
    ) %>%
    # dplyr::mutate(expected_fdr = as.character(expected_fdr)) %>%
    ggplot(mapping = aes(x = expected_fdr, 
                         y = observed_fdr,
                         color = condition_on, 
                         shape = condition_on)) +   
    geom_point() + 
    geom_line() +
    geom_abline(yintercept=0, slope=1) +
    facet_wrap(~wrapFacetLabels( gold_standard_name ), scales = "free_y" )+
    ggtitle("Calibration when correcting for possible confounders", "Real data") +
    scale_x_continuous(breaks = (0:2)/2, limits = 0:1) 
  ggsave("figures/confounders.svg", width = 9, height = 3)
  ggsave("figures/confounders.pdf", width = 9, height = 3)
}

# Fig <ecoli> F) confounder adjustment power
{
  calibration_gs %>%
    subset(
      !omit_knockout_evaluation_samples &
        knockoff_method == "glasso_1e-04" &
        seed==1 &
        !address_genetic_perturbations &
        gold_standard_name %in% most_important_gold_standards
    ) %>%
    ggplot(mapping = aes(x = expected_fdr, 
                         y = log10(num_discoveries+1),
                         color = condition_on, 
                         shape = condition_on)) +   
    geom_point() + 
    facet_wrap(~wrapFacetLabels( gold_standard_name ) )+
    ggtitle("Power when correcting for possible confounders", "Real data") +
    scale_x_continuous(breaks = (0:2)/2, limits = 0:1) + 
    ylab("Log10( 1 + number of discoveries )")
  ggsave("figures/confounders_power.svg", width = 7, height = 2.5)
  ggsave("figures/confounders_power.pdf", width = 7, height = 2.5)
}

# Fig <ecolisupp> A) tsnes
{
  # Visualize data and knockoffs: does the important structure seem to be captured?
  which_to_tsne = conditions_with_summaries %>% 
    subset(
      T &
        knockoff_method != "NA" & 
        seed==1 & 
        condition_on == "none"  & 
        omit_knockout_evaluation_samples == F & 
        do_simulate_y
    )
  matrices_to_tsne = data.table::rbindlist( c( 
    knockoffs[which_to_tsne$index_in_list_of_knockoffs], 
    list( data.frame( ecoli_tf_expression ) )
  ) ) 
  which_to_tsne %<>% rbind( ., "real data" )
  set.seed(0) # make tsne repeatable
  tsnes = tsne::tsne(matrices_to_tsne)
  facet_order = c(
    "real data",
    "sample",
    "glasso_1e-04",
    "glasso_0.001",
    "glasso_0.01",
    "shrinkage", 
    "mixture",
    "permuted"
  )
  tsnes %>% 
    as.data.frame() %>%
    set_colnames(c("V1", "V2")) %>% 
    mutate(knockoff_method = factor(which_to_tsne$knockoff_method %>% rep( each = 805 ), levels = facet_order)) %>%
    mutate(cluster = clusters$cluster %>% rep(times = dplyr::n()/805) %>% as.character()) %>%
    subset(condition_on == "none" & !address_genetic_perturbations & seed==1) %>%
    ggplot() +
    geom_point(aes(x = V1,  y = V2, color = cluster)) +
    coord_fixed() +
    xlab("tsne1") + ylab("tsne2") +
    facet_wrap(~knockoff_method, nrow = 2) +
    ggtitle("Joint embedding of E. coli TF expression and knockoffs", "Real data")
  ggsave(("figures/tsnes.pdf"), width = 10, height =6)
  ggsave(("figures/tsnes.svg"), width = 10, height =6)
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
      tibble::rownames_to_column(var = "regulator") %>%
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
    dplyr::mutate(source_network = gsub("_tu_augmented|_", " ", source_network)) %>%
    dplyr::mutate(check_network = gsub("_tu_augmented|_", " ", check_network)) %>%
    dplyr::mutate(source_network = gsub("new test", "validation", source_network)) %>%
    dplyr::mutate(check_network = gsub("new test", "validation", check_network)) %>%
    ggplot() + 
    geom_bar(stat = "identity", aes(x = regulator, y = -n_missing), fill = "red") + 
    geom_bar(stat = "identity", aes(x = regulator, y =  n_present), fill ="black") + 
    coord_flip() + 
    ggtitle("Gold standard power", subtitle = "Edges from ... ") + 
    facet_grid( check_network ~ source_network, scales = "free", switch = "y") +
    ylab("Number of hypotheses that are \n(not confirmed | confirmed)") + 
    xlab("Edges compared to ...") 
  
  ggsave("figures/gs_comparison.pdf", width = 10, height = 7)
  ggsave("figures/gs_comparison.svg", width = 10, height = 7)
}









