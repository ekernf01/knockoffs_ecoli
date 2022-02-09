# setwd("~/Desktop/jhu/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v26/")
suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("magrittr")
  library("glasso")
  library("limma")
  
})
set.seed(0)
sink("sessionInfo.txt")
sessionInfo()
sink()

# for eric laptop
DATALAKE = "~/Desktop/jhu/research/datalake"
# for AWS
if(!dir.exists(DATALAKE)){
  DATALAKE = "~/datalake"
}
# otherwise, complain.
if(!dir.exists(DATALAKE)){
  stop("Datalake not found. Place it in '~/datalake' or Modify `dream5_ecoli_setup.R`.\n")
}


# Set to T to run a fast test of the code. Results will not be biologically 
# meaningful; this is for software testing. 
test_mode = F

# We'll run several experiments with different handling of perturbations,
# confounding influences, and knockoff construction.
conditions = rbind(
  # study random variation in knockoff construction
  expand.grid(
    knockoff_type = c("sample"),
    address_genetic_perturbations = c( F ),
    condition_on = c( "none" ),
    shrinkage_param = NA,
    seed = 2:5
  ),
  # study knockoff generation methods
  expand.grid(
    knockoff_type = c("naive", "shrinkage", "mixture"),
    address_genetic_perturbations = c( F ),
    condition_on = c( "none"),
    shrinkage_param = NA,
    seed = 1
  ),
  # study other knockoff generation methods
  expand.grid(
    knockoff_type = c( "glasso"),
    address_genetic_perturbations = c( F ),
    condition_on = c( "none"),
    shrinkage_param = c(0.1, 0.01, 0.001),
    seed = 1
  ),
  # study confounder handling
  expand.grid(
    knockoff_type = c("sample"),
    address_genetic_perturbations = c( F, T ),
    condition_on = c(  "none",
                       "pert_labels",
                       "pert_labels_plus_pca10",
                       "pert_labels_plus_pca20",
                       "pert_labels_plus_pca30",
                       "pert_labels_plus_pca50" ),
    shrinkage_param = NA,
    seed = 1
  ),
  # study confounder handling
  expand.grid(
    knockoff_type = c("glasso"),
    address_genetic_perturbations = c( F, T ),
    condition_on = c(  "none",
                       "pert_labels",
                       "pert_labels_plus_pca10",
                       "pert_labels_plus_pca20",
                       "pert_labels_plus_pca30",
                       "pert_labels_plus_pca50" ),
    shrinkage_param = 0.001,
    seed = 1
  )
)
write.csv(conditions, "experiments_to_run.csv")


if(test_mode){
  # return random results for fast testing of the code
  fast_lasso_penalty = function(X, X_k, y) {
    cat(".")
    rnorm(ncol(X))
  }
} else {
  # Compute just part of the lasso path for faster results
  dfmax = 21
  fast_lasso_penalty = function(X, X_k, y) {
    cat(".")
    suppressWarnings(
      knockoff::stat.lasso_lambdasmax(
        y   = y,
        X   = X,
        X_k = X_k,
        dfmax = dfmax
      )
    )
  }
}

# Load data
cat("\nPrepping data.\n")
{
  withr::with_dir(
    file.path(DATALAKE, "dream5/DREAM5_network_inference_challenge/Network3"),
    {
      ecoli_expression      = read.table("input data/net3_expression_data.tsv", header = T)
      ecoli_metadata        = read.table("input data/net3_chip_features.tsv", header = T, comment.char = "")
      ecoli_tf              = read.table("input data/net3_transcription_factors.tsv")
      ecoli_network_dream5  = read.table("gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv")
      ecoli_anonymization   = read.table("anonymization/net3_gene_ids.tsv")
    }
  )
  ecoli_dream5_new_test = read.csv("../DREAM5_ecoli_53_novel_hypotheses.csv", header = T, comment.char = "#")
  
  # Make a couple useful alterations to metadata
  ecoli_metadata$Time_Rank = ave(ecoli_metadata$Time, ecoli_metadata$X.Experiment,
                                 FUN = function(x) rank(x, na.last = "keep")/max(rank(x)) )
  ecoli_metadata$Perturbations[is.na(ecoli_metadata$Perturbations)] = "none"

  # De-anonymize genes
  ecoli_anonymization_by_name = setNames(toupper(ecoli_anonymization$V2), ecoli_anonymization$V1)
  ecoli_network_dream5 %<>% set_colnames(c("Gene1", "Gene2", "is_confirmed"))
  ecoli_network_dream5 %<>% subset(is_confirmed>0)
  ecoli_network_dream5 %<>% dplyr::mutate(Gene1_name = ecoli_anonymization_by_name[Gene1] )
  ecoli_network_dream5 %<>% dplyr::mutate(Gene2_name = ecoli_anonymization_by_name[Gene2] )
  ecoli_network_dream5$Gene1 = NULL
  ecoli_network_dream5$Gene2 = NULL
  
  # Don't worry about capitalization
  ecoli_dream5_new_test$Gene1_name %<>% toupper
  ecoli_dream5_new_test$Gene2_name %<>% toupper
}

cat("\nPrepping knockout-based gold standard.\n")
{
  ecoli_metadata = ecoli_metadata %>% 
    dplyr::group_by(X.Experiment) %>% 
    dplyr::summarise(
      HasControl = any( is.na(DeletedGenes)), 
      HasKO = any(!is.na(DeletedGenes)),
      HasReplication = max(Repeat)>1, 
      HasNoKnownConfounding = length(unique(Perturbations))==1 
    ) %>%
    merge(ecoli_metadata, ., all.x = T) %>% 
    dplyr::mutate(HasKO = !is.na(HasKO) & HasKO) %>% 
    dplyr::mutate(HasControl = !is.na(HasControl) & HasControl) %>% 
    dplyr::mutate(DeletedGenesNoNA = ifelse(is.na(DeletedGenes), "", DeletedGenes)) 
  ecoli_network_ko = list()
  for(e in unique(subset(ecoli_metadata, HasKO & HasControl & HasReplication & HasNoKnownConfounding, select = "X.Experiment", drop = T))){
    current_expression =    t(ecoli_expression[ecoli_metadata$X.Experiment==e,])
    current_metadata =
      ecoli_metadata %>%
      subset(X.Experiment==e) %>%
      dplyr::mutate(DeletedGenesNoNA = factor(DeletedGenesNoNA, levels = sort(unique(DeletedGenesNoNA))) )
    design.trt=model.matrix(~DeletedGenesNoNA, data = current_metadata)
    fitTrtMean <- lmFit(current_expression, design.trt)
    fit.contrast=contrasts.fit(fitTrtMean, coefficients = colnames(design.trt)[-1])
    efit.contrast=eBayes(fit.contrast)
    ecoli_network_ko[[as.character(e)]] = 
      efit.contrast$p.value %>%
      as.data.frame() %>%
      dplyr::add_rownames(var = "Gene2_name") %>%
      tidyr::pivot_longer(cols = !Gene2_name, 
                          names_prefix = "DeletedGenesNoNA",
                          names_to = "Gene1_name", 
                          values_to = "p_value_KO") %>%
      mutate(X.Experiment = e)
  }
  ecoli_network_ko %<>% data.table::rbindlist()
  ecoli_network_ko = ecoli_network_ko[,c("Gene1_name", "Gene2_name", "p_value_KO", "X.Experiment")]
  ecoli_network_ko %<>% subset(!grepl(",", Gene1_name)) #Ignore double knockouts for now
  ecoli_network_ko$q_value_KO = p.adjust(ecoli_network_ko$p_value_KO, method = "fdr")
  ecoli_network_ko$Gene1_name = ecoli_anonymization_by_name[ecoli_network_ko$Gene1_name] %>% toupper
  ecoli_network_ko$Gene2_name = ecoli_anonymization_by_name[ecoli_network_ko$Gene2_name] %>% toupper
  ecoli_network_ko %>% 
    mutate(targetIsDecoy = grepl("DECOY", Gene2_name)) %>%
    ggplot()  +
    geom_histogram(aes(x = p_value_KO, fill = targetIsDecoy), bins = 100)  + 
    facet_wrap(~targetIsDecoy, ncol = 1, scales = "free_y") +
    ggtitle( "Decoy gene differential expression")
}

cat("\nPrepping ChIP gold standard.\n")
# Load curated interactions from regulonDB v10.9
{
  tidy_regulondb = function(ecoli_network_regulondb10_9) {
    ecoli_network_regulondb10_9 %<>%
      subset(V8=="gene" & V7=="tf") %>%
      extract(c(2, 4)) %>%
      set_colnames(c("Gene1_name", "Gene2_name"))
    ecoli_network_regulondb10_9[["Gene1_name"]] %<>% toupper
    ecoli_network_regulondb10_9[["Gene2_name"]] %<>% toupper
    ecoli_network_regulondb10_9
  }
  ecoli_network_regulondb10_9 =
    read.table(
      sep = "\t",
      file.path(DATALAKE, "modern_ecoli/regulondb_10.9/full/genetic_network.txt")
    ) 
  ecoli_network_regulondb10_9_additional_binding_support =
    ecoli_network_regulondb10_9 %>%
    subset(grepl("sequence|binding|SELEX|mutation", ignore.case = T, V6)) %>% # Evidence of direct binding
    tidy_regulondb
  ecoli_network_regulondb10_9 %<>% tidy_regulondb
}

# Load ChIP-seq data from regulonDB v10.9
{
  schema = read.table(sep = ":", strip.white = T,comment.char = "", row.names = 1, text =
                        "#Column 1: RegulonDB Dataset Identifier
    #Column 2: Dataset Type
    #Column 3: Dataset Name
    #Column 4: Transcription Factor Name
    #Column 5: Effect
    #Column 6: Regulated Object
    #Column 7: TFBs Absolute genome left position
    #Column 8: TFBs Absolute genome right position
    #Column 9: DNA sequence
    #Column 10: Growth Condition
    #Column 11: Reference PMID")$V2
  ecoli_network_chip =
    list.files(
      file.path(DATALAKE,"modern_ecoli/regulondb_10.9/highthroughput/"),
      full.names = T
    ) %>%
    lapply(read.table, sep = "\t") %>%
    data.table::rbindlist() %>%
    set_colnames(schema) %>%
    as.data.frame()
  # MelR regulates MelAB. Missing from RegulonDB but there is high quality data from an included study.
  # See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1308901/
  ecoli_network_chip %<>% rbind(read.csv(row.names = 1, check.names = F, text = 
    "RegulonDB Dataset Identifier,Dataset Type,Dataset Name,Transcription Factor Name,Effect,Regulated Object,TFBs Absolute genome left position,TFBs Absolute genome right position,DNA sequence,Growth Condition,Reference PMID
    1,<NA>,TF CHIP,<NA>,MelR,<NA>,MelA,<NA>,<NA>,<NA>,<NA>,16301522
    2,<NA>,TF CHIP,<NA>,MelR,<NA>,MelB,<NA>,<NA>,<NA>,<NA>,16301522"
  ))
  # Filter out some samples
  ecoli_network_chip_metadata = read.csv(file.path(DATALAKE,"modern_ecoli/regulondb_10.9/table_chip_controls.csv"))
  ecoli_network_chip %<>% merge(ecoli_network_chip_metadata, by.x = "Reference PMID", by.y = "PMID")
  ecoli_network_chip %<>% subset( 
    "keep" == Recommendation |  # keep whole studies
      mapply(grepl, `Transcription Factor Name`, Recommendation, ignore.case = T) # keep only certain regulators from a couple studies
  )
  # Tidy to match other gold standards
  ecoli_network_chip %<>%
    extract(c("Transcription Factor Name", "Regulated Object")) %>%
    set_colnames(c("Gene1_name", "Gene2_name")) %>%
    dplyr::mutate(is_confirmed = T)
  ecoli_network_chip = ecoli_network_chip[!duplicated(ecoli_network_chip),]

  # Fix some genes with multiple names
  ecoli_network_chip$Gene1_name[ecoli_network_chip$Gene1_name=="Sigma70"] = "rpoD"
  ecoli_network_chip$Gene1_name[ecoli_network_chip$Gene1_name=="Sigma38"] = "rpoS"
  ecoli_network_chip$Gene1_name[ecoli_network_chip$Gene1_name=="Sigma32"] = "rpoH"
  ecoli_network_chip$Gene1_name[ecoli_network_chip$Gene1_name=="Cra"    ] = "fruR"
  ecoli_network_chip$Gene1_name[ecoli_network_chip$Gene1_name=="H-NS"   ] = "hns"
  ecoli_network_chip$Gene1_name[ecoli_network_chip$Gene1_name=="IHF"    ] = "ihfA"
  # Fix some annotations like "WITHIN bcsE" to "bcsE"
  ecoli_network_chip$Gene2_name %<>% gsub("WITHIN |END OF |PSEUDOGENE|-[0-9]", "", ignore.case = T, .)
  # Fix multi-gene entries
  is_ecoli_operon = function(operon){
    is_lower = function(x) x==tolower(x)
    is_upper = function(x) x==toupper(x)
    nchar(operon)>4 && is_lower(substr(operon, 1, 3)) && is_upper(substr(operon, start = 4, stop = nchar(operon)))
  } 
  separate_ecoli_operon = function(operon){
    if(is_ecoli_operon(operon)){
      prefix = substr(operon, 1, 3)
      suffixes = substr(operon, start = 4, stop = nchar(operon)) %>% strsplit(split = "") %>% extract2(1)
      gene_list = paste(paste0(prefix, suffixes), collapse = ",")
      return(gene_list)
    } else {
      return(operon)
    }
  }
  ecoli_network_chip %<>% dplyr::mutate(Gene2_name = sapply(Gene2_name, separate_ecoli_operon))
  ecoli_network_chip %<>% tidyr::separate_rows(Gene2_name) 
  # Remove errant suffixes like -1 and -2
  ecoli_network_chip %<>% subset(nchar(Gene2_name)>=3)
  ecoli_network_chip =
    ecoli_network_chip %>%
    subset(Gene1_name=="ihfA") %>%
    mutate(Gene1_name="ihfB") %>%
    rbind(ecoli_network_chip)
  ecoli_network_chip[["Gene1_name"]] %<>% toupper
  ecoli_network_chip[["Gene2_name"]] %<>% toupper
}

# Augment the ChIP benchmark with all genes of any targeted transcription unit.
{
  ecoli_tu = readr::read_delim(
    file.path(DATALAKE, "modern_ecoli/transcription_units/All_transcription_units_of_E._coli_K-12_substr._MG1655.txt"),
    delim = "\t", col_types = readr::cols()
  )
  colnames(ecoli_tu)[[2]] = "Gene2_name"
  ecoli_tu[["Gene2_name"]] %<>% gsub(pattern = "3'ETS-<i>leuZ</i>", replacement = "3ETSleuZ", fixed = T)
  ecoli_tu %<>% tidyr::separate_rows(Gene2_name, )
  ecoli_tu[["Gene2_name"]] %<>% toupper
  # Tricky pair of merges here.
  # 1. Add the transcription unit for each target gene
  # 2. Add the other target genes for each transcription unit
  ecoli_network_tu_augmented = merge(ecoli_network_chip, ecoli_tu, by = "Gene2_name", all.x = T)
  ecoli_network_tu_augmented = merge(ecoli_network_tu_augmented, ecoli_tu, by = "Transcription-Units", all.y = T)
  ecoli_network_tu_augmented %<>% dplyr::mutate(not_in_original = Gene2_name.x != Gene2_name.y)
  ecoli_network_tu_augmented[["Gene2_name"]] = ecoli_network_tu_augmented[["Gene2_name.y"]]
  ecoli_network_tu_augmented[["Gene2_name.x"]] = NULL
  ecoli_network_tu_augmented[["Gene2_name.y"]] = NULL
  doops = duplicated(ecoli_network_tu_augmented[c("Gene1_name","Gene2_name")])
  ecoli_network_tu_augmented = ecoli_network_tu_augmented[!doops, ]
  ecoli_network_tu_augmented %<>% subset(!is.na(Gene1_name))
  #Gold standard sizes with and without augmentation:
  dim(ecoli_network_chip)
  dim(ecoli_network_tu_augmented)
  table(ecoli_network_chip$Gene1_name) %>% sort
  table(ecoli_network_tu_augmented$Gene1_name) %>% sort
}
AVAILABLE_GOLD_STANDARDS = c( "dream5", "curated", "chip", "chip_augmented", "knockout" )

# Check out the expression data briefly
ecoli_expression[1:4, 1:4]
dim(ecoli_expression)
ecoli_tf_expression = ecoli_expression[ecoli_tf[[1]]] %>% as.matrix
clusters = kmeans(ecoli_expression, centers = 7)

# Look at genetic perturbations' effect on expression
try({
  perturbations = data.frame(gene = NA, direction = NA, mean_control = 0, mean_perturbed = 0)
  i = 0
  for(g in unique(ecoli_metadata$OverexpressedGenes)){
    if(is.na(g)){
      next
    }
    i = i + 1
    idx = ecoli_metadata$OverexpressedGenes==g
    idx[is.na(idx)] = F
    perturbations[i, "gene"] = g
    perturbations[i, "direction"] = "overexpression"
    perturbations[i, "mean_control"]   = mean(ecoli_expression[!idx, g])
    perturbations[i, "mean_perturbed"] = mean(ecoli_expression[ idx, g])
  }
  for(g in unique(ecoli_metadata$DeletedGenes)){
    if(is.na(g)){
      next
    }
    i = i + 1
    idx = ecoli_metadata$DeletedGenes==g
    idx[is.na(idx)] = F
    perturbations[i, "gene"] = g
    perturbations[i, "direction"] = "deletion"
    g = strsplit(g, ",")[[1]]
    perturbations[i, "mean_control"]   = mean( as.matrix(ecoli_expression)[!idx, g])
    perturbations[i, "mean_perturbed"] = mean( as.matrix(ecoli_expression)[ idx, g])
  }
  ggplot(perturbations) +
    geom_point(aes(x = mean_control, y = mean_perturbed, colour = direction)) +
    geom_abline(aes(slope = 1, intercept = 0)) +
    ggtitle("Genetic perturbations in E coli")
  ggsave("genetic_perturbations.pdf", width = 6, height = 6)
})

# This will be used later to check the results.
check_against_gold_standards = function(DF){
  for( gold_standard_name in AVAILABLE_GOLD_STANDARDS ){
    if(gold_standard_name=="dream5"){
      gold_standard = ecoli_network_dream5
    } else if(gold_standard_name=="knockout"){
      gold_standard = ecoli_network_ko
    } else if(gold_standard_name=="curated"){
      gold_standard = ecoli_network_regulondb10_9
    } else if(gold_standard_name=="chip"){
      gold_standard = ecoli_network_chip
    } else if(gold_standard_name=="chip_augmented"){
      gold_standard = ecoli_network_tu_augmented
    } else {
      stop("Gold standard by that name not found.\n")
    }
    gold_standard = gold_standard[c("Gene1_name", "Gene2_name")]
    gold_standard[["is_confirmed"]] = T
    dir.create(gold_standard_name, recursive = T, showWarnings = F)
    withr::with_dir(gold_standard_name, {
      # Show degree distribution of gold standard
      gold_standard %>%
        extract2("Gene1_name") %>%
        table %>%
        hist(40, main = "TF promiscuity", xlab = "Number of targets per TF")
      gold_standard %>%
        extract2("Gene2_name") %>%
        table %>%
        hist(main = "Target promiscuity", xlab = "Number of TFs per target")
    })
    # Include forwards & backwards edges -- we don't expect to get the directionality right here.
    gold_standard_reversed = gold_standard
    gold_standard_reversed[1:2] = gold_standard[2:1]
    gold_standard_symmetric = rbind(gold_standard, gold_standard_reversed)
    # If forwards is true and backwards is false, keep only the trues in both directions.
    gold_standard_symmetric %<>% dplyr::arrange(desc(is_confirmed))
    gold_standard_symmetric = gold_standard_symmetric[!duplicated(gold_standard_symmetric$Gene1,
                                                                  gold_standard_symmetric$Gene2),]
    # This makes true edges correct; the rest is NA and fails to distinguish between
    # unknowns and known negatives.
    DF_plus_gs = merge(DF, gold_standard_symmetric, by = c("Gene1_name", "Gene2_name"), all.x = T, all.y = F)
    # This treats everything as a known negative.
    DF_plus_gs[["is_confirmed"]][is.na(DF_plus_gs[["is_confirmed"]])] = F
    # In ChIP data, we can distinguish between unknowns (not chipped) and known
    # negatives (chipped and not bound).
    if(grepl("chip", gold_standard_name)){
      DF_plus_gs[["is_confirmed"]][ !( DF_plus_gs[["Gene1_name"]] %in% gold_standard[["Gene1_name"]] ) ] = NA
    }
    # Omit edges where that target gene is NEVER present in the gold standard.
    # This could be viewed as overly generous to us but it makes sense if there are naming inconsistencies.
    DF_plus_gs[["is_confirmed"]][ !( DF_plus_gs[["Gene2_name"]] %in% gold_standard[["Gene2_name"]] ) ] = NA
    DF_plus_gs[["is_confirmed"]] %>% table
    DF_plus_gs[["is_confirmed"]] %>% is.na %>% table
    dim(DF_plus_gs)

    # Compute calibration vs gold standard
    {
      calibration = data.frame(nominal_fdr = (c(1:100)/100) %>% c(quantile(DF$q, probs = c(1:100)/100)),
                               empirical_fdr = NA,
                               num_discoveries = NA)
      i = 0
      for(fdr in calibration$nominal_fdr){
        i = i + 1
        calibration$empirical_fdr[i] =
          DF_plus_gs %>%
          subset(q<fdr) %>%
          extract2("is_confirmed") %>%
          mean(na.rm = T) %>%
          subtract(1, .)
        calibration$num_discoveries[i] =
          DF_plus_gs %>%
          subset(q<fdr) %>%
          extract2("is_confirmed") %>%
          is.na %>%
          not %>%
          sum
      }
      # Add simple binomial standard errors
      calibration %<>% mutate( moe_95 = 1.96 * sqrt( empirical_fdr * ( 1 - empirical_fdr ) / num_discoveries ) )
    }
    withr::with_dir(gold_standard_name, {
      # Save all results for follow-up
      write.csv(DF_plus_gs, ("results_with_evaluation.csv"))
      write.csv(calibration, ("ecoli_calibration.csv"))
      # Make them smaller; they're ~100MB each
      try({
        system("gzip -f results_with_evaluation.csv")
        system("gzip -f ecoli_calibration.csv")
      })
      # Plot and save calibration
      ggplot(calibration) +
        geom_point(aes(x = nominal_fdr, y = empirical_fdr)) +
        geom_errorbar(aes(x = nominal_fdr,
                          ymin = pmax(0, empirical_fdr - moe_95),
                          ymax = pmin(1, empirical_fdr + moe_95)
        )) +
        ggtitle("Calibration versus gold standard") +
        scale_y_continuous(limits = 0:1)
      ggsave(("calibration.pdf"), width = 4, height = 4)
      # Plot and save PR curve
      ggplot(calibration) +
        geom_point(aes(x = num_discoveries/nrow(gold_standard_symmetric),
                       y = 1-empirical_fdr)) +
        ggtitle("Calibration versus gold standard") +
        xlab("Recall") +
        ylab("Precision") +
        scale_y_continuous(limits = 0:1)
      ggsave(("pr.pdf"), width = 4, height = 4)

      # How many discoveries did we make per gene compared to the gold standard?
      merge(
        DF_plus_gs %>%
          dplyr::group_by(Gene2_name) %>%
          dplyr::summarize(n_targets = length(Gene1_name)) %>%
          dplyr::ungroup(),
        DF_plus_gs %>%
          subset(q<0.1) %>%
          dplyr::group_by(Gene2_name) %>%
          dplyr::summarize(n_targets = length(Gene1_name)) %>%
          dplyr::ungroup(),
        by = "Gene2_name"
      ) %>%
        ggplot(aes(x = n_targets.x, y = n_targets.y )) +
        geom_point() +
        geom_smooth(method = "lm", formula = y ~ x) +
        ggtitle("Target in-degree") +
        xlab("Gold standard") +
        ylab("Knockoff q<0.1")
      ggsave(("Target promiscuity.pdf"))

      merge(
        DF_plus_gs %>%
          dplyr::group_by(Gene1_name) %>%
          dplyr::summarize(n_targets = length(Gene2_name)) %>%
          dplyr::ungroup(),
        DF_plus_gs %>%
          subset(q<0.1) %>%
          dplyr::group_by(Gene1_name) %>%
          dplyr::summarize(n_targets = length(Gene2_name)) %>%
          dplyr::ungroup(),
        by = "Gene1_name"
      ) %>% ggplot(aes(x = n_targets.x, y = n_targets.y )) +
        geom_point() +
        geom_smooth(method = "lm", formula = y ~ x) +
        ggtitle("TF out-degree") +
        xlab("Gold standard") +
        ylab("Knockoff q<0.1")
      ggsave(("TF promiscuity.pdf"))
    })
  }
}

# How powerful is the ChIP-seq?
add_totals = function(X, name = "Mean", FUN = colMeans) {
  totals = t(matrix(FUN(X)))
  totals %>% 
    set_colnames(colnames(X)) %>%
    set_rownames(name) %>%
    rbind(X)
}

