suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("magrittr")
  library("glasso")
  library("mclust")
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

GAUSSIAN_KNOCKOFF_STRATEGIES = c("glasso", "sample", "shrinkage", "bagging")

# Set to T to run a fast test of the code. Results will not be biologically 
# meaningful; this is for software testing. 
test_mode = F

# We'll run several experiments with different handling of perturbations,
# confounding influences, and knockoff construction.
{
  conditions = rbind(
    # study random variation in knockoff construction
    expand.grid(
      knockoff_type = c("sample"),
      address_genetic_perturbations = c( F ),
      condition_on = c( "none" ),
      shrinkage_param = NA,      
      proportion_removed = c(0),
      do_simulate_y = F,
      seed = 1:5
    ),
    # study knockoff generation methods
    expand.grid(
      knockoff_type = c("sample", "naive", "shrinkage", "mixture"),
      address_genetic_perturbations = c( F ),
      condition_on = c( "none"),
      shrinkage_param = NA,
      proportion_removed = c(0, 0.1),
      do_simulate_y = T,
      seed = 1
    ),
    # study other knockoff generation methods
    expand.grid(
      knockoff_type = c( "glasso"),
      address_genetic_perturbations = c( F ),
      condition_on = c( "none"),
      shrinkage_param = c(0.01, 0.001, 0.0001),
      proportion_removed = c(0, 0.1),
      do_simulate_y = T,
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
      proportion_removed = c(0),
      do_simulate_y = F,
      seed = 1
    )
  )
  # Save time by not re-doing settings shared across experiments.
  # change the seed explicitly if you want replicates.
  conditions = conditions[!duplicated(conditions),] 
  # Currently auto-remove outliers only for Gaussians (not mixtures)
  conditions = subset(conditions, proportion_removed==0 | (knockoff_type %in% GAUSSIAN_KNOCKOFF_STRATEGIES ) )
  write.csv(conditions, "experiments_to_run.csv")
}

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
ecoli_networks = list()
{
  withr::with_dir(
    file.path(DATALAKE, "dream5/DREAM5_network_inference_challenge/Network3"),
    {
      ecoli_expression      = read.table("input data/net3_expression_data.tsv", header = T)
      ecoli_metadata        = read.table("input data/net3_chip_features.tsv", header = T, comment.char = "")
      ecoli_tf              = read.table("input data/net3_transcription_factors.tsv")
      ecoli_networks$dream5 = read.table("gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv")
      ecoli_anonymization   = read.table("anonymization/net3_gene_ids.tsv")
    }
  )
  # Make a couple useful alterations to metadata
  ecoli_metadata$Time_Rank = ave(ecoli_metadata$Time, ecoli_metadata$X.Experiment,
                                 FUN = function(x) rank(x, na.last = "keep")/max(rank(x)) )
  ecoli_metadata$Perturbations[is.na(ecoli_metadata$Perturbations)] = "none"

  # De-anonymize genes
  ecoli_networks$dream5 %<>% set_colnames(c("Gene1", "Gene2", "is_confirmed"))
  ecoli_anonymization_by_name = setNames(toupper(ecoli_anonymization$V2), ecoli_anonymization$V1)
  ecoli_networks$dream5 %<>% dplyr::mutate(Gene1_name = ecoli_anonymization_by_name[Gene1] )
  ecoli_networks$dream5 %<>% dplyr::mutate(Gene2_name = ecoli_anonymization_by_name[Gene2] )
  ecoli_networks$dream5$Gene1 = NULL
  ecoli_networks$dream5$Gene2 = NULL
  
  # Standardize representation
  ecoli_networks$dream5[["is_confirmed"]] %<>% equals(1)
  ecoli_networks$dream5_new_test$Gene1_name %<>% toupper
  ecoli_networks$dream5_new_test$Gene2_name %<>% toupper
  
  # Restrict to positive findings only
  ecoli_networks$dream5 %<>% subset(is_confirmed)
}

cat("\nPrepping gold standards based on new experiments in DREAM5.\n")
{
  ecoli_networks$dream5_new_test = read.csv("../DREAM5_ecoli_53_novel_hypotheses.csv", header = T, comment.char = "#")
  ecoli_networks$dream5_new_test$is_confirmed %<>% equals("yes")
  ecoli_networks$dream5_new_test$Gene1_name %<>% toupper
  ecoli_networks$dream5_new_test$Gene2_name %<>% toupper
}

cat("\nPrepping M3Dknockout-based gold standard.\n")
{
  # Find single-knockout experiments with no aliasing
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
  # Test for differential expression with limma
  ecoli_networks$M3Dknockout = list()
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
    ecoli_networks$M3Dknockout[[as.character(e)]] = 
      efit.contrast$p.value %>%
      as.data.frame() %>%
      dplyr::add_rownames(var = "Gene2_name") %>%
      tidyr::pivot_longer(cols = !Gene2_name, 
                          names_prefix = "DeletedGenesNoNA",
                          names_to = "Gene1_name", 
                          values_to = "p_value_KO") %>%
      mutate(X.Experiment = e)
  }
  ecoli_networks$M3Dknockout %<>% data.table::rbindlist()
  ecoli_networks$M3Dknockout = ecoli_networks$M3Dknockout[,c("Gene1_name", "Gene2_name", "p_value_KO", "X.Experiment")]
  ecoli_networks$M3Dknockout %<>% subset(!grepl(",", Gene1_name)) #Ignore double knockouts for now
  ecoli_networks$M3Dknockout$Gene1_name = ecoli_anonymization_by_name[ecoli_networks$M3Dknockout$Gene1_name] %>% toupper
  ecoli_networks$M3Dknockout$Gene2_name = ecoli_anonymization_by_name[ecoli_networks$M3Dknockout$Gene2_name] %>% toupper
  ecoli_networks$M3Dknockout %<>% 
    group_by(Gene1_name, Gene2_name) %>% 
    summarize(p_value_KO = poolr::fisher(p_value_KO)$p)
  ecoli_networks$M3Dknockout$q_value_KO = p.adjust(ecoli_networks$M3Dknockout$p_value_KO, method = "fdr")
  ecoli_networks$M3Dknockout %<>% mutate(targetIsDecoy = grepl("DECOY", Gene2_name))
  ecoli_networks$M3Dknockout %>%
    ggplot()  +
    geom_histogram(aes(x = p_value_KO), bins = 100)  + 
    facet_wrap(~ifelse(targetIsDecoy, "Decoy", "Gene"), ncol = 1, scales = "free_y") +
    xlab("P value") +
    ggtitle( "Differential expression after knockout")
  ggsave("knockout_differential_expression_inflation.svg", width = 5, height = 5)
  ggsave("knockout_differential_expression_inflation.pdf", width = 5, height = 5)
  ecoli_networks$M3Dknockout$is_confirmed = NA
  ecoli_networks$M3Dknockout$is_confirmed[ecoli_networks$M3Dknockout$q_value_KO < 0.05] = T
  ecoli_networks$M3Dknockout %<>% subset(!is.na(is_confirmed))
  ecoli_networks$M3Dknockout %<>% subset(!targetIsDecoy)
}

cat("\nPrepping regulonDB 10.9 gold standard.\n")
# Load curated interactions from regulonDB v10.9
{
  tidy_regulondb = function(regulondb10_9) {
    regulondb10_9 %<>%
      subset(V8=="gene" & V7=="tf") %>%
      extract(c(2, 4)) %>%
      set_colnames(c("Gene1_name", "Gene2_name"))
    regulondb10_9[["Gene1_name"]] %<>% toupper
    regulondb10_9[["Gene2_name"]] %<>% toupper
    regulondb10_9
  }
  ecoli_networks$regulondb10_9 =
    read.table(
      sep = "\t",
      file.path(DATALAKE, "modern_ecoli/regulondb_10.9/full/genetic_network.txt")
    ) 
  ecoli_networks$regulondb10_9_additional_binding_support =
    ecoli_networks$regulondb10_9 %>%
    subset(grepl("sequence|binding|SELEX|mutation", ignore.case = T, V6)) %>% # Evidence of direct binding
    tidy_regulondb
  ecoli_networks$regulondb10_9 %<>% tidy_regulondb
  ecoli_networks$regulondb10_9$is_confirmed = T
  ecoli_networks$regulondb10_9_additional_binding_support$is_confirmed = T
}

cat("\nPrepping RegulonDB ChIP gold standard.\n")
# Load ChIP data from regulonDB v10.9
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
  ecoli_networks$chip =
    list.files(
      file.path(DATALAKE,"modern_ecoli/regulondb_10.9/highthroughput/"),
      full.names = T
    ) %>%
    lapply(read.table, sep = "\t") %>%
    data.table::rbindlist() %>%
    set_colnames(schema) %>%
    as.data.frame()
  # MelR regulates MelAB. This is missing from RegulonDB but there is high quality data from an included study.
  # See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1308901/
  ecoli_networks$chip %<>% rbind(read.csv(row.names = 1, check.names = F, text = 
                                            "RegulonDB Dataset Identifier,Dataset Type,Dataset Name,Transcription Factor Name,Effect,Regulated Object,TFBs Absolute genome left position,TFBs Absolute genome right position,DNA sequence,Growth Condition,Reference PMID
    1,<NA>,TF CHIP,<NA>,MelR,<NA>,MelA,<NA>,<NA>,<NA>,<NA>,16301522
    2,<NA>,TF CHIP,<NA>,MelR,<NA>,MelB,<NA>,<NA>,<NA>,<NA>,16301522"
  ))
  # Filter out some samples
  ecoli_chip_metadata = read.csv(file.path(DATALAKE,"modern_ecoli/regulondb_10.9/table_chip_controls.csv"))
  ecoli_networks$chip %<>% merge(ecoli_chip_metadata, by.x = "Reference PMID", by.y = "PMID")
  ecoli_networks$chip %<>% subset( 
    "keep" == Recommendation |  # keep whole studies
      mapply(grepl, `Transcription Factor Name`, Recommendation, ignore.case = T) # keep only certain regulators from a couple studies
  )
  # Tidy to match other gold standards
  ecoli_networks$chip %<>%
    extract(c("Transcription Factor Name", "Regulated Object")) %>%
    set_colnames(c("Gene1_name", "Gene2_name")) %>%
    dplyr::mutate(is_confirmed = T)
  ecoli_networks$chip = ecoli_networks$chip[!duplicated(ecoli_networks$chip),]
  
  # Fix some genes with multiple names
  ecoli_networks$chip$Gene1_name[ecoli_networks$chip$Gene1_name=="Sigma70"] = "rpoD"
  ecoli_networks$chip$Gene1_name[ecoli_networks$chip$Gene1_name=="Sigma38"] = "rpoS"
  ecoli_networks$chip$Gene1_name[ecoli_networks$chip$Gene1_name=="Sigma32"] = "rpoH"
  ecoli_networks$chip$Gene1_name[ecoli_networks$chip$Gene1_name=="Cra"    ] = "fruR"
  ecoli_networks$chip$Gene1_name[ecoli_networks$chip$Gene1_name=="H-NS"   ] = "hns"
  ecoli_networks$chip$Gene1_name[ecoli_networks$chip$Gene1_name=="IHF"    ] = "ihfA"
  # Fix some annotations like "WITHIN bcsE" to "bcsE"
  ecoli_networks$chip$Gene2_name %<>% gsub("WITHIN |END OF |PSEUDOGENE|-[0-9]", "", ignore.case = T, .)
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
  ecoli_networks$chip %<>% dplyr::mutate(Gene2_name = sapply(Gene2_name, separate_ecoli_operon))
  ecoli_networks$chip %<>% tidyr::separate_rows(Gene2_name) 
  # Remove errant suffixes like -1 and -2
  ecoli_networks$chip %<>% subset(nchar(Gene2_name)>=3)
  ecoli_networks$chip =
    ecoli_networks$chip %>%
    subset(Gene1_name=="ihfA") %>%
    mutate(Gene1_name="ihfB") %>%
    rbind(ecoli_networks$chip)
  ecoli_networks$chip[["Gene1_name"]] %<>% toupper
  ecoli_networks$chip[["Gene2_name"]] %<>% toupper
}

cat("\nPrepping RegulonDB knockout gold standard.\n")
# Load ChIP data from regulonDB v10.9
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
  ecoli_networks$regulonDB_knockout =
    list.files(
      file.path(DATALAKE,"modern_ecoli/regulondb_10.9/highthroughput2/"),
      full.names = T
    ) %>%
    lapply(read.table, sep = "\t") %>%
    data.table::rbindlist() %>%
    set_colnames(schema) %>%
    as.data.frame()
  # Tidy to match other gold standards
  ecoli_networks$regulonDB_knockout %<>%
    extract(c("Transcription Factor Name", "Regulated Object")) %>%
    set_colnames(c("Gene1_name", "Gene2_name")) %>%
    dplyr::mutate(is_confirmed = T)
  ecoli_networks$regulonDB_knockout = ecoli_networks$regulonDB_knockout[!duplicated(ecoli_networks$regulonDB_knockout),]
  ecoli_networks$regulonDB_knockout %<>% subset(Gene1_name != "nd")
  ecoli_networks$regulonDB_knockout %<>% subset(Gene2_name != "nd")
  # Fix some genes with multiple names
  ecoli_networks$regulonDB_knockout$Gene1_name[ecoli_networks$regulonDB_knockout$Gene1_name=="Sigma70"] = "rpoD"
  ecoli_networks$regulonDB_knockout$Gene1_name[ecoli_networks$regulonDB_knockout$Gene1_name=="Sigma38"] = "rpoS"
  ecoli_networks$regulonDB_knockout$Gene1_name[ecoli_networks$regulonDB_knockout$Gene1_name=="Sigma32"] = "rpoH"
  ecoli_networks$regulonDB_knockout$Gene1_name[ecoli_networks$regulonDB_knockout$Gene1_name=="Cra"    ] = "fruR"
  ecoli_networks$regulonDB_knockout$Gene1_name[ecoli_networks$regulonDB_knockout$Gene1_name=="H-NS"   ] = "hns"
  ecoli_networks$regulonDB_knockout$Gene1_name[ecoli_networks$regulonDB_knockout$Gene1_name=="IHF"    ] = "ihfA"
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
  ecoli_networks$regulonDB_knockout %<>% dplyr::mutate(Gene2_name = sapply(Gene2_name, separate_ecoli_operon))
  ecoli_networks$regulonDB_knockout %<>% tidyr::separate_rows(Gene2_name) 
  ecoli_networks$regulonDB_knockout[["Gene1_name"]] %<>% toupper
  ecoli_networks$regulonDB_knockout[["Gene2_name"]] %<>% toupper
}

# Augment the ChIP benchmark with all genes of any targeted transcription unit.
augment_gold_standard = function(gold_standard){
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
  gold_standard_augmented = merge(gold_standard, ecoli_tu, by = "Gene2_name", all.x = T)
  gold_standard_augmented = merge(gold_standard_augmented, ecoli_tu, by = "Transcription-Units", all.y = T)
  gold_standard_augmented %<>% dplyr::mutate(not_in_original = Gene2_name.x != Gene2_name.y)
  gold_standard_augmented[["Gene2_name"]] = gold_standard_augmented[["Gene2_name.y"]]
  gold_standard_augmented[["Gene2_name.x"]] = NULL
  gold_standard_augmented[["Gene2_name.y"]] = NULL
  doops = duplicated(gold_standard_augmented[c("Gene1_name","Gene2_name")])
  gold_standard_augmented = gold_standard_augmented[!doops, ]
  gold_standard_augmented %<>% subset(!is.na(Gene1_name))
  
  # Notify user of augmentation effects
  cat("Augmenting a gold standard with additional targets in the same operon as a known target.\n")
  cat("Dimensions and out-degrees before and after augmentation:\n")
  dim(gold_standard) %>% cat("\n")
  dim(gold_standard_augmented) %>% cat("\n")
  table(gold_standard$Gene1_name) %>% sort %>% print
  table(gold_standard_augmented$Gene1_name) %>% sort %>% print

  gold_standard_augmented
}
ecoli_networks$chip_tu_augmented = augment_gold_standard(ecoli_networks$chip)
ecoli_networks$M3Dknockout_tu_augmented = augment_gold_standard(gold_standard = ecoli_networks$M3Dknockout)
ecoli_networks$regulonDB_knockout_tu_augmented = augment_gold_standard(gold_standard = ecoli_networks$regulonDB_knockout)

# Take a network and fill in missing entries as negatives (A does not regulate B).
# For e.g. ChIP, this is well motivated as long as it's sensitive enough.
# For curated collections with positives only, it's not as well motivated but we have no other choice. 
# For collections already having explicit negatives, default behavior is to return it unaltered. 
add_negatives_if_none = function(DF, possible_targets = unique(DF$Gene2_name), warn = T, force = F){
  if(!force & any(!DF$is_confirmed)){
    if(warn){ warning("Explicit negatives already present. Skipping.\n") }
    return(DF)
  }
  negatives = expand.grid(
    Gene1_name = unique(DF[["Gene1_name"]]), 
    Gene2_name = unique(possible_targets)
  )
  DF %<>% merge(negatives, all.y = T)
  DF[["is_confirmed"]][is.na(DF[["is_confirmed"]])] = F
  DF
}
# Test
# add_negatives(data.frame(Gene1_name = 1:2, Gene2_name = 1:2, is_confirmed = T))
ecoli_networks %<>% lapply(add_negatives_if_none, warn = F, force = F)



{
  cat("\nPrepping gold standards based on intersection of knockout and ChIP data.\n")
  # label a pair as:
  # positive if in both
  # negative if in neither
  # unknown otherwise
  fill_in_confident_results = function(is_confirmed_chip, is_confirmed_knockout){
    do_one = function(n){
      if(n==2){
        return(T)
      } else if(n==1){
        return(NA)
      } else {
        return(F)
      }
    }
    sapply(is_confirmed_chip + is_confirmed_knockout, do_one)
  }
  ecoli_networks$chip_intersection_M3Dknockout = merge(
    ecoli_networks$chip_tu_augmented, 
    ecoli_networks$M3Dknockout_tu_augmented,
    by = c("Gene1_name", "Gene2_name"),
    suffixes = c("_chip", "_knockout")
  ) %>%
    mutate(is_confirmed = fill_in_confident_results(is_confirmed_chip, is_confirmed_knockout))  %>% 
    subset(!is.na(is_confirmed)) 
  
  cat("Gold standard size:\n")
  ecoli_networks$M3Dknockout_tu_augmented %>% 
    extract(c("Gene1_name", "is_confirmed")) %>%
    table %>% 
    print
  
  ecoli_networks$chip_intersection_RegulonDB_knockout = merge(
    ecoli_networks$chip_tu_augmented, 
    ecoli_networks$regulonDB_knockout_tu_augmented,
    by = c("Gene1_name", "Gene2_name"),
    suffixes = c("_chip", "_knockout")
  ) %>%
    mutate(is_confirmed = fill_in_confident_results(is_confirmed_chip, is_confirmed_knockout))  %>% 
    subset(!is.na(is_confirmed)) 
  
  cat("Gold standard size:\n")
  ecoli_networks$chip_intersection_RegulonDB_knockout %>% 
    extract(c("Gene1_name", "is_confirmed")) %>%
    table %>% 
    print
}  

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

#' Helper function for below. Given a table with q-values and T/F/NA ground truth, compute calibration at various FDR's. 
#' 
compute_calibration_with_unknowns = function(DF_plus_gs){
  calibration = data.frame(nominal_fdr = (c(1:100)/100) %>% c(quantile(DF_plus_gs$q, probs = c(1:100)/100, na.rm = T)),
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
  return(calibration)
}

#' This handles a difference between a causal DAG and the corresponding MRF structure (the moral graph).
#' Edges in the MRF can occur between parents & children (good) or spouses (bad), For us, all TF-TF pairs 
#' are potential spouses, so we have the option to omit them from benchmarks. 
#' 
omit_possible_spouses = function(DF_plus_gs){
  subset(DF_plus_gs, !(Gene1 %in% ecoli_tf) | !(Gene2 %in% ecoli_tf)  )
}

#' Check calibration against all gold standards, accounting for unknown edge orientation, possible spouses, and incompleteness of gold standards.
#' This will be used later to score the results.
#' 
check_against_gold_standards = function(DF){
  for( gold_standard_name in names(ecoli_networks) ){
    # Check gold standard and set things up
    cat("Evaluating against ", gold_standard_name, " gold standard\n")
    gold_standard = ecoli_networks[[gold_standard_name]]
    stopifnot(all(c("Gene1_name", "Gene2_name", "is_confirmed") %in% colnames(gold_standard)))
    stopifnot("Missing values in gold standards are not allowed."=!any(is.na(gold_standard$is_confirmed)))
    gold_standard = gold_standard[,c("Gene1_name", "Gene2_name", "is_confirmed")]
    dir.create(gold_standard_name, recursive = T, showWarnings = F)
    
    # Include forwards & backwards edges -- we don't expect to get the directionality right here.
    gold_standard_reversed = gold_standard
    gold_standard_reversed[1:2] = gold_standard[2:1]
    gold_standard_symmetric = rbind(gold_standard, gold_standard_reversed)
    # If forwards is true and backwards is false, keep only the trues in both directions.
    gold_standard_symmetric %<>% dplyr::arrange(desc(is_confirmed))
    gold_standard_symmetric = gold_standard_symmetric[!duplicated(gold_standard_symmetric$Gene1_name,
                                                                  gold_standard_symmetric$Gene2_name),]
    
    calibration_across_spouse_handling = vector("list", 2)
    for(do_omit_possible_spouses in c(F, T)){
      DF_plus_gs = merge(DF, gold_standard_symmetric, by = c("Gene1_name", "Gene2_name"), all.x = F, all.y = F)
      if(do_omit_possible_spouses){
        DF_plus_gs = omit_possible_spouses(DF_plus_gs)
      }
      # Computing q-values and then subsetting is not the same as subsetting and then computing q-values.
      # We want them w.r.t. the subset we are able to evaluate. 
      DF_plus_gs$q = rlookc::knockoffQvals(DF_plus_gs$knockoff_stat)
      calibration_across_spouse_handling[[do_omit_possible_spouses+1]] = compute_calibration_with_unknowns(DF_plus_gs)
      calibration_across_spouse_handling[[do_omit_possible_spouses+1]][["do_omit_possible_spouses"]] = do_omit_possible_spouses
    }
    calibration = Reduce(rbind, calibration_across_spouse_handling)


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

#' Add a row to a table, showing the means of the input columns.
#' 
#' 
add_totals = function(X, name = "Mean", FUN = colMeans) {
  totals = t(matrix(FUN(X)))
  totals %>% 
    set_colnames(colnames(X)) %>%
    set_rownames(name) %>%
    rbind(X)
}


#' Take a network and make it "unbalanced" by removing
#' different fractions of positive or negative examples.
#'
makeUnbalancedNetwork = function(network, prop_positives, prop_negatives, seed = 0){
  network %<>% subset(!is.na(is_confirmed))
  if(all(network$is_confirmed, na.rm = T)){
    stop("Cannot make it unbalanced with no negative examples.\n")
  }
  set.seed(seed)
  include_if_positive = rbinom(prob = prop_positives, size = 1, n = nrow(network))
  include_if_negative = rbinom(prob = prop_negatives, size = 1, n = nrow(network))
  network$include = ifelse( network$is_confirmed, include_if_positive, include_if_negative) %>% as.logical
  network %<>% subset(include)
  network
}
# Quick test
# makeUnbalancedNetwork(data.frame(is_confirmed = rep(0:1, 20)), 0.5, 0.5)

# Deploy on networks with large but incomplete positive and negative sets
for(bias in c("positive", "negative")){
  for(network in c("chip_intersection_M3Dknockout", "chip_intersection_RegulonDB_knockout", "regulondb10_9")){
    prop_positives = ifelse(bias == "positive", 1, 0.5)
    prop_negatives = ifelse(bias == "negative", 1, 0.5)
    ecoli_networks[[paste0(network, "_biased_", bias)]] = 
      makeUnbalancedNetwork( ecoli_networks[[network]], 
                             prop_positives, prop_negatives) 
  }
}

# Compare basic descriptive properties of the gold standards used in this study. 
{
  showNetworkProperties = function(network, name){
    return(
      Reduce(
        f = rbind,
        list(
          data.frame(
            count = network %>%
              dplyr::group_by(Gene2_name) %>% 
              dplyr::summarise(degree = sum(is_confirmed)) %>%
              extract2("degree"),
            summary = "outDegree", 
            network = name
          ),
          data.frame(
            count = network %>%
              dplyr::group_by(Gene1_name) %>% 
              dplyr::summarise(degree = sum(is_confirmed)) %>%
              extract2("degree"),
            summary = "inDegree", 
            network = name
          ),
          totals = 
            network$is_confirmed %>% 
            table %>%
            setNames(c("Negatives", "Positives"))  %>%
            as.data.frame %>%
            setNames(c("summary", "count")) %>% 
            dplyr::mutate(network = name)
        )
      )
    )
  }
  
  networkProperties = mapply(showNetworkProperties, ecoli_networks, names(ecoli_networks), SIMPLIFY = F) %>%
    Reduce(f = rbind) %>%
    subset(count != 0)
  networkProperties %>% 
    dplyr::mutate(summary = factor(summary, levels = c("inDegree", "outDegree", "Positives", "Negatives"))) %>%
    dplyr::group_by(network, summary) %>% 
    dplyr::summarise_all(median) %>% 
    tidyr::pivot_wider(names_from = summary, values_from = count) %>%
    write.csv("GoldStandardProperties.csv")
}

#' Test gaussian mixture models on the M3 data.
#' 
compare_mixture_models = function(X, output_name = "mclust_bic") {
  mclust_bic = list(
    log2_data     = mclust::mclustBIC(X %>% log2, G = c(1:10, 15, 20, 25, 50, 100)),
    original_data = mclust::mclustBIC(X,          G = c(1:10, 15, 20, 25, 50, 100))
  )
  saveRDS(mclust_bic, paste0(output_name, ".Rdata"))
  write.csv(
    mclust_bic %>%
      lapply(divide_by, 1e3) %>%
      lapply(round) %>% 
      mapply(cbind, ., transformation = names(mclust_bic), SIMPLIFY = F) %>% 
      Reduce(f=rbind),
    paste0(output_name, ".csv")
  )
}

#' Estimate covariance many different ways. 
#'
#' Usually called via make_knockoffs, not directly. For 'naive' and 'mixtures' methods, returns NULL.
#' For "sample", "bagging", "shrinkage", "glasso", returns a P by P covariance matrix estimate; P=ncol(X).
estimate_covariance = function(X, method){
  if (method == "sample"){
    Sigma = cov(X)
  } else if (method == "bagging"){
    Sigma = boot::boot(X, statistic = function(original, idx) c(cov(original[idx,])), R = 50)
    Sigma = matrix(colSums(Sigma$t), ncol = ncol(X))
  } else if (method == "shrinkage"){
    if(is.na(shrinkage_param)){
      off_diagonal_shrinkage = corpcor::estimate.lambda(X)
    } else {
      off_diagonal_shrinkage = shrinkage_param
    }
    Sigma = as(corpcor::cor.shrink(X), "matrix")
  } else if (method == "glasso"){
    glasso_out = glasso::glasso(cov(X), rho = shrinkage_param, penalize.diagonal = F)
    Sigma = glasso_out$w
    write.table( mean(glasso_out$wi!=0), ("glasso_wi_fraction_nonzero.csv"))
  } else if (method %in% c("mixture", "naive")) {
    return(NULL)
  } else {
    stop(paste0("Invalid 'method' arg ", method, " in estimate_covariance."))
  }
  return(Sigma)
}

#' Set knockoffs equal to original data for the most difficult observations. Then re-fit. 
#'
#' @param X data
#' @param proportion_removed fraction to treat as outliers.
#' @param ... passed to estimate_covariance 
#' 
#' Usually called via make_knockoffs, not directly.
#' 
make_knockoffs_remove_outliers = function(X, proportion_removed, method, Sigma = NULL){
  if(is.null(Sigma)){
    Sigma = estimate_covariance(X, method = method)
  }
  if(proportion_removed==0){
    Xin = X
    Sigma_in = Sigma
  } else {
    # remove outliers, re-estimate covariance
    Z = scale(X, center = T, scale = F)
    mahalanobis = diag(Z %*% solve(Sigma, t(Z)))
    outliers = mahalanobis>quantile(mahalanobis, 1-proportion_removed)
    Xout = X[outliers,]
    Xin = X[!outliers,]
    Sigma_in = estimate_covariance(Xin, method)
  }

  # For TF targets, use fast leave-one-out knockoffs (LOOKs).
  # This output contains all of them in a compact form.
  # Realize them later on using rlookc::formOneLook.
  looks_compact_inliers = rlookc::generateLooks(
    Xin,
    mu = colMeans(Xin),
    Sigma = Sigma_in,
    statistic = fast_lasso_penalty,
    output_type = "knockoffs_compact",
    vars_to_omit = seq(ncol(X))
  )

  # That was just the inliers. For the outliers, we'll copy the original data. 
  # Now we need to intersperse them properly.
  # Empty containers in the right shape
  if(proportion_removed==0){
    return(looks_compact_inliers)
  } else {
    looks_compact = list(
      knockoffs = matrix(0, nrow = nrow(X), ncol = ncol(X)),
      updates =
        rep(
          x = list( list(
            mean_update_left  = matrix(0, nrow = nrow(X)),
            mean_update_right = matrix(0, ncol = ncol(X) - 1),
            sqrt_cov_update   = matrix(0, ncol = ncol(X) - 1),
            random_update     = matrix(0, nrow = nrow(X))
          ) ),
          length.out = ncol(X)
        ),
      vars_to_omit = 1:ncol(X)
    )
    # Fill it in 
    looks_compact[["knockoffs"]][!outliers,] = looks_compact_inliers$knockoffs
    looks_compact[["knockoffs"]][outliers,] = X[outliers,]
    for(left_out in names(looks_compact[["updates"]])){
      for(update_type in c( "mean_update_left",
                            "random_update")  ){
        # Inliers are taken care of by previous code.
        looks_compact          [["updates"]][left_out][update_type][!outliers,] = 
          looks_compact_inliers[["updates"]][left_out][update_type] 
        # Outliers need no update.
        # Copying part of the data is copying part of the data whether or not a given var is included. 
        looks_compact          [["updates"]][left_out][update_type][outliers,] = 0
      }     
      for(update_type in c( "mean_update_right", 
                            "sqrt_cov_update" )  ){
        looks_compact          [["updates"]][left_out][update_type] = 
          looks_compact_inliers[["updates"]][left_out][update_type] 
      }
    }
  }
  return(looks_compact)
}

#' This is the main user-facing code to make knockoffs during each experiment.
#'
#' Takes care of correcting for PC's and making knockoffs under all different types of models and covariance estimates. 
#' 
make_knockoffs = function(ecoli_tf_expression_working_copy, confounder_indicators, knockoff_type, proportion_removed = 0, Sigma = NULL){
  X = cbind(ecoli_tf_expression_working_copy, confounder_indicators)
  if(knockoff_type == "mixture"){
    ecoli_tf_expression_working_copy_knockoffs =
      rlookc::computeGaussianMixtureKnockoffs(
        X,
        do_high_dimensional = T,
        hard_assignments = clusters$cluster
      )
  } else if(knockoff_type == "naive"){
    ecoli_tf_expression_working_copy_knockoffs = apply(X, 2, sample, replace = F)
  } else if(knockoff_type %in% GAUSSIAN_KNOCKOFF_STRATEGIES){
    # These methods differ only in the covariance estimate Sigma, done above, but
    # otherwise they're the same.
    looks_compact = make_knockoffs_remove_outliers(X, proportion_removed = proportion_removed, method = knockoff_type, Sigma = Sigma)
    ecoli_tf_expression_working_copy_knockoffs = looks_compact$knockoffs
  } else {
    stop("Knockoff type not recognized.\n")
  }
  # Sometimes, dedicated efficient code to make LOOKs is not yet available.
  # So that the code downstream can remain the same, I pretend to make
  # LOOKs inside the same data structure used above, even when they just contain
  # the full knockoffs unmodified.
  if(knockoff_type %in% c("naive", "mixture")){
    looks_compact = list(
      knockoffs = ecoli_tf_expression_working_copy_knockoffs,
      updates =
        rep(
          x = list( list(
            mean_update_left  = matrix(0, nrow = 805),
            mean_update_right = matrix(0, ncol = ncol(X) - 1),
            sqrt_cov_update   = matrix(0, ncol = ncol(X) - 1),
            random_update     = matrix(0, nrow = 805)
          ) ),
          length.out = 334
        ),
      vars_to_omit = 1:334
    )
  }
  
  # If confounders were included in knockoff generation, strip them out.
  ecoli_tf_expression_working_copy_knockoffs = ecoli_tf_expression_working_copy_knockoffs[,1:334]
  list(ecoli_tf_expression_working_copy_knockoffs, looks_compact)
}

# One experiment, e.g. "gaussian knockoffs with glasso 0.001 covariance estimate, correcting for 20 principal components 
# and with no special handling of genetic perturbations."
do_one = function(condition_index, omit_knockout_evaluation_samples = F, reuse_results = F){
  stopifnot(condition_index %in% seq_along(conditions[[1]]))
  withr::with_message_sink(file.path("logs", condition_index), {
    withr::with_output_sink(file.path("logs", condition_index), {
      cat("\n======================================\n")
      cat("Experiment", condition_index, "of", nrow(conditions), '\n')
      cat(  "======================================\n")
      attach(conditions[condition_index,], warn.conflicts = F)
      tryCatch(set.seed(seed), error = function(e){warning("Random seed not found and not set.\n")})
      working_directory = ""
      conditions$omit_knockout_evaluation_samples = omit_knockout_evaluation_samples
      for(j in seq_along(conditions[condition_index,])){
        prefix = paste0( colnames(conditions)[[j]], "=", conditions[condition_index, j] )
        working_directory %<>% file.path(prefix, .)
      }
      
      ecoli_metadata_working_copy = ecoli_metadata
      ecoli_expression_working_copy = ecoli_expression
      ecoli_tf_expression_working_copy = ecoli_tf_expression
      if(omit_knockout_evaluation_samples){
        samples_to_keep = !(ecoli_metadata$X.Experiment %in% ecoli_networks$knockout$X.Experiment)
        ecoli_metadata_working_copy      %<>% extract(samples_to_keep, )
        ecoli_expression_working_copy    %<>% extract(samples_to_keep, )
        ecoli_tf_expression_working_copy %<>% extract(samples_to_keep, )
        stopifnot(nrow(ecoli_metadata_working_copy)==nrow(ecoli_expression_working_copy))
      } 
      
      cat("\nTesting condition:", working_directory, "\n")
      dir.create(working_directory, recursive = T, showWarnings = F)
      withr::with_dir(
        working_directory,
        tryCatch( expr = {
          if( address_genetic_perturbations ){
            # When predicting targets of G, if G was knocked out, set expression to 0.
            # If G was overexpressed, leave expression as is.
            # Do this before constructing knockoffs.
            for(i in seq_along(ecoli_metadata_working_copy[["DeletedGenes"]])){
              g = ecoli_metadata_working_copy[["DeletedGenes"]][[i]]
              if(!is.na(g)){
                ecoli_expression_working_copy[i, strsplit(g, ",")[[1]]] = 0
              }
            }
            # When predicting regulators of G, omit or impute any expression value that has been intervened on.
            # To do this easily later, here's a little function.
            
            #' Determine samples to omit.
            #'
            #' Input is a gene label (matched to metadata, so anonymous as in "G2388", not 2388 or "RELA")
            #' Output is a set of samples where that gene is intervened upon.
            #' These are useless for determining regulators of the target gene.
            get_samples_to_include = function(target_name){
              stopifnot(grepl("^G", target_name)) # supposed to be anonymized names, not non-anonymized
              !or(
                ecoli_metadata_working_copy$DeletedGenes %>%
                  strsplit(",") %>%
                  sapply(function(x) is.element(target_name, x)),
                ecoli_metadata_working_copy$OverexpressedGenes %>%
                  strsplit(",") %>%
                  sapply(function(x) is.element(target_name, x))
              )
            }
          } else {
            get_samples_to_include = function(x) rep(T, nrow(ecoli_metadata_working_copy))
          }
          
          # Construct indicators for confounders
          {
            perturbations = model.matrix( ~Perturbations+0, ecoli_metadata_working_copy)
            perturbations %<>% extract(,-which(colnames(perturbations) == "Perturbationsnone"))
            pca_results = irlba::prcomp_irlba(ecoli_expression_working_copy, n = 50)
            if(condition_on == "none"){
              confounder_indicators = NULL
            } else if(condition_on == "pert_labels"){
              confounder_indicators = perturbations
            } else if(grepl("^pert_labels_plus_pca", condition_on)){
              n_pc = 10
              maybe_n_pc = gsub("^pert_labels_plus_pca", "", condition_on) %>% as.numeric
              if(!is.na(maybe_n_pc)){
                n_pc = maybe_n_pc
              }
              confounder_indicators = cbind( pca_results$x[,1:n_pc], perturbations)
            } 
          }
          
          # Make knockoffs. 
          cat("\nConstructing knockoffs...\n")
          l = make_knockoffs(ecoli_tf_expression_working_copy = ecoli_tf_expression_working_copy,
                             confounder_indicators = confounder_indicators,
                             knockoff_type = knockoff_type, 
                             proportion_removed = proportion_removed)
          print("... knockoffs ready.")
          ecoli_tf_expression_working_copy_knockoffs = l[[1]]
          looks_compact = l[[2]]
          rm(l)
          write.csv(ecoli_tf_expression_working_copy_knockoffs, "knockoffs.csv")
          
          # Check power: how different are knockoffs from original features?
          cat("\nChecking correlation of variables vs their knockoffs.")
          cor_with_knockoffs = cor(ecoli_tf_expression_working_copy_knockoffs, ecoli_tf_expression_working_copy) %>% diag
          data.frame( cor_with_knockoffs,
                      is_decoy = grepl("DECOY", ecoli_anonymization_by_name[colnames(ecoli_tf_expression_working_copy)])) %>%
            ggplot() +
            geom_histogram(aes(x = cor_with_knockoffs, fill = is_decoy, bins = 30)) +
            ggtitle("Correlation of each TF with its knockoff")
          ggsave(("correlations_with_knockoffs.pdf"), width = 5, height = 3)
          
          # Check calibration for batch adjustment.
          # No association should be detected with factors explicitly conditioned on.        
          cat("\nChecking association with known confounders.")
          if(condition_on!="none"){
            confounder_association_qvals =
              rlookc::knockoffQvals(
                cor(confounder_indicators, ecoli_tf_expression_working_copy) -
                  cor(confounder_indicators, ecoli_tf_expression_working_copy_knockoffs)
              )
            pdf("confounder_association_qvals.pdf", width = 5, height = 2.5)
            hist(confounder_association_qvals, 40, xlim = 0:1)
            dev.off()
          }
          
          # Check calibration with simulated Y
          if(do_simulate_y){
            cat("\nChecking calibration\n")
            Sigma = estimate_covariance(cbind(ecoli_tf_expression_working_copy, confounder_indicators), method = knockoff_type)
            calibration_typical = rlookc::calibrate__simulateY(
              cores = 1,
              rng_seed = seed,
              n_sim = 1000,
              X = ecoli_tf_expression_working_copy %>% sweep(2, colMeans(ecoli_tf_expression_working_copy), "-"),
              knockoffs = replicate(10, simplify = F, {
                make_knockoffs(ecoli_tf_expression_working_copy = ecoli_tf_expression_working_copy,
                               confounder_indicators = confounder_indicators,
                               knockoff_type = knockoff_type, 
                               Sigma = Sigma,
                               proportion_removed = proportion_removed)[[1]] %>%
                  sweep(2, colMeans(ecoli_tf_expression_working_copy), "-") 
              }),
              statistic = fast_lasso_penalty,
              plot_savepath = ("average_case_calibration.pdf")
            )
            write.csv(calibration_typical$calibration$fdr, ("average_case_calibration.csv"))
          }
          # Deploy on non-TF targets
          if(!reuse_results){
            cat("\nComputing LASSO paths. This will take a long time.\n")
            target_genes = colnames(ecoli_expression_working_copy) %>% setdiff(ecoli_tf[[1]])
            w = list()
            for(i in seq_along(target_genes)){
              target = target_genes[[i]]
              if(i %% 100 == 0){
                cat( "\n  Target", i, "reached. Remaking knockoffs. ")
                l = make_knockoffs(ecoli_tf_expression_working_copy = ecoli_tf_expression_working_copy,
                                   confounder_indicators = confounder_indicators,
                                   knockoff_type = knockoff_type, 
                                   proportion_removed = proportion_removed)
                ecoli_tf_expression_working_copy_knockoffs = l[[1]]
                looks_compact = l[[2]]
                rm(l)
              }
              samples_to_include = get_samples_to_include(target)
              w[[target]] = fast_lasso_penalty(
                y   = ecoli_expression_working_copy[samples_to_include,target],
                X   = ecoli_tf_expression_working_copy[samples_to_include,],
                X_k = ecoli_tf_expression_working_copy_knockoffs[samples_to_include,]
              )
            }
            length(w[[1]])
            length(w)
            # Deploy on TF targets
            cat("\n")
            w_tf = list()
            for(i in seq_along(ecoli_tf[[1]])){
              target = ecoli_tf[[1]][[i]]
              if(i %% 100 == 0){
                cat( "\n  Target", i, "reached. Remaking knockoffs. ")
                l = make_knockoffs(ecoli_tf_expression_working_copy = ecoli_tf_expression_working_copy,
                                   confounder_indicators = confounder_indicators,
                                   knockoff_type = knockoff_type)
                ecoli_tf_expression_working_copy_knockoffs = l[[1]]
                looks_compact = l[[2]]
                rm(l)
              }
              knockoffs_minus_i = rlookc::formOneLook(
                knockoffs = looks_compact$knockoffs,
                updates = looks_compact$updates,
                vars_to_omit = looks_compact$vars_to_omit,
                k = i
              )[,seq(ncol(ecoli_tf_expression_working_copy)-1)]
              samples_to_include = get_samples_to_include(target)
              w_tf[[target]] = fast_lasso_penalty(
                y   = ecoli_expression_working_copy[samples_to_include,target],
                X   = ecoli_tf_expression_working_copy[samples_to_include,-i],
                X_k = knockoffs_minus_i[samples_to_include,]
              )
            }
            length(w_tf[[1]])
            length(w_tf)
            # Assemble results nicely
            DF_tf_target = DF_tf_tf = list()
            for(k in seq_along(ecoli_tf[[1]])){
              DF_tf_target[[k]] = data.frame(
                Gene1 = ecoli_tf[[1]][-k],
                Gene2 = ecoli_tf[[1]][ k],
                knockoff_stat = w_tf[[k]]
              )
            }
            for(k in seq_along(target_genes)){
              DF_tf_tf[[k]] = data.frame(
                Gene1 = ecoli_tf[[1]],
                Gene2 = target_genes[ k],
                knockoff_stat = w[[k]]
              )
            }
            DF = data.table::rbindlist(c(DF_tf_tf, DF_tf_target))
            signedsqrt = function(x) sign(x)*sqrt(abs(x))
            hist(DF$knockoff_stat %>% signedsqrt, 40)
            DF$knockoff_stat %>% median
            write.csv(DF, ("results.csv"))
          } else {
            DF = read.csv("results.csv")
          }
          
          # Add qvals and real gene names
          dim(DF)
          DF$q = DF$knockoff_stat %>% rlookc::knockoffQvals(offset = 0)
          DF$Gene1_name = ecoli_anonymization_by_name[DF$Gene1]
          DF$Gene2_name = ecoli_anonymization_by_name[DF$Gene2]
          DF %<>% dplyr::mutate(Gene1_is_decoy = grepl("DECOY", Gene1_name))
          DF %<>% dplyr::mutate(Gene2_is_decoy = grepl("DECOY", Gene2_name))
          
          # How many discoveries do we make?
          {
            pdf(("qval_distribution.pdf"))
            plot(ecdf(DF$q), xlim = 0:1)
            dev.off()
          }
          # Check against gold standards
          check_against_gold_standards(DF)
        }, 
        error=function(cond) {
          cat("Error message received:")
          cat(cond$message)
          return()
        }
        )
      )
      sink()
    })
  })
}

