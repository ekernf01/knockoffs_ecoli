# This script defines a list "ecoli_networks" containing our gold standards as dataframes.
# It also defines functions compute_calibration_with_unknowns and check_against_gold_standards.
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

# Locate data
# for eric laptop
DATALAKE = "~/Desktop/jhu/research/datalake"
# for AWS
if(!dir.exists(DATALAKE)){
  DATALAKE = "~/datalake"
}
# otherwise, complain and quit.
if(!dir.exists(DATALAKE)){
  stop("Datalake not found. Place it in '~/datalake' or modify `dream5_ecoli_setup.R`.\n")
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
    fitTrtMean <- limma::lmFit(current_expression, design.trt)
    fit.contrast=limma::contrasts.fit(fitTrtMean, coefficients = colnames(design.trt)[-1])
    efit.contrast=limma::eBayes(fit.contrast)
    ecoli_networks$M3Dknockout[[as.character(e)]] =
      efit.contrast$p.value %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Gene2_name") %>%
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
  # The ihfA antibody captures ihfB too
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
  
  cat("Gold standard size:\n")
  ecoli_networks$chip_intersection_M3Dknockout %>%
    extract(c("Gene1_name", "is_confirmed")) %>%
    table %>%
    print
}

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
      # Computing fdr and then subsetting is not the same as subsetting and then computing fdr.
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

