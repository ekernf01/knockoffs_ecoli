# This script defines our E coli experiments.

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

GAUSSIAN_KNOCKOFF_STRATEGIES = c("glasso", "sample", "shrinkage", "bagging")

# We'll run several experiments with different handling of perturbations,
# confounding influences, and knockoff construction.
{
  conditions = rbind(
    # study non-knockoff methods
    expand.grid(
      knockoff_type = c("NA", "NA"),
      fdr_control_method = c("gaussian_mirror", "GeneNet"),
      address_genetic_perturbations = c( F ),
      condition_on = c( "none"),
      shrinkage_param = NA,
      do_simulate_y = T,
      seed = 1,
      omit_knockout_evaluation_samples = c(T, F)
    ),
    # study knockoff generation methods
    expand.grid(
      knockoff_type = c("sample", "permuted", "shrinkage", "mixture"),
      fdr_control_method = c("knockoffs", "knockoffs", "knockoffs", "knockoffs"),
      address_genetic_perturbations = c( F ),
      condition_on = c( "none"),
      shrinkage_param = NA,
      do_simulate_y = T,
      seed = 1,
      omit_knockout_evaluation_samples = c(T, F)
    ),
    # study another knockoff generation method with different tuning params
    expand.grid(
      fdr_control_method = "knockoffs",
      knockoff_type = c( "glasso"),
      address_genetic_perturbations = c( F ),
      condition_on = c( "none"),
      shrinkage_param = c(0.01, 0.001, 0.0001),
      do_simulate_y = T,
      seed = 1,
      omit_knockout_evaluation_samples = c(T, F)
    ),
    # study random variation in knockoff construction
    expand.grid(
      fdr_control_method = c("knockoffs"),
      knockoff_type = c("sample"),
      address_genetic_perturbations = c( F ),
      condition_on = c( "none" ),
      shrinkage_param = NA,
      do_simulate_y = F,
      seed = 1:5,
      omit_knockout_evaluation_samples = c(T, F)
    ),
    # study confounder handling
    expand.grid(
      fdr_control_method = "knockoffs",
      knockoff_type = c("glasso"),
      address_genetic_perturbations = c( F, T ),
      condition_on = c(  "none",
                         "pert_labels",
                         "pert_labels_plus_pca10",
                         "pert_labels_plus_pca20",
                         "pert_labels_plus_pca30",
                         "pert_labels_plus_pca50" ),
      shrinkage_param = 0.0001,
      do_simulate_y = F,
      seed = 1,
      omit_knockout_evaluation_samples = c(T, F)
    )
  )
  # Save time by not re-doing settings shared across experiments.
  # change the seed explicitly if you want replicates.
  conditions = conditions[!duplicated(conditions),]
  write.csv(conditions, "experiments_to_run.csv")
}


# Compute just part of the lasso path for faster results
dfmax = 21
fast_lasso_penalty = function(X, X_k, y) {
  suppressWarnings(
    knockoff::stat.lasso_lambdasmax(
      y   = y,
      X   = X,
      X_k = X_k,
      dfmax = dfmax
    )
  )
}

# Check out the expression data briefly
ecoli_expression[1:4, 1:4]
dim(ecoli_expression)
ecoli_tf_expression = ecoli_expression[ecoli_tf[[1]]] %>% as.matrix
clusters = kmeans(ecoli_expression, centers = 7)
table(clusters$cluster)
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

#' This handles a difference between a causal DAG and the corresponding MRF structure (the moral graph).
#' Edges in the MRF can occur between parents & children (good) or spouses (bad), For us, all TF-TF pairs
#' are potential spouses, so we have the option to omit them from benchmarks.
#'
omit_possible_spouses = function(DF_plus_gs){
  subset(DF_plus_gs, !(Gene1 %in% ecoli_tf) | !(Gene2 %in% ecoli_tf)  )
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
#' Usually called via make_knockoffs, not directly. For 'permuted' and 'mixture' methods, returns NULL.
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
  } else if (method %in% c("mixture", "permuted")) {
    return(NULL)
  } else {
    stop(paste0("Invalid 'method' arg ", method, " in estimate_covariance."))
  }
  return(Sigma)
}

#' Set knockoffs equal to original data for the most difficult observations. Then re-fit.
#'
#' @param X data
#' @param ... passed to estimate_covariance
#'
#' Usually called via make_knockoffs, not directly.
#'
make_gaussian_knockoffs = function(X, method, Sigma = NULL){
  if(is.null(Sigma)){
    Sigma = estimate_covariance(X, method = method)
  }
  Xin = X
  Sigma_in = Sigma
  looks_compact = rlookc::generateLooks(
    Xin,
    mu = colMeans(Xin),
    Sigma = Sigma_in,
    statistic = fast_lasso_penalty,
    output_type = "knockoffs_compact",
    vars_to_omit = seq(ncol(X))
  )
  return(looks_compact)
}

#' This is the main user-facing code to make knockoffs during each experiment.
#'
#' Takes care of correcting for PC's and making knockoffs under all different types of 
#' models and covariance estimates.
#'
make_knockoffs = function(ecoli_tf_expression_working_copy, confounder_indicators, knockoff_type, Sigma = NULL){
  X = cbind(ecoli_tf_expression_working_copy, confounder_indicators)
  if(knockoff_type == "mixture"){
    ecoli_tf_expression_working_copy_knockoffs =
      rlookc::computeGaussianMixtureKnockoffs(
        X,
        do_high_dimensional = T,
        hard_assignments = clusters$cluster
      )
  } else if(knockoff_type == "permuted"){
    ecoli_tf_expression_working_copy_knockoffs = apply(X, 2, sample, replace = F)
  } else if(knockoff_type %in% GAUSSIAN_KNOCKOFF_STRATEGIES){
    # These methods differ only in the covariance estimate Sigma, done above, but
    # otherwise they're the same.
    looks_compact = make_gaussian_knockoffs(X, method = knockoff_type, Sigma = Sigma)
    ecoli_tf_expression_working_copy_knockoffs = looks_compact$knockoffs
  } else {
    stop("Knockoff type not recognized.\n")
  }
  # If confounders were included in knockoff generation, strip them out.
  ecoli_tf_expression_working_copy_knockoffs = ecoli_tf_expression_working_copy_knockoffs[,1:334]
  return(ecoli_tf_expression_working_copy_knockoffs)
}

#' Run one experiment.
#'
#' One experiment, e.g. "gaussian knockoffs with glasso 0.001 covariance estimate, correcting for 20 principal components
#' and with no special handling of genetic perturbations."
do_one = function(condition_index, reuse_results = F, test_mode = F){
  cat("\n======================================\n")
  cat("Experiment", condition_index, "of", nrow(conditions), '\n')
  cat(  "======================================\n")
  attach(conditions[condition_index,], warn.conflicts = F)
  tryCatch(set.seed(seed), error = function(e){warning("Random seed not found and not set.\n")})
  
  # Set up the data
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
  
  {
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
    
    if(fdr_control_method=="knockoffs"){
      # Make knockoffs.
      cat("\nConstructing knockoffs...\n")
      ecoli_tf_expression_working_copy_knockoffs = make_knockoffs(ecoli_tf_expression_working_copy = ecoli_tf_expression_working_copy,
                                                                  confounder_indicators = confounder_indicators,
                                                                  knockoff_type = knockoff_type)
      print("... knockoffs ready.")
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
    }  
    # Check calibration with simulated Y
    if(do_simulate_y){
      cat("\nChecking calibration\n")
      if(fdr_control_method=="knockoffs"){
        Sigma = estimate_covariance(cbind(ecoli_tf_expression_working_copy, confounder_indicators), method = knockoff_type)
        calibration_typical = rlookc::calibrate__simulateY(
          cores = 1,
          rng_seed = seed,
          n_sim = ifelse(test_mode, 5, 1000),
          active_set_size = 1,#pmax(1, rpois(2, n = 1000)),
          X = ecoli_tf_expression_working_copy %>% sweep(2, colMeans(ecoli_tf_expression_working_copy), "-"),
          knockoffs = replicate(10, simplify = F, {
            make_knockoffs(ecoli_tf_expression_working_copy = ecoli_tf_expression_working_copy,
                           confounder_indicators = confounder_indicators,
                           knockoff_type = knockoff_type,
                           Sigma = Sigma) %>%
              sweep(2, colMeans(ecoli_tf_expression_working_copy), "-")
          }),
          statistic = fast_lasso_penalty,
          plot_savepath = ("average_case_calibration.pdf")
        )
      } else {
        if(fdr_control_method=="gaussian_mirror"){
          library("parallel")
          alternative_method = function(y, X) {
            rlookc::calibrate__getQvals(GM::gm(
              y=y, 
              x=X, 
              do_simultaneous = T, 
              ncores = 15
            )$gm_statistics)
          }
        } else if(fdr_control_method=="GeneNet"){
          m.pcor = GeneNet::ggm.estimate.pcor(ecoli_tf_expression_working_copy)
          alternative_method = function(y, X){
            m.pcor = GeneNet::ggm.estimate.pcor(cbind(y,X), lambda = attr(m.pcor, "lambda"), verbose = F)
            fdr = fdrtool::fdrtool(m.pcor[1,-1], statistic = "correlation", verbose = F) %>%
              extract2("qval")
            y_idx = ncol(X)+1
            # fdr is not preserved under arbitrary subsetting, but local fdr (1 - posterior prob) is.
            # So, after subsetting, recompute fdr as cumulative average of local fdr.
            fdr %>% 
              subset( (node2==y_idx)) %>% 
              dplyr::arrange(-prob) %>%
              dplyr::mutate(qval = dplyr::cummean(1-prob)) %>%
              dplyr::arrange(node1) %>% 
              extract2("qval") 
          }
        } else {
          stop("fdr_control_method must be 'knockoffs' or 'GeneNet' or 'gaussian_mirror'.")
        }
        calibration_typical = rlookc::calibrate__simulateY(
          cores = 1,
          rng_seed = seed,
          n_sim = ifelse(test_mode, 5, 1000),
          active_set_size = pmax(1, rpois(2, n = 1000)),
          X = ecoli_tf_expression_working_copy %>% sweep(2, colMeans(ecoli_tf_expression_working_copy), "-"),
          alternative_method = alternative_method,
          plot_savepath = ("average_case_calibration.pdf")
        )
      }
      write.csv(calibration_typical$calibration$fdr, ("average_case_calibration.csv"))
      write.csv(calibration_typical$calibration$p_true_discoveries, ("average_case_power.csv"))
    }
    # Deploy on non-TF targets only, to avoid "spouse" problems as described in the paper
    if(!reuse_results){
      cat("\nComputing regulators for all non-TF genes. This will take a long time.\n")
      target_genes = colnames(ecoli_expression_working_copy) %>% setdiff(ecoli_tf[[1]])
      w = list()
      if(test_mode){target_genes = target_genes[1:4]}
      for(i in seq_along(target_genes)){
        cat(".")
        target = target_genes[[i]]
        samples_to_include = get_samples_to_include(target)
        if(fdr_control_method=="gaussian_mirror"){
          library("parallel")
          w[[target]] = GM::gm(
            y   = ecoli_expression_working_copy[samples_to_include,target],
            x   = ecoli_tf_expression_working_copy[samples_to_include,],
            do_simultaneous = T, 
            ncores = 15
          )$gm_statistics
        } else if(fdr_control_method=="GeneNet"){
          y   = ecoli_expression_working_copy[samples_to_include,target]
          X   = ecoli_tf_expression_working_copy[samples_to_include,]
          w[[target]] = GeneNet::ggm.estimate.pcor(cbind(y,X), verbose = F)[1, -1]
          rm(y)
          rm(X)
        } else if(fdr_control_method=="knockoffs"){
          if(i %% 100 == 0){
            cat( "\n  Target", i, "reached. Remaking knockoffs. ")
            ecoli_tf_expression_working_copy_knockoffs = 
              make_knockoffs(
                ecoli_tf_expression_working_copy = ecoli_tf_expression_working_copy,
                confounder_indicators = confounder_indicators,
                knockoff_type = knockoff_type
              )
          }
          w[[target]] = fast_lasso_penalty(
            y   = ecoli_expression_working_copy[samples_to_include,target],
            X   = ecoli_tf_expression_working_copy[samples_to_include,],
            X_k = ecoli_tf_expression_working_copy_knockoffs[samples_to_include,]
          )
        }
      }
      # Assemble results nicely
      DF = list()
      for(target in target_genes){
        DF[[target]] = data.frame(
          Gene1 = colnames(ecoli_tf_expression_working_copy),
          Gene2 = target,
          knockoff_stat = w[[target]]
        )
      }
      DF = data.table::rbindlist(DF)
      signedsqrt = function(x) sign(x)*sqrt(abs(x))
      hist(DF$knockoff_stat %>% signedsqrt, 40)
      DF$knockoff_stat %>% median
      write.csv(DF, ("results.csv"))
    } else {
      DF = read.csv("results.csv")
    }
    
    # Add qvals and real gene names
    dim(DF)
    
    if(fdr_control_method=="gaussian_mirror"){
      DF$q = DF$knockoff_stat %>% rlookc::knockoffQvals(offset = 0)
    } else if(fdr_control_method=="GeneNet"){
      DF$q = DF$knockoff_stat %>%       
        fdrtool::fdrtool(statistic = "correlation", verbose = F) %>%
        extract2("qval")
    } else if(fdr_control_method=="knockoffs"){
      DF$q = DF$knockoff_stat %>% rlookc::knockoffQvals(offset = 0)
    }
    DF$Gene1_name = ecoli_anonymization_by_name[DF$Gene1]
    DF$Gene2_name = ecoli_anonymization_by_name[DF$Gene2]
    DF %<>% dplyr::mutate(Gene1_is_decoy = grepl("DECOY", Gene1_name))
    DF %<>% dplyr::mutate(Gene2_is_decoy = grepl("DECOY", Gene2_name))
    
    # How many discoveries do we make?
    try({
      pdf(("qval_distribution.pdf"))
      plot(ecdf(DF$q), xlim = 0:1)
      dev.off()
    })
    # Check against gold standards
    check_against_gold_standards(DF)
    
  }
  
}

do_one_safe = function(condition_index, reuse_results = F, test_mode = F){
  stopifnot("'condition_index' is bigger than 'nrow(conditions)'."=condition_index %in% seq_along(conditions[[1]]))
  # Make the output folder
  working_directory = ""
  for(j in seq_along(conditions[condition_index,])){
    prefix = paste0( colnames(conditions)[[j]], "=", conditions[condition_index, j] )
    working_directory %<>% file.path(prefix, .)
  }
  dir.create(working_directory, recursive = T, showWarnings = F)
  cat("\nTesting condition:", working_directory, "\n")
  withr::with_message_sink(
    file.path("logs", condition_index), 
    {
      withr::with_output_sink(
        file.path("logs", condition_index),
        {
          withr::with_dir(
            working_directory,
            code = {
              
              tryCatch( 
                do_one(condition_index, reuse_results = F, test_mode = F),
                error=function(cond) {
                  cat("Error message received:")
                  cat(cond$message)
                  return()
                }
              )
            }
          )
        }
      )
    }
  )
}

