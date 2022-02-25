source("../dream5_ecoli_setup.R")

make_knockoffs = function(ecoli_tf_expression_working_copy, confounder_indicators, knockoff_type){
  X = cbind(ecoli_tf_expression_working_copy, confounder_indicators)
  if(knockoff_type == "mixture"){
    ecoli_tf_expression_working_copy_knockoffs =
      rlookc::computeGaussianMixtureKnockoffs(
        X,
        do_high_dimensional = T,
        hard_assignments = clusters$cluster
      )
  } else if(knockoff_type == "naive"){
    # Assumes variables are in columns! Otherwise you're
    # scrambling all genes within a sample, not all samples within a gene.
    ecoli_tf_expression_working_copy_knockoffs = apply(X, 2, sample, replace = F)
  } else if (knockoff_type=="sample"){
    Sigma = cov(X)
  } else if (knockoff_type == "shrinkage"){
    if(is.na(shrinkage_param)){
      off_diagonal_shrinkage = corpcor::estimate.lambda(X)
    } else {
      off_diagonal_shrinkage = shrinkage_param
    }
    Sigma = as(corpcor::cor.shrink(X), "matrix")
  } else if (knockoff_type == "glasso"){
    glasso_out = glasso::glasso(cov(X), rho = shrinkage_param, penalize.diagonal = F)
    Sigma = glasso_out$w
    write.table( mean(glasso_out$wi!=0), ("glasso_wi_fraction_nonzero.csv"))
  }
  
  if(knockoff_type %in% c("glasso", "sample", "shrinkage")){
    # These methods differ only in the covariance estimate Sigma, done above, but
    # otherwise they're the same.
    ecoli_tf_expression_working_copy_knockoffs = rlookc::computeGaussianKnockoffs(X, Sigma = Sigma)
    # For TF targets, use fast leave-one-out knockoffs (LOOKs).
    # This output contains all of them in a compact form.
    # Realize them later on using rlookc::formOneLook.
    looks_compact = rlookc::generateLooks(
      X,
      mu = colMeans(X),
      Sigma = Sigma,
      statistic = fast_lasso_penalty,
      output_type = "knockoffs_compact",
      vars_to_omit = seq(ncol(ecoli_tf_expression_working_copy))
    )
  } else if(knockoff_type %in% c("naive", "mixture")){
    # Sometimes, efficient code to make LOOKs is not available.
    # So that the code downstream can remain the same, I pretend to make
    # LOOKs inside the same data structure used above, even when they just contain
    # the full knockoffs unmodified.
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
  } else {
    warning("Knockoff type not recognized.\n")
    break
  }
  # If confounders were included in knockoff generation, strip them out.
  ecoli_tf_expression_working_copy_knockoffs = ecoli_tf_expression_working_copy_knockoffs[,1:334]
  list(ecoli_tf_expression_working_copy_knockoffs, looks_compact)
}

dir.create("logs")
do_one = function(condition_index, omit_knockout_evaluation_samples = F, reuse_results = F){
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
      try({
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
        
        
        # Make knockoffs. Check fit.
        cat("\nMaking knockoffs.\n")
        l = make_knockoffs(ecoli_tf_expression_working_copy = ecoli_tf_expression_working_copy,
                           confounder_indicators = confounder_indicators,
                           knockoff_type = knockoff_type)
        ecoli_tf_expression_working_copy_knockoffs = l[[1]]
        looks_compact = l[[2]]
        rm(l)
        write.csv(ecoli_tf_expression_working_copy_knockoffs, "knockoffs.csv")
        
        # Check power: how different are knockoffs from original features?
        cor_with_knockoffs = cor(ecoli_tf_expression_working_copy_knockoffs, ecoli_tf_expression_working_copy) %>% diag
        data.frame( cor_with_knockoffs,
                    is_decoy = grepl("DECOY", ecoli_anonymization_by_name[names(cor_with_knockoffs)])) %>%
          ggplot() +
          geom_histogram(aes(x = cor_with_knockoffs, fill = is_decoy, bins = 30)) +
          ggtitle("Correlation of each TF with its knockoff")
        ggsave(("correlations_with_knockoffs.pdf"), width = 5, height = 3)
        
        # Check calibration for batch adjustment.
        # No association should be detected with factors explicitly conditioned on.
        if(condition_on!="none"){
          confounder_association_qvals =
            rlookc::knockoffQvals(
              cor(confounder_indicators, ecoli_tf_expression_working_copy) -
                cor(confounder_indicators, ecoli_tf_expression_working_copy_knockoffs)
            )
          pdf(("confounder_association_qvals.pdf"), width = 5, height = 2.5)
          hist(confounder_association_qvals, 40, xlim = 0:1)
          dev.off()
        }
        
        # Check calibration
        cat("\nChecking calibration\n")
        calibration_typical = rlookc::simulateY(
          n_sim = 100,
          X = ecoli_tf_expression_working_copy %>% sweep(2, colMeans(ecoli_tf_expression_working_copy), "-"),
          knockoffs = ecoli_tf_expression_working_copy_knockoffs %>% sweep(2, colMeans(ecoli_tf_expression_working_copy), "-") ,
          statistic = fast_lasso_penalty,
          plot_savepath = ("average_case_calibration.pdf")
        )
        write.csv(calibration_typical$calibration$fdr, ("average_case_calibration.csv"))
        
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
                                 knockoff_type = knockoff_type)
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
      })
    )
    sink()
  })
}


# # For interactive use, if you only update the gold standards and not the knockoff procedure, 
# # you can re-run this with "reuse_results = T".
parallel::mclapply(nrow(conditions):1, do_one, omit_knockout_evaluation_samples = F, reuse_results = F, mc.cores = parallel::detectCores())
parallel::mclapply(nrow(conditions):1, do_one, omit_knockout_evaluation_samples = T, reuse_results = F, mc.cores = parallel::detectCores())

# This is for internal use to experiment fast on a particularly important case.
# do_one(34, omit_knockout_evaluation_samples = T, reuse_results = T)
