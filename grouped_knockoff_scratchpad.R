# This is a work in progress, hence commented out. Sorry. 
# Group TFs whose effects you can't distinguish.
# n_groups = 1
# set.seed(0)
# gene_groups = kmeans(
#   centers = n_groups, 
#   ecoli_tf_expression %>%
#     sweep(2, colMeans(ecoli_tf_expression), "-") %>%
#     sweep(2, apply(ecoli_tf_expression, 2, sd),  "/") %>%
#     t
# )$cluster
# names(gene_groups) = colnames(ecoli_tf_expression)

n_groups = 334
gene_groups = seq_along(ecoli_tf_expression)

# Make knockoffs. Check fit.
{
  X = cbind(ecoli_tf_expression, confounder_indicators)
  # Put all confounders in one group; they'll be discarded anyway.
  groups_including_confounders = gene_groups %>% c(rep(n_groups + 1, ncol(confounder_indicators)))
  # The rlookc code wants a list of variable indices with length equal to the number of groups.
  # Not a list of group indices with length equal to the number of variables!
  # Sorry.
  groups_reformatted = lapply(seq(n_groups+1), function(g) which(g==groups_including_confounders))
  t1 = Sys.time()
  ecoli_tf_expression_knockoffs = rlookc::computeGaussianKnockoffs(X)
  t2 = Sys.time()
  print(t2 - t1)
  # For TF targets, use fast leave-one-out knockoffs (LOOKs).
  # This output contains all of them in a compact form.
  # Realize them later on using rlookc::formOneLook.
  looks_compact = rlookc::generateLooks(
    X,
    mu = colMeans(X),
    Sigma = cov(X),
    statistic = knockoff::stat.lasso_lambdasmax, 
    output_type = "knockoffs_compact", 
    vars_to_omit = seq(ncol(ecoli_tf_expression))
  )  
  t3 = Sys.time()
  print(t3 - t2)
  
}
ecoli_tf_expression_knockoffs = ecoli_tf_expression_knockoffs[,1:334]
