library(GeneNet)
dir.create("genets_benchmark")
X = ecoli_tf_expression %>% sweep(2, colMeans(ecoli_tf_expression), "-")
m.pcor = GeneNet::ggm.estimate.pcor(X)
genenets_fdr_one_target = function(y, X){
  m.pcor = GeneNet::ggm.estimate.pcor(cbind(X,y), lambda = attr(m.pcor, "lambda"), verbose = F)
  sink()
  fdr = GeneNet::network.test.edges(m.pcor, plot = F, verbose = F)
  sink()
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


calibration_nonlinear = rlookc::calibrate__simulateY(
  n_sim = 1000,
  X = X,
  rng_seed = 1,
  alternative_method = genenets_fdr_one_target,
  plot_savepath = ("genets_benchmark/calibration_simulated_y_nonlinear.pdf")
)
saveRDS(calibration_nonlinear, "genets_benchmark/calibration_simulated_y_nonlinear.Rdata")

calibration_linear = rlookc::calibrate__simulateY(
  n_sim = 1000,
  X = X,
  FUN = function(x) x,
  rng_seed = 1,
  alternative_method = genenets_fdr_one_target,
  plot_savepath = ("genets_benchmark/calibration_simulated_y_linear.pdf")
)

saveRDS(calibration_linear, "genets_benchmark/calibration_simulated_y_linear.Rdata")
