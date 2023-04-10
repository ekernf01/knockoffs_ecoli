# create local user library path (not present by default)
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

# Install packages as they were in late 2021
install.packages("versions", lib = Sys.getenv("R_LIBS_USER"))
library(versions)
versions::install.dates(
  pkgs = c("tidyverse", "magrittr", "glasso", "ggplot2", "knockoff", "irlba", "data.table",
           "tsne", "FNN", "BiocManager", "poolr", "svglite", "GeneNet", "mclust"),
  dates = "2021-11-05",
  lib = Sys.getenv("R_LIBS_USER")
)


BiocManager::install("limma", version = "3.14")
# Install our package
install.packages("~/rlookc", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))

# Install BINCO
install.packages("https://cran.r-project.org/src/contrib/Archive/BINCO/BINCO_0.1-1.tar.gz", from = "source", repo = NULL)