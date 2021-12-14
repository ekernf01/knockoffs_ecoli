# create local user library path (not present by default)
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

# Install packages as they were in late 2021
install.packages("versions", lib = Sys.getenv("R_LIBS_USER"))
library(versions)
versions::install.dates(
  pkgs = c("tidyverse", "magrittr", "glasso", "ggplot2", "knockoff", "irlba", "data.table", "umap"),
  dates = "2021-11-05", 
  lib = Sys.getenv("R_LIBS_USER")
)
# Install our package
install.packages("~/rlookc", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))
