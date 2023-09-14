# create local user library path (not present by default)
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
# Install packages as they were in late 2021
install.packages("https://cran.r-project.org/src/contrib/Archive/remotes/remotes_2.4.2.tar.gz", 
                 lib = Sys.getenv("R_LIBS_USER"),
                 source = TRUE,
                 repos = NULL)
pv  = read.table(header = 1, text = 
                   "package version
tidyverse 1.3.1
magrittr 2.0.1
glasso 1.11
ggplot2 3.3.5
knockoff 0.3.3
irlba 2.3.3
data.table 1.14.2
tsne 0.1.3
FNN 1.1.3
BiocManager 1.30.16
poolr 1.0.0
svglite 2.0.0
vita 1.0.0
GeneNet 1.2.15
mclust 5.4.7"
)
for(i in rownames(pv)){
  remotes::install_version(
    pv[i, "package"], 
    version = pv[i, "version"], 
    lib = Sys.getenv("R_LIBS_USER"), 
    upgrade = "never", 
    quiet = TRUE, 
    repos = "https://cloud.r-project.org/"
  )
}
BiocManager::install(c("limma", "GENIE3", "doRNG"), version = "3.14")

# Install our package
install.packages("~/rlookc", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))

# Install BINCO
install.packages("https://cran.r-project.org/src/contrib/Archive/BINCO/BINCO_0.1-1.tar.gz", from = "source", repo = NULL)
