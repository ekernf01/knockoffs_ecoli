
# Install our package
Rscript -e 'install.packages("rlookc", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))'

# Retrieve the datasets used.
echo "Fetching data..."
mkdir ~/datalake
cd ~/datalake
curl https://zenodo.org/record/6573413/files/modern_ecoli.zip -o modern_ecoli.zip 
curl https://zenodo.org/record/6573413/files/dream5.zip -o dream5.zip
unzip modern_ecoli.zip
unzip dream5.zip
cd /

# Enter the ecoli demo repo.
cd knockoffs_ecoli
mkdir v28
cd v28

# Run tests. 
Rscript ../dream5_ecoli.R 
# Mop up (easy to run this locally too to refine plots)
Rscript ../dream5_ecoli_plots.R