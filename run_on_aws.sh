# This script runs the E coli knockoff GRN experiments assuming
# a blank Ubuntu 18.04 as a starting point. (I used Ubuntu 18 on AWS.)

# Retrieve the ecoli demo repo
git clone https://github.com/ekernf01/knockoffs_ecoli.git

# Install aws cli, build-essential (c compiler), git, curl, and R v4
echo "Installing software..."
sudo apt-get update -qq
sudo apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev libfontconfig1-dev build-essential awscli gdebi-core
R_VERSION=4.1.2
curl -O https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-${R_VERSION}_1_amd64.deb
sudo gdebi r-${R_VERSION}_1_amd64.deb
sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript
R --version # 4.1.2 desired

# Retrieve the datasets used.
echo "Fetching data..."
wget https://zenodo.org/record/6573413/files/modern_ecoli.zip
wget https://zenodo.org/record/6573413/files/dream5.zip
unzip modern_ecoli.zip
unzip dream5.zip

# Retrieve our package
git clone https://github.com/ekernf01/rlookc.git

# Enter the ecoli demo repo.
cd knockoffs_ecoli

# Change this if you want to run a new set of conditions without overwriting. 
mkdir v28 && cd v28

# Install some R packages
Rscript ../dream5_ecoli_install.R
Rscript ../dream5_ecoli_install.R # it needs to refresh for some reason. Don't worry if you see "no package called versions"; that resolves the second time.
# Run tests. 
nohup Rscript ../dream5_ecoli.R 
# Mop up (easy to run this locally too to refine plots)
Rscript ../dream5_ecoli_plots.R

# Optional step: export your results to S3
aws configure #put your credentials
aws s3 sync .. s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/dream5_sa_ec/ecoli
# # Then, on laptop:
# aws s3 sync s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v28 v28
