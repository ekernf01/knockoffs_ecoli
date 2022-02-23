# This script runs the E coli knockoff GRN experiments assuming
# a blank Ubuntu 18.04 as a starting point.

# For now, I use this interactively to put in passwords etc, but it's almost fully automated.

# Retrieve the ecoli demo repo
git clone https://github.com/ekernf01/knockoffs_ecoli.git
# Run the rest all at once
# source knockoffs_ecoli/run_on_aws.sh

# Install aws cli, build-essential, git, and R v4
echo "Installing software..."
sudo apt-get update
sudo apt-get install -y build-essential
sudo apt-get install -y awscli
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
sudo apt update
sudo apt install -y libssl-dev libcurl4-openssl-dev libxml2-dev
sudo apt install -y r-base r-base-dev r-base-core
R --version #should be 4.1.2 
aws configure

# Retrieve the datasets used.
echo "Fetching data..."
aws s3 cp --recursive s3://cahanlab/eric.kernfeld/datalake/dream5       datalake/dream5
aws s3 cp --recursive s3://cahanlab/eric.kernfeld/datalake/modern_ecoli datalake/modern_ecoli

# Retrieve our package
git clone https://github.com/ekernf01/rlookc.git

# Enter the ecoli demo repo.
cd knockoffs_ecoli

# Change this if you want to run a new set of conditions
mkdir v28 && cd v28

# Install some R packages
Rscript ../dream5_ecoli_install.R
Rscript ../dream5_ecoli_install.R # it needs to refresh for some reason. Don't worry if you see "no package called versions"; that resolves the second time.
# Run tests. 
nohup Rscript ../dream5_ecoli.R 
# Mop up (easy to run this locally too to refine plots)
Rscript ../dream5_ecoli_plots.R
# export results 
aws s3 sync .. s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/dream5_sa_ec/ecoli
# # Then, on laptop:
# aws s3 sync s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v27 v27
