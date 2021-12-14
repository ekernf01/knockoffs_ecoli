# This script runs the E coli knockoff GRN experiments assuming
# a blank Ubuntu 18.04 as a starting point.

# For now, I use this interactively to put in passwords etc, but it's almost fully automated.

# Retrieve the ecoli demo repo
git clone https://ekernf01@bitbucket.org/ekernf01/transcriptome_knockoffs.git

# Install aws cli, build-essential, git, and R v4
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
aws s3 cp --recursive s3://cahanlab/eric.kernfeld/datalake/dream5       datalake/dream5
aws s3 cp --recursive s3://cahanlab/eric.kernfeld/datalake/modern_ecoli datalake/modern_ecoli

# Retrieve our package
<<<<<<< HEAD
git clone https://github.com/ekernf01/rlookc.git
=======
git clone https://ekernf01@bitbucket.org/ekernf01/rlookc.git
>>>>>>> cfd97fdcac23b2bf0a52cf5c2362563ffde59565

# Enter the ecoli demo repo.
cd transcriptome_knockoffs/applications/dream5_sa_ec/ecoli/

# Change this if you want to run a new set of conditions
mkdir v26 && cd v26

# Install some R packages
Rscript ../dream5_ecoli_install.R
Rscript ../dream5_ecoli_install.R # it needs to refresh for some reason

# Run tests. 
<<<<<<< HEAD
nohup Rscript ../dream5_ecoli.R 

# Mop up (easy to run this locally too to refine plots)
=======
nohup Rscript ../dream5_ecoli.R &

# Mop up (easy to run this locally too)
>>>>>>> cfd97fdcac23b2bf0a52cf5c2362563ffde59565
Rscript ../dream5_ecoli_plots.R

# export results 
aws s3 sync .. s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/dream5_sa_ec/ecoli

# # Then, on laptop:
<<<<<<< HEAD
# mkdir v26 && cd v26
=======
# mkdir v24 && cd v24
>>>>>>> cfd97fdcac23b2bf0a52cf5c2362563ffde59565
# aws s3 sync s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/dream5_sa_ec/ecoli/v24 .
