FROM ubuntu:20.04

# Install aws cli, build-essential (c compiler), curl, git, and R v4
RUN apt-get update -qq
RUN apt-get install -yq libicu-dev 
RUN apt-get install -yq libc6-dev 
RUN apt-get install -yq libc6 
RUN apt-get install -yq libfreetype6-dev 
RUN apt-get install -yq libexpat1-dev 
RUN apt-get install -yq libssl-dev 
RUN apt-get install -yq libcurl4-openssl-dev 
RUN apt-get install -yq libxml2-dev 
RUN apt-get install -yq libfontconfig1-dev 
RUN apt-get install -yq build-essential 
RUN apt-get install -yq awscli 
RUN apt-get install -yq gdebi-core 
RUN apt-get install -yq git
RUN apt-get install -yq curl
RUN R_VERSION=4.1.2 && curl -O https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-${R_VERSION}_1_amd64.deb 
RUN R_VERSION=4.1.2 && gdebi -n r-${R_VERSION}_1_amd64.deb
RUN R_VERSION=4.1.2 && ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R && ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript
RUN R --version


# Get our code
RUN git clone https://github.com/ekernf01/knockoffs_ecoli.git
RUN git clone https://github.com/ekernf01/rlookc.git

# Install the R packages we need
RUN Rscript knockoffs_ecoli/dream5_ecoli_install.R
RUN Rscript knockoffs_ecoli/dream5_ecoli_install.R


