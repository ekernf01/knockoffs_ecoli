# Retrieve the datasets used.
echo "Fetching data..."
mkdir ~/datalake
cd ~/datalake
wget https://zenodo.org/record/6573413/files/modern_ecoli.zip
wget https://zenodo.org/record/6573413/files/dream5.zip
unzip modern_ecoli.zip
unzip dream5.zip

# At this point, your directory structure should look like this.
#
#
# |-- dream5
# |   |-- DREAM5_network_inference_challenge
# |   |   |-- Evaluation\ scripts
# |   |   |-- Network1
# |   |   |-- Network2
# |   |   |-- Network3
# |   |   |-- Network4
# |   |   `-- README.txt
# |   |-- NIHMS420148-supplement-1.pdf
# |   |-- NIHMS420148-supplement-2.zip
# |   |-- NIHMS420148-supplement-3.xls
# |   |-- NIHMS420148-supplement-5.zip
# |   |-- NIHMS420148-supplement-6.zip
# |   |-- NIHMS420148-supplement-7.xls
# |   |-- NIHMS420148-supplement-8.xls
# |   `-- nihms420148.pdf
# |-- dream5.zip
# |-- modern_ecoli
# |   |-- m3d
# |   |   `-- E_coli_v4_Build_6
# |   |-- regulondb_10.9
# |   |   |-- README.md
# |   |   |-- curated
# |   |   |-- full
# |   |   |-- highthroughput
# |   |   |-- highthroughput2
# |   |   `-- table_chip_controls.csv
# |   `-- transcription_units
# |       |-- All_transcription_units_of_E._coli_K-12_substr._MG1655.txt
# |       `-- README.txt
# `-- modern_ecoli.zip


# Enter the ecoli demo repo.
cd /
cd knockoffs_ecoli
mkdir v35
cd v35
mkdir logs

# Run tests. 
Rscript ../dream5_ecoli.R  &> logs/analysis.txt
# Mop up (easy to run this locally too to refine plots)
Rscript ../dream5_ecoli_plots.R &> logs/plots.txt
