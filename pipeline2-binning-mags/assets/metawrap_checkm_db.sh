#!/usr/bin/env bash


mkdir -p /fs/project/PAS1117/bioinformatic_tools/checkm1
cd /fs/project/PAS1117/bioinformatic_tools/checkm1

wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz

tar -xzf checkm_data_2015_01_16.tar.gz
rm checkm_data_2015_01_16.tar.gz