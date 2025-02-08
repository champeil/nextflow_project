#!/bin/bash 
# this script is for configuring conda environment
# author: laojp
# time: 2025.02.08
# platform: SYSUCC bioinformatics platform

echo -e "Download the latest version of miniconda"
wget -c -t 0 https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

echo -e "\tInstall miniconda"
# default to accept licence, but not overwrite through -f, and -p to specify the installation path
bash Miniconda3-latest-Linux-x86_64.sh -b 

echo -e "\tinitialize conda"
~/miniconda3/bin/conda init
source ~/.bashrc

echo -e "done"