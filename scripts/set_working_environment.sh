#!/bin/bash

# install conda
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-py312_24.5.0-0-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh

# set conda channels
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment
conda create -n sORF_project
conda activate sORF_project

# install packages (specific versions)

conda install bioconda::fastqc=0.12.1
conda install bioconda::multiqc=1.25.2
conda install bioconda::fastp=0.24.0
conda install bioconda::hisat2
conda install bioconda::stringtie=2.2.3
conda install bioconda::gffread=0.9.12
conda install bioconda::gffcompare=0.12.6
conda install bioconda::mmseqs2

# install samtools from source
mkdir bin
cd bin
wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2
tar -vxjf samtools-1.3.tar.bz2
cd samtools-1.3
make

# set environment for MiPepid
git clone https://github.com/MindAI/MiPepid/

conda deactivate
conda create -n mipepid
conda activate mipepid

conda install scikit-learn=0.19.1
conda install conda-forge::biopython
conda install anaconda::pandas