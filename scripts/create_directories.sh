#!/bin/bash

mkdir smORF_project
cd smORF_project

# create directory for reference data (annotation, genome, proteome)
mkdir reference

# create directory for storing RNA-seq data 
mkdir raw_data

# create directories for quality control (qc)
mkdir qc
mkdir qc/fastqc_raw
mkdir qc/fastqc_trimmed

#create directory for trimming
mkdir trimmed_data

#create directory for mapping
mkdir mapping

#create directory for transcriptome annotation
mkdir stringtie

# create directory for MMSeqs2 
mkdir MMSEQ_SEARCH

#create directory for validation with MS data
mkdir MS_data


