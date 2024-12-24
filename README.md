# plant_smORFs

![Pipeline_image](https://github.com/user-attachments/assets/835587cf-cad3-4385-9779-ce67d2e30abd)

We provide a step-by-step pipeline for the identification of potentially coding smORFs in plant transcriptomic data. 

Our pipeline contains four major steps, which are:
- preprocessing of raw transcriptomic data;
- transcriptome assembly;
- prediction and filtration of smORFs;
- validation using mass-spectrometry data.


Here, we provide the following scripts:
- #### create directories.sh
creates smORF_project with all directories that are necessary for the analysis.


- #### set_working_environment.sh
provides guidelines on how to install software and manage working environment using conda.


- #### mapping.sh
example script that automates mapping using Hisat2 (Kim et al., 2015) for multiple samples.


- #### stringtie_assembly.sh
example script that automates transcript assembly using Stringtie2 (Pertea et al., 2015) for multiple samples.


- #### smorf_filtration.py
filters smORFs predicted by MiPepid based on coding potential and sequence length.

- #### mmseq_filtration.py
filters smORFs based on search against annotated proteins using MMSeqs2 (Steinegger & Söding, 2017).
