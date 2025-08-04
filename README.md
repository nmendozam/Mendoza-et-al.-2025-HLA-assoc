# <p align="center"> HLA-DRB1*15:01 is associated with a reduced likelihood of longevity in northern European men </p>

This repository contains the code and data used for the analysis of the paper with the same title. Here, an association test of imputed HLA alleles revealed a sex-specific association with HLA-DRB1*15:01

## Requirements
Long processes or highly computationally intensive tasks were run on the HPC cluster at the Institute of Clinical Molecular Biology (IKMB) in Kiel, Germany. Therefore, the scripts with the extension `.job` are SLURM job scripts and can be run on an HPC cluster by executing `sbatch file.job`, but they can also be run locally by executing `bash file.job`. The scripts are designed to be run in a Linux environment, and they assume that the necessary software and dependencies are installed on the system. To run the code, you will need the following software:
- Python 3.8 or higher
- R 4.5.1
- plink 1.9

And the following Python packages:
- [alleleTools](https://github.com/nmendozam/alleleTools)
- [PyHLA](https://github.com/felixfan/PyHLA)
- [HLA-TAPAS](https://github.com/immunogenomics/HLA-TAPAS)
- pandas==2.2.3
- numpy==2.3.1
- scikit-learn==1.16.0

## Repository Structure
Each folder corresponds to a specific analysis step or data processing task. The main folders are:

 - `1.QC/` contains the quality control script.
 - `2.Imputation/` has the imputation of HLA alleles using HLA-TAPAS and 1000 Genomes reference panel.
 - `3.association/` has the association tests performed with PyHLA
 - `4.SNPassociation/` contains the association of the MS related variant rs9267649 with longevity.
 - `5.Danish/`
 - `6.UKBiobank/` includes the replication in the UK Biobank
 - `7.Validation/` is the comparison of HLA-seq genotyping with imputed HLA alleles to assess the reliability of the imputation approach.
 - `8.ukb_apoe_n_hla/` interaction analysis of HLA-DRB1*15:01 with APOE
 - `9.immunogenicity/` contains the analysis of immunogenicity association with longevity.