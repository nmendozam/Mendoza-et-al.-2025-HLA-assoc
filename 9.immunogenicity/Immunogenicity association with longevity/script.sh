# #!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"
conda activate allele_tools

# INPUT="1.IEDB_peptides_len13-25_DRB1 - v2.xlsx"
# INPUT="1.IEDB_peptides_len13-25_DRB1_source_human.xlsx"
INPUT="1.IEDB_peptides_len13-25_DRB1_source_human + negative.xlsx"

# Extract only important epitopes
python 00_Input/extract_epitopes.py "00_Input/$INPUT" "01_preparing files/Selection of epitopes.xlsx"

# pre-process epitopes for TLimmuno2
python "01_preparing files/pre-process.py" "01_preparing files/Selection of epitopes.xlsx"  "01_preparing files/HLA-seq_allele_freq.csv" "02_immunogenicity_prediction/immunogenicity_input.csv"

# run immunogenicity predictions
conda activate TLimmuno2
cd 02_immunogenicity_prediction/TLimmuno2/
echo $PWD
# Important note for the command below: The library paths are added to enable the use
# of CUDA nodes in the cau cluster. TLimmuno is a custom conda environment that goes
# beyond what was indicated in the initial instructions of the package, again to
# enable CUDA nodes.
    srun -p gpu --gpus-per-node=2 --time=02:00:00 --mem=64G --cpus-per-task=8 \
        python Python/TLimmuno2.py\
            --mode file\
            --intdir "../immunogenicity_input.csv"\
            --outdir ".."
cd ../..
conda activate allele_tools
 
# post-process predictions for the statistic model
python 02_immunogenicity_prediction/post-process.py "02_immunogenicity_prediction/result.csv" "02_immunogenicity_prediction/Ageing.2field.pyhla" "03_logistic_model/DRB1"  


# # Run the statistic model
# # Rscript  "03_logistic_model/DRB1.csv"  "03_logistic_model/DRB1_1.csv"


mkdir -p "output/$INPUT"
cp 03_logistic_model/DRB1* "output/$INPUT/"
cp "01_preparing files/Selection of epitopes.xlsx" "output/$INPUT/"
