#!/bin/bash
module add gcc12-env/12.1.0 && module add miniconda3

OUT_FATHER="father_test_with_pc.white.assoc"
OUT_MOTHER="mother_test_with_pc.white.assoc"

# Generate phenotype files for both parents
echo "Generating phenotype files..."
conda activate allele_tools
cd phenotype
python main.py

cd ..
# Integrate phenotype in an allele table
# for fathers
altools convert ukb2allele HLA_imputed_data.csv \
    --phenotype phenotype/father_lli.with_ethnicity.csv \
    --remove_pheno_zero \
    --output output/father_lli.ethnicity.alt
# for mothers
altools convert ukb2allele HLA_imputed_data.csv \
    --phenotype phenotype/mother_lli.with_ethnicity.csv \
    --remove_pheno_zero \
    --output output/mother_lli.ethnicity.alt

# Perform association
# for fathers
conda activate py_hla
python PyHLA/ukb_pc.py \
    --input output/father_lli.ethnicity.alt \
    --summary -d 4 \
    --out output/father_lli.ethnicity.summary
awk '/^Population/ { matched = 1 } matched' output/father_lli.ethnicity.summary \
    | awk '/\*/ {print $1}' \
    | grep -v "DRB1\*15:01" > output/exclude_alleles.fathers.txt
srun -c 16 --mem=64g --time=20:00 \
    python PyHLA/PyHLA.py \
        --input output/father_lli.ethnicity.alt \
        --assoc \
        --out $OUT_FATHER \
        --covar ukb_pc.tsv \
        --test logistic \
        --model additive \
        --exclude output/exclude_alleles.fathers.txt \
        --covar-name pc_1,pc_2,pc_3
echo "Association result long-lived fathers"
cat <(head -n 1 $OUT_FATHER) <(grep "DRB1.*15:01" $OUT_FATHER)
# for mothers
conda activate py_hla
conda activate py_hla
python PyHLA/ukb_pc.py \
    --input output/mother_lli.ethnicity.alt \
    --summary -d 4 \
    --out output/mother_lli.ethnicity.summary
awk '/^Population/ { matched = 1 } matched' output/mother_lli.ethnicity.summary \
    | awk '/\*/ {print $1}' \
    | grep -v "DRB1\*15:01" > output/exclude_alleles.mothers.txt
srun -c 16 --mem=64g --time=20:00 \
    python PyHLA/PyHLA.py \
        --input output/mother_lli.ethnicity.alt \
        --assoc \
        --out $OUT_MOTHER \
        --covar ukb_pc.tsv \
        --test logistic \
        --model additive \
        --exclude output/exclude_alleles.mothers.txt \
        --covar-name pc_1,pc_2,pc_3
echo "Association result long-lived mothers"
cat <(head -n 1 $OUT_MOTHER) <(grep "DRB1.*15:01" $OUT_MOTHER)

