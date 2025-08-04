#!/usr/bin/env bash
eval "$(conda shell.bash hook)"
conda activate TLimmuno2

cd TLimmuno2
python Python/TLimmuno2.py --mode file --intdir ../filtered_epitopes.csv --outdir ..