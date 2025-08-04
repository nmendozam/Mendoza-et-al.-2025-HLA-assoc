#!/usr/bin/env python3
import argparse
import re

import pandas as pd

parser = argparse.ArgumentParser(description='Extract the desired epitopes from IEDB export')
parser.add_argument(
        'input',
        type=str,
        help='Excel file exported from IEDB with single header, containing the epitopes.',
)
parser.add_argument(
        'output',
        type=str,
        help='Output file in excel format "xlsx"',
)
args = parser.parse_args()


df = pd.read_excel(args.input, sheet_name="Sheet1")
df = df.dropna(subset=["Epitope - Source Molecule"])

desired_epitopes = [
    "Myelin basic protein",
    "Calreticulin",
    "Apolipoprotein B-100",
    "Cofilin-1",

]
search_patterns = [re.escape(m) for m in desired_epitopes]

name_match = df["Epitope - Molecule Parent"].str.contains('|'.join(search_patterns), case=False).fillna(False)


subset = df.loc[
        name_match,
        ["Epitope - Source Molecule", "Epitope - Name"]
]

subset = subset.rename(columns={
        "Epitope - Source Molecule": "Source",
        "Epitope - Name": "Epitope",
        }
)

subset = subset.drop_duplicates(subset=["Source", "Epitope"])

# Limit epitopes longer than 21, because TLimmuno2 only accepts epitopes up to 21 amino acids
epitope_len = subset["Epitope"].apply(lambda e: len(e.split(" ")[0]))
subset = subset[epitope_len <= 21]


subset.to_excel(args.output, index=False)
