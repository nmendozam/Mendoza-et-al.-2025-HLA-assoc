"""
Read results from NetMHCIIpan and get the epitopes
with a high rank only for HLA-DRB1*15:01. This
will later be used to predict immunogenicity.
"""

import argparse

import pandas as pd

# %%

parser = argparse.ArgumentParser(
    description="Filter NetMHCIIpan results for HLA-DRB1*15:01"
)
parser.add_argument(
    "--input",
    type=str,
    default="3590022_NetMHCIIpan.xlsx",
    help="Input file with NetMHCIIpan results",
)
parser.add_argument(
    "--output",
    type=str,
    default="filtered_epitopes.csv",
    help="Output file with filtered epitopes",
)

args = parser.parse_args()

df = pd.read_excel(args.input, engine="openpyxl", header=[0, 1])

df.head()
# %%
high_rank = df[df[("DRB1_1501", "Rank")] <= 0.05]
# %%
# filter ranks higher than DRB1_1501's rank
alleles = high_rank.columns.get_level_values(0).unique()
alleles = [
    allele for allele in alleles if allele != "DRB1_1501" and "Unnamed" not in allele
]

filtered = high_rank.copy()
for allele in alleles:
    filtered = filtered[filtered[(allele, "Rank")] >= filtered[("DRB1_1501", "Rank")]]



# %% Make to columns one with Core and the other with DRB1_1501
final = filtered["DRB1_1501"].sort_values(by="Score", ascending=False)

final["allele"] = "DRB1_1501"
final = final[["Core", "allele"]]

# Remove core duplicates
final = final.drop_duplicates(subset=["Core"])

print(final)

final.to_csv(args.output, index=False, header=False, sep=",")


# %%
