#%% imports
import os
import pandas as pd

gene_names = [
    "A", "A_1", "B", "B_1", "C", "C_1", "DPA1", "DPA1_1", "DPB1", "DPB1_1", "DQA1", "DQA1_1", "DQB1", "DQB1_1", "DRB1", "DRB1_1"
]

# %% Get the data
random_seq_file_name = 'Accuracy measurement/random_seq.csv'
random_chip_file_name = 'Accuracy measurement/random_chip.csv'

random_seq = pd.read_csv(random_seq_file_name, index_col=0)
random_chip = pd.read_csv(random_chip_file_name, index_col=0)


# %% Liquify keep source index on melt
def liquify(df, value_name="allele"):
    df = df.melt(value_vars=gene_names, var_name="gene", value_name=value_name, ignore_index=False).reset_index()
    return df.rename(columns={0: "Patient-ID"})

melt_seq = liquify(random_seq, "allele_seq")
melt_chip = liquify(random_chip, "allele_chip")


# %% Calculate accuracy
def normalize_resolution(allele_df, allele_column_name, allele_resolution=1):
    for column in allele_column_name:
        allele_df[column] = allele_df[column].map(lambda x: x if x != x else ":".join(x.split(":")[0:allele_resolution]))
    return allele_df

# %% Generate concordance table for DRB1
drb1 = pd.merge(
    random_seq[["DRB1", "DRB1_1"]],
    random_chip[["DRB1", "DRB1_1"]],
    left_index=True,
    right_index=True,
    how="inner",
    suffixes=("", "_chip"))

drb1 = drb1.rename(columns={"DRB1": "seq", "DRB1_1": "seq_1", "DRB1_chip": "chip", "DRB1_1_chip": "chip_1"})
drb1 = normalize_resolution(drb1, ["seq", "seq_1", "chip", "chip_1"], allele_resolution=3)
drb1.fillna("nan", inplace=True)

# swap alleles that require it
need_swap = ((drb1.seq != drb1.chip) & (drb1.seq_1 != drb1.chip_1)
    &
    ((drb1.seq == drb1.chip_1) | (drb1.seq_1 == drb1.chip)))
temp = drb1.seq.copy()
drb1.loc[need_swap, "seq"] = drb1.seq_1[need_swap]
drb1.loc[need_swap, "seq_1"] = temp[need_swap]

# Compute concordance
table = pd.crosstab(drb1["chip"], drb1["seq"], colnames=["HLA-seq"], rownames=["Immunochip"])
table2 = pd.crosstab(drb1["chip_1"], drb1["seq_1"], colnames=["HLA-seq"], rownames=["Immunochip"])
table = table.add(table2, fill_value=0)

table.to_csv("Accuracy measurement/DRB1_concordance.csv")

# table["DRB1*15:01:01"]
# %% Calculate statistics

drb1 = pd.merge(
    random_seq[["DRB1", "DRB1_1"]],
    random_chip[["DRB1", "DRB1_1"]],
    left_index=True,
    right_index=True,
    how="inner",
    suffixes=("", "_chip"))

drb1 = drb1.rename(columns={"DRB1": "seq", "DRB1_1": "seq_1", "DRB1_chip": "chip", "DRB1_1_chip": "chip_1"})
drb1 = normalize_resolution(drb1, ["seq", "seq_1", "chip", "chip_1"], allele_resolution=3)
drb1.fillna("nan", inplace=True)

# swap alleles that require it
need_swap = ((drb1.seq != drb1.chip) & (drb1.seq_1 != drb1.chip_1)
    &
    ((drb1.seq == drb1.chip_1) | (drb1.seq_1 == drb1.chip)))
temp = drb1.seq.copy()
drb1.loc[need_swap, "seq"] = drb1.seq_1[need_swap]
drb1.loc[need_swap, "seq_1"] = temp[need_swap]


allele = "DRB1\\*15:01:01"
c = drb1.chip.str.contains(allele)
c1 = drb1.chip_1.str.contains(allele)
s = drb1.seq.str.contains(allele)
s1 = drb1.seq_1.str.contains(allele)

TP = (c & s).sum() + (c1 & s1).sum()
TN = (~c & ~s).sum() + (~c1 & ~s1).sum()
FP = (c & ~s).sum() + (c1 & ~s1).sum()
FN = (~c & s).sum() + (~c1 & s1).sum()

accuracy = (TP + TN) / (TP + TN + FP + FN)
precision = TP / (TP + FP)
recall = TP / (TP + FN)
specificity = TN / (TN + FP)
f1 = 2 * (precision * recall) / (precision + recall)

print("TP: %d, FP: %d, TN: %d, FN: %d" % (TP, FP, TN, FN))

print("precission: %.02f%%" % (100 * precision))
print("accuracy: %.02f%%" % (100 * accuracy))
print("recall: %.02f%%" % (100 * recall))
print("specificity: %.02f%%" % (100 * specificity))
