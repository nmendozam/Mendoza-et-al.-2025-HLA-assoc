import pandas as pd

# Read phenotype file
pheno = pd.read_csv('Parental_age_phenotypes.csv')
pheno.drop('Ethnicity', axis=1, inplace=True)
pheno.drop('parental LLI', axis=1, inplace=True)

# Read ethnicity file
ethnicity = pd.read_csv('Self_Etnicity_datafile.tsv', sep='\t')
ethnicity.columns = ['eid', 'ethnicity']

# Merge both files
merged = pheno.merge(ethnicity, on='eid')
print(merged.head())

merged.to_csv('family_phenotypes_with_ethnicity.tsv', sep='\t', index=False)

def LL_parent_phenotype(person, parent_kind: str):
    """Find long-live parents and cases. Parent kind is either fathers or mothers"""
    # If not white or british
    if not (person.ethnicity == 1001 or person.ethnicity == 1):
        return 0 # Exclude 

    # Look for long-lived parents
    if person[f"{parent_kind}_age_at_death"] >= 95 or person[f"{parent_kind}_age"] >= 95:
        return 2 # Case

    # If age at death is unknown exclude 
    if person[f"{parent_kind}_age_at_death"] == 0:
        return 0 # Exclude

    # Look for long-lived parents
    if person[f"{parent_kind}_age_at_death"] >= 90 or person[f"{parent_kind}_age"] >= 90:
        return 0 # Exclude

    # Exclude young parental deaths
    if person[f"{parent_kind}_age_at_death"] <= 40 or person[f"{parent_kind}_age"] <= 40:
        return 0 # Exclude

    return 1 # Everyone else is control

# Calculate long-lived father categorical phenotype
father_LLI = merged.loc[:, ["eid", "fid" , "na", "na.1", "sex"]]
father_LLI["Pheno"] = merged.apply(lambda row: LL_parent_phenotype(row, "fathers"), axis=1)
father_LLI.to_csv("father_lli.with_ethnicity.csv", index=False, header=False, sep=" ")

# Calculate long-lived mother categorical phenotype
mother_LLI = merged.loc[:, ["eid", "fid" , "na", "na.1", "sex"]]
mother_LLI["Pheno"] = merged.apply(lambda row: LL_parent_phenotype(row, "mothers"), axis=1)
mother_LLI.to_csv("mother_lli.with_ethnicity.csv", index=False, header=False, sep=" ")

