import math
import argparse

import pandas as pd


def get_column_names(df):
    """
    Gets the gene names from the first row and returns
    a list with the names for each column. This method
    assumes that the df is in pyHLA format
    """
    genes = df.iloc[1,2:].apply(lambda x: x.split("*")[0]).to_list()
    genes_index = list()
    for gene in genes:
        if gene in genes_index:
            genes_index.append(gene + "_1")
        else:
            genes_index.append(gene)
    return ["ID", "LLI"] + genes_index


def norm_resolution(df, allele_resolution=2):
    """
    Normalizes the allele resolutions to the indicated level
    """
    df = df.fillna("")
    df.iloc[:,2:] = df.iloc[:,2:].map(lambda x: ":".join(x.split(":")[0:allele_resolution]))
    return df


parser = argparse.ArgumentParser(description='Prepare immunogenicity scores for statistical analysis')
parser.add_argument(
        'input',
        type=str,
        help='CSV with the resulting scores from TLimmuno2',
        default='result.csv',
)
parser.add_argument(
        'genotype_file',
        type=str,
        help='input file with the genotype data of the individuals in pyHLA format',
        default = 'Ageing.2field.pyhla',
)
parser.add_argument(
        'output_basename',
        type=str,
        help='basename of the output files of this script',
        default='DRB1'
)

args = parser.parse_args()

img = pd.read_csv(args.input)

img['BindingAffinity(nM)'] = img['prediction']#.apply(lambda x: 50000 ** (1 - x))

# pivot the data
group = img.groupby(['pep', 'HLA']).last().reset_index()
binding_affinity = group.pivot(index=['HLA'], columns='pep', values='BindingAffinity(nM)')
binding_affinity.index = binding_affinity.index.str.replace(r'(.*)_(\d{2})(\d{2})', r'\1*\2:\3', regex=True)


genotypes = pd.read_csv(args.genotype_file, sep='\t', header=None)
genotypes.columns = get_column_names(genotypes)
genotypes = norm_resolution(genotypes, allele_resolution=2)
genotypes.head()

genotypes = genotypes[['ID', 'LLI', 'DRB1', 'DRB1_1']]

binding_affinity.head()




# Merge the genotypes with the binding scores
genotypes['HLA'] = genotypes['DRB1']
result = pd.merge(genotypes, binding_affinity, left_on='HLA', right_index=True)
genotypes['HLA'] = genotypes['DRB1_1']
result_1 = pd.merge(genotypes, binding_affinity, left_on='HLA', right_index=True)
result.head()
result_1.head() 

# mean of values from result and result_1 and drop the HLA column
result = result.drop(columns=['HLA'])
result_1 = result_1.drop(columns=['HLA'])
result = result.set_index('ID')
result_1 = result_1.set_index('ID')

# save the result
result.to_csv(args.output_basename + '.csv')
result_1.to_csv(args.output_basename + '_1.csv')
