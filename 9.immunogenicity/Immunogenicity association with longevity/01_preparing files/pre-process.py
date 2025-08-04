import itertools

import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='prepare peptide-allele pairs for TLimmuno')
parser.add_argument(
        'input_epitopes',
        type=str,
        help='Selection of epitopes. Excel file with Epitope name and sequence',
)
parser.add_argument(
        'input_alleles',
        type=str,
        help='Selection of alleles. CSV file with the alleles you want to be tested. Only alleles from the DRB locus will be included in the final file.',
)
parser.add_argument(
        'output',
        type=str,
        help='Path to the output csv file',
        default='immunogenicity_input.csv',
)
args = parser.parse_args()


epitopes = pd.read_excel(args.input_epitopes)
alleles = pd.read_csv(args.input_alleles)

DRB1 = alleles[alleles.Allele.str.contains('DRB1')].Allele.apply(lambda x: x.replace('*', '_').replace(':', '')).to_list()
epitopes = epitopes.Epitope.apply(lambda x: x.split(' ')[0]).to_list()

epitopes = list(set(epitopes)) # remove duplicates

df = pd.DataFrame(list(itertools.product(epitopes, DRB1)),columns=['Epitope','allele'])
df.to_csv(args.output,index=False,header=False)

df.head()



# python Python/TLimmuno2.py --mode file --intdir Immunogenicity scores/immunogenicity_input.csv --outdir .
