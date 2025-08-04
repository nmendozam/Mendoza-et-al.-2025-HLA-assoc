#%% import libraries
import pandas as pd
import scipy

# Input parameters
vcf_file = "merged.vcf"
snp_ids = [
    # "chr19:44908684:SG", # This is for APOE-ε4
    "chr19:44908822:SG", # This is for APOE-ε2
    "chr6:32445274:SG",
]
allele_code = ["1", "0"] # 1 is the minor allele, 0 is the major allele

# Functions
def read_vcf(file_name):
    num_header = 0
    with open(file_name) as file:
        for line in file.readlines():
            if line.startswith("##"):
                num_header += 1
            else:
                break
    vcf = pd.read_csv(file_name, sep="\t", skiprows=num_header)
    return vcf.rename({"#CHROM": "CHROM"}, axis=1)

def count_allele_interact(allele):
    pp = int(allele.all(axis=1).sum())
    pn = int((allele.iloc[:,0] & ~allele.iloc[:,1]).sum())
    np = int((~allele.iloc[:,0] & allele.iloc[:,1]).sum())
    nn = int((~allele.iloc[:,0] & ~allele.iloc[:,1]).sum())
    return [pp, pn, np, nn]


def simpleTest(data, test = 'fisher'):
    '''
    This function performs a simple association test
    using either Fisher's exact test or Chi-squared test.
    '''
    pvalue = 'NA'
    if test == 'fisher':
        try:
            OR, pvalue = scipy.stats.fisher_exact(data)
            # Calculate confidence interval

        except:
            pvalue = 'NA'
    elif test == 'chisq':
        try:
            chi2, pvalue, dof, expected = scipy.stats.chi2_contingency(data)
        except:
            pvalue = 'NA'
    return pvalue
def tenTests(x, y, test = 'fisher'):
    '''
    Ten interaction tests for two alleles A and B.
    The tests are based on the 2x2 contingency table of counts of the two alleles.
    The implementation was taken from PyHLA
    '''
    n1 = n2 = n3 = n4 = 0
    # test 1
    n1 = x[0] + x[1]
    n2 = x[2] + x[3]
    n3 = y[0] + y[1]
    n4 = y[2] + y[3]
    data = [[n1, n2], [n3, n4]]
    OR1 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p1 = simpleTest(data, test)
    # test2
    n1 = x[0] + x[2]
    n2 = x[1] + x[3]
    n3 = y[0] + y[2]
    n4 = y[1] + y[3]
    data = [[n1, n2], [n3, n4]]
    OR2 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p2 = simpleTest(data, test)
    # test 3
    n1 = x[0]
    n2 = x[2]
    n3 = y[0]
    n4 = y[2]
    data = [[n1, n2], [n3, n4]]
    OR3 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p3 = simpleTest(data, test)
    # test 4
    n1 = x[1]
    n2 = x[3]
    n3 = y[1]
    n4 = y[3]
    data = [[n1, n2], [n3, n4]]
    OR4 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p4 = simpleTest(data, test)
    # test5
    n1 = x[0]
    n2 = x[1]
    n3 = y[0]
    n4 = y[1]
    data = [[n1, n2], [n3, n4]]
    OR5 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p5 = simpleTest(data, test)
    # test6
    n1 = x[2]
    n2 = x[3]
    n3 = y[2]
    n4 = y[3]
    data = [[n1, n2], [n3, n4]]
    OR6 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p6 = simpleTest(data, test)
    # test7
    n1 = x[1]
    n2 = x[2]
    n3 = y[1]
    n4 = y[2]
    data = [[n1, n2], [n3, n4]]
    OR7 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p7 = simpleTest(data, test)
    # test8
    n1 = x[0]
    n2 = x[3]
    n3 = y[0]
    n4 = y[3]
    data = [[n1, n2], [n3, n4]]
    OR8 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p8 = simpleTest(data, test)
    # calculate confidence intervals
    import math
    se = math.sqrt(1.0/n1 + 1.0/n2 + 1.0/n3 + 1.0/n4)
    l95 = math.exp(math.log(OR8) - 1.96 * se)
    u95 = math.exp(math.log(OR8) + 1.96 * se)
    print("Test 8: P=%4.2e OR=%0.2f 95%% CI: %0.2f-%0.2f" % (p8, OR8, l95, u95))

    # test9
    n1 = x[0]
    n2 = x[1]
    n3 = x[2]
    n4 = x[3]
    data = [[n1, n2], [n3, n4]]
    OR9 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p9 = simpleTest(data, test)
    # test10
    n1 = y[0]
    n2 = y[1]
    n3 = y[2]
    n4 = y[3]
    data = [[n1, n2], [n3, n4]]
    OR10 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p10 = simpleTest(data, test)
    return [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,OR1,OR2,OR3,OR4,OR5,OR6,OR7,OR8,OR9,OR10]

#%% Read file

vcf = read_vcf(vcf_file)
vcf.drop(
    ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
    axis=1,
    inplace=True,
)

# %% Read vcf file and phenotype file
vcf = vcf.reset_index().set_index("ID")
vcf = vcf.transpose()
vcf = vcf.drop("index")

phenotype_file = "02_parental_age_phenotypes/father_age.categorical.tsv.csv"
column_names = [
    "IID", "FID", "X", "XX", "SEX", "CASE"
] 
phe = pd.read_csv(phenotype_file, sep=" ", names=column_names)

ctl_ids = phe.IID[phe.CASE == 1]
case_ids = phe.IID[phe.CASE == 2]
case_ids = list(set(case_ids.astype(str)).intersection(set(vcf.index)))
ctl_ids = list(set(ctl_ids.astype(str)).intersection(set(vcf.index)))

# %% Perform interaction association test
def get_allele_carriers(vcf, snp_ids, allele_codes):
    carriers = pd.DataFrame(index=vcf.index, columns=snp_ids)
    for snp_id, allele_code in zip(snp_ids, allele_codes):
        carriers[snp_id] = vcf[snp_id].str.contains(allele_code)

    return carriers

carriers = get_allele_carriers(vcf, snp_ids, allele_code)


case_carriers = carriers.loc[list(case_ids), snp_ids]
ctl_carriers = carriers.loc[list(ctl_ids), snp_ids]

x = count_allele_interact(case_carriers) # cases
y = count_allele_interact(ctl_carriers) # controls

results = tenTests(x, y, 'fisher')

# Both test 3 and test 4 are significant: A is associated with the disease independently of B.
if results[2] <= 0.05 and results[3] <= 0.05:
    print("A is associated with the disease independently of B.")
    print(results[2], results[3])
# Both test 5 and test 6 are significant: B is associated with the disease independently of A.
if results[4] <= 0.05 and results[5] <= 0.05:
    print("B is associated with the disease independently of A.")
    print(results[4], results[5])
# Both test 3 and test 5 are significant
if results[2] <= 0.05 and results[4] <= 0.05:
    print("A and B show interaction.")
    print(results[2], results[4])
# Test 7 is significant
if results[6] <= 0.05:
    print("Difference between A and B is associated with the disease.")
    # print(results[6])
    print(f"\tp-value: {results[6]} OR: {results[16]}")
# Test 8 is significant
if results[7] <= 0.05:
    print("A and B have combined action.")
    print("P=%4.2e" % results[7], "OR=%0.2f" % results[17])
# Test 9 is significant
if results[8] <= 0.05:
    print("A and B are in LD in cases.")
    print(results[8])
# Test 10 is significant
if results[9] <= 0.05:
    print("A and B are in LD in controls.")
    print(results[9])

print(results)
# %%
print("Perform regular test for each SNP")
# Perform regular test for each SNP
for snp_id in snp_ids:
    cases = case_carriers[snp_id].astype(bool)
    controls = ctl_carriers[snp_id].astype(bool)
    n1 = cases.sum()
    n2 = (~cases).sum()
    n3 = controls.sum()
    n4 = (~controls).sum()
    data = [[n1, n2], [n3, n4]]
    OR = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p = simpleTest(data, 'fisher')

    # Calculate confidence interval
    import math
    se = math.sqrt(1.0/n1  + 1.0/n2 +  1.0/n3 + 1.0/n4)
    l95 = math.exp(math.log(OR) - 1.96 * se)
    u95 = math.exp(math.log(OR) + 1.96 * se)
    

    print(f"\tSNP: {snp_id}")
    print(f"\t\tOdds ratio: {OR} p-value: {p}")
    print(f"\t\t95% CI: {l95} - {u95}")
    print("\n")

# %%
