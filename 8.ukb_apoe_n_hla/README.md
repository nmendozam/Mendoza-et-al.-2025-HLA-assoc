From UK Biobank whole genome sequencing data the genotyping of
ApoE and HLA-DRB1\*15:01 was extracted
```bash
bcftools view -r chr19:44908684,chr19:44908822 ukb_wgs/chr19.ukb23374_v1_cfilter_restored_mac4.phased.bcf > apoe.vcf
bgzip -c apoe.vcf > apoe.vcf.gz
bcftools index apoe.vcf.gz
bcftools view -r chr6:32445274 ukb_wgs/chr6.ukb23374_v1_cfilter_restored_mac4.phased.bcf > hla.vcf
bgzip -c hla.vcf > hla.vcf.gz
bcftools index hla.vcf.gz
bcftools concat *gz > merged.vcf
```

# Results of interaction analysis
After running the analysis of interaction between APOE4 and HLA-DRB1\*15:01.
This interaction test is described in PyHLA but based on other to previous papers
```
A is APOE2 and B is HLA-DRB1*15:01
     p-val     OR  Test
 0.0039588   1.11  [1] A associated?
 0.0001979   0.88  [2] B associated?
 0.1326132   1.12  [3] A associated in B-positives?
 0.0138869   1.10  [4] A associated in B-negatives?
 0.1737417   0.89  [5] B associated in A-positives?
 0.0005926   0.88  [6] B associated in A-negatives?
 6.451e-06   1.24  [7] Difference between A and B association?
 0.9721910   0.99  [8] Combined A-B association?
 0.9663468   0.99  [9] Linkage disequilibrium in cases
 0.0763062   0.98  [10] Linkage disequilibrium in controls

A is APOE4 and B is HLA-DRB1*15:01
     p-val     OR  Test
 4.173e-13   0.79  [1] A associated?
 0.0001979   0.88  [2] B associated?
 5.466e-05   0.77  [3] A associated in B-positives?
 1.788e-09   0.80  [4] A associated in B-negatives?
 0.0210421   0.86  [5] B associated in A-positives?
 0.0037236   0.90  [6] B associated in A-negatives?
 0.0159020   0.89  [7] Difference between A and B association?
 1.540e-10   0.69  [8] Combined A-B association?
 0.7967601   0.98  [9] Linkage disequilibrium in cases
 0.0011641   1.02  [10] Linkage disequilibrium in controls
```


# For the epistasis analysis with plink

```bash
cut -d ' ' -f 1,2,7 02_parental_age_phenotypes/father_age_death.categorical.tsv.csv > father_age_death.categorical.txt
python interaction_assoc.py
plink --vcf merged.vcf --make-pheno father_age_death.categorical.txt '*' --make-bed --out merged

```