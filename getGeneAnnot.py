import allel
import pandas as pd

#path = '/storage/resources/datasets/gtex/53844/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WGS_SNP_CNV.genotype-calls-vcf.c1/GTEx_Analysis_20150112_WholeGenomeSeq_VarSitesAnnot.vcf.gz'
path = '/storage/cynthiawu/trans_eQTL/Annotations/GTEx_Analysis_20150112_WholeGenomeSeq_VarSitesAnnot_filtered.recode.vcf'
df = allel.vcf_to_dataframe(path, fields=['variants/CHROM', 'variants/POS', 'variants/CSQ'])
print('VCF file read')
df['CSQ'] = df['CSQ'].str.split('|').str[3]
print('Gene annotation included')
df.to_csv('/storage/cynthiawu/trans_eQTL/Annotations/snp-geneAnnot.csv', index=False, sep = '\t')
