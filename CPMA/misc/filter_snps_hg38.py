import argparse
import numpy as np
import sqlite3
import pandas as pd

def getRegions(input, snps_file, output, chr):
    #genotype_all = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/GTExNormalizedSNPGenotypes_chr1_samplename_inter.table', sep='\t')
    genotype_all = pd.read_csv(input, sep='\t')
    coding_regions = pd.read_csv('/storage/polo/GTEx_v8_processed_genotype/coding_bed/hg38_coding_exon.bed', skiprows=1, header=None, sep='\t')
    coding_regions = coding_regions.loc[coding_regions[0].isin([f'chr{chr}'])]
    #coding_regions = coding_regions.loc[coding_regions[0].isin([f'chr{chr}', 'chr1_gl000191_random', 'chr1_gl000192_random'])]
    coding_regions.columns = ['chr', 'start', 'end', 'region', '0', 'strand']
    genotype_all = genotype_all.join(genotype_all['chrom_start'].str.split('_', 1, expand=True).rename(columns={0:'chr', 1:'pos'}))
    #print(genotype_all.columns)
    genotype_all['pos'] = genotype_all['pos'].astype(int)
    print('tables read')

    #Make the db in memory
    conn = sqlite3.connect(':memory:')
    #write the tables
    coding_regions.to_sql('coding', conn, index=False)
    genotype_all.to_sql('genotype', conn, index=False)
    print('tables saved in sql')

    qry = '''
        select 
        *
        from
            genotype join coding on
            pos between start and end
       '''
    df = pd.read_sql_query(qry, conn)

    df = df.drop(columns=['chr', 'pos', 'start', 'end', 'region', '0', 'strand'])
    df = df.drop_duplicates()
    print('tables finished')

    # /storage/cynthiawu/trans_eQTL/GTex_filteredsnps_pos.txt
    with open(snps_file) as f:
        snp_pos = f.read().splitlines() 

    df = df[df['chrom_start'].isin(snp_pos)]
    #df.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/GTExNormalizedSNPGenotypes_chr1_samplename_inter_coding.table', sep='\t')
    df.to_csv(output, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input genotype file")
    parser.add_argument("-s", "--snps_file", required=True, help="Input file with the formatted snp position that passed filters")
    parser.add_argument("-o", "--output", required=True, help="Output file with genotype for coding regions")
    parser.add_argument("-c", "--chr", required=True, help="Chromosome Number")
    params = parser.parse_args()
    getRegions(params.input, params.snps_file, params.output, params.chr)


if __name__ == "__main__":
    main()
