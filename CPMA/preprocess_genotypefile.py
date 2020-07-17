import argparse
import pandas as pd

def preprocessFiles(geno_file, inter_file, geno_out):
    #genotype = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/GTExNomalizedSNPGenotypes_chr1_first10_samplename.table', sep='\t')
    #expression = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/Clean_expression_first10.tsv', sep='\t')

    genotype = pd.read_csv(geno_file, sep='\t')
    with open(inter_file) as f:
        intersect = f.read().splitlines()
    #change column naming
    genotype['start'] = genotype['start'].astype(str)
    genotype['chrom_start'] = genotype[['chrom', 'start']].agg('_'.join, axis=1)
    genotype.columns = genotype.columns.str.replace('-', '.')

    #find intersection between genotype and expression samples
    geno_col = set(genotype.columns)
    chrom_start = list(genotype['chrom_start'])
    geno_inter = genotype[intersect]
    geno_inter.insert(0, 'chrom_start', chrom_start)
    #geno_inter.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/GTExNomalizedSNPGenotypes_chr1_first10_samplename_intersect.table', sep='\t', index=False)
    geno_inter.to_csv(geno_out, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genotype", required=True, help="Input genotype file")
    parser.add_argument("-i", "--inter_samples", required=True, help="Input file with intersecting samples")
    parser.add_argument("-o", "--geno_out", required=True, help="Output genotype file with intersecting samples")
    params = parser.parse_args()
    preprocessFiles(params.genotype, params.inter_samples, params.geno_out)


if __name__ == "__main__":
    main()
