import argparse
import pandas as pd

def preprocessFiles(geno_file, express_file, cov_file, geno_out, express_out, cov_out):
    #genotype = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/GTExNomalizedSNPGenotypes_chr1_first10_samplename.table', sep='\t')
    #expression = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/Clean_expression_first10.tsv', sep='\t')

    genotype = pd.read_csv(geno_file, sep='\t')
    expression = pd.read_csv(express_file, sep='\t')
    #change column naming
    genotype['start'] = genotype['start'].astype(str)
    genotype['chrom_start'] = genotype[['chrom', 'start']].agg('_'.join, axis=1)
    genotype.columns = genotype.columns.str.replace('-', '.')

    #find intersection between genotype and expression samples
    geno_col = set(genotype.columns)
    exp_col = set(expression.columns)
    intersect = geno_col.intersection(exp_col)
    chrom_start = list(genotype['chrom_start'])
    geno_inter = genotype[intersect]
    geno_inter.insert(0, 'chrom_start', chrom_start)
    exp_inter = expression[intersect]
    #geno_inter.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/GTExNomalizedSNPGenotypes_chr1_first10_samplename_intersect.table', sep='\t', index=False)
    geno_inter.to_csv(geno_out, sep='\t', index=False)
    #exp_inter.to_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/Clean_expression_first10_intersect.tsv', sep='\t')
    exp_inter.to_csv(express_out, sep='\t')
    #covariates = pd.read_csv('/storage/cynthiawu/trans_eQTL/gtex650_transpose_index_relabel_tab_samplename.pca', sep='\t')
    covariates = pd.read_csv(cov_file, sep='\t')
    cov_inter = covariates[intersect]
    #cov_inter.to_csv('/storage/cynthiawu/trans_eQTL/gtex650_transpose_index_relabel_tab_samplename_inter.pca', sep='\t')
    cov_inter.to_csv(cov_out, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genotype", required=True, help="Input genotype file")
    parser.add_argument("-e", "--expression", required=True, help="Input expression file")
    parser.add_argument("-c", "--covariates", required=True, help="Input covariates file")
    parser.add_argument("-o", "--geno_out", required=True, help="Output genotype file with intersecting samples")
    parser.add_argument("-p", "--express_out", required=True, help="Output expression file with intersecting samples")
    parser.add_argument("-q", "--cov_out", required=True, help="Output covariates file with intersecting samples")
    params = parser.parse_args()
    preprocessFiles(params.genotype, params.expression, params.covariates, params.geno_out, params.express_out, params.cov_out)


if __name__ == "__main__":
    main()
