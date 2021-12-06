import argparse
import pandas as pd

def intersectFiles(geno_file, express_file, cov_file, geno_out, express_out, cov_out):
    #genotype = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/GTExNomalizedSNPGenotypes_chr1_first10_samplename.table', sep='\t')
    #expression = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/Clean_expression_first10.tsv', sep='\t')

    genotype = pd.read_csv(geno_file, sep='\t', index_col=0)
    #expression = pd.read_csv(express_file, sep='\t', index_col=0)
    expression = pd.read_csv(express_file, index_col=0)
    
    covariates = pd.read_csv(cov_file, sep='\t', index_col=0)
    # Should include only European samples if European covariates file is chosen
    
    #change column naming
    expression.columns = expression.columns.str.replace('.', '-')

    cov_col = set(covariates.columns)
    #find intersection between genotype and expression samples
    geno_col = set(genotype.columns)
    exp_col = set(expression.columns)
    intersect = (geno_col.intersection(exp_col)).intersection(cov_col)

    # get intersecting columns and write to file
    geno_inter = genotype[intersect]
    geno_inter.to_csv(geno_out, sep='\t')

    exp_inter = expression[intersect]
    exp_inter.to_csv(express_out, sep='\t')
    
    cov_inter = covariates[intersect]
    cov_inter.to_csv(cov_out, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genotype", required=True, help="Input genotype file")
    parser.add_argument("-e", "--expression", required=True, help="Input expression file")
    parser.add_argument("-c", "--covariates", required=True, help="Input covariates file")
    parser.add_argument("-i", "--geno_out", required=True, help="Output genotype file of intersecting samples")
    parser.add_argument("-p", "--express_out", required=True, help="Output expression file with intersecting samples")
    parser.add_argument("-q", "--cov_out", required=True, help="Output covariates file with intersecting samples")
    params = parser.parse_args()
    intersectFiles(params.genotype, params.expression, params.covariates, params.geno_out, params.express_out, params.cov_out)


if __name__ == "__main__":
    main()
