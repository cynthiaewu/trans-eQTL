import pandas as pd
import argparse


def filter_genes(input, output):
    gene_annot = pd.read_csv('/storage/cynthiawu/trans_eQTL/gene_annotations_gencode.v19.genes.v7.csv', sep='\t', header=None)
    #expression = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr2/Clean_expression_chr2_inter.tsv', sep='\t')
    expression = pd.read_csv(input, sep='\t')
    
    p_coding = []
    for index, row in gene_annot.iterrows():
        if row[1] == 'protein_coding':
            p_coding.append(row[0])  
    
    expression_pcoding = expression[expression['Unnamed: 0'].isin(p_coding)] 
    expression_pcoding_mean = expression_pcoding.loc[(expression_pcoding != 0).any(axis=1)]
    expression_pcoding_mean.to_csv(output, index=False, header=True, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input expression file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput filtered expression file, only protein coding genes, median rpkm > 0")
    params = parser.parse_args()
    filter_genes(params.input, params.output)


if __name__ == "__main__":
    main()
