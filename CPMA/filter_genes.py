import pandas as pd
import numpy as np
import argparse


def filter_genes(input, output):
    gene_annot = pd.read_csv('/storage/cynthiawu/trans_eQTL/gene_annotations_gencode.v19.genes.v7.csv', sep='\t', header=None)
    #expression = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr2/Clean_expression_chr2_inter.tsv', sep='\t')
    expression = pd.read_csv(input, sep='\t')
    
    #Get protein coding genes.
    p_coding = []
    for index, row in gene_annot.iterrows():
        if row[1] == 'protein_coding':
            p_coding.append(row[0])  
    
    expression_pcoding = expression[expression['Unnamed: 0'].isin(p_coding)] 
    expression_pcoding_mean = expression_pcoding.loc[(expression_pcoding != 0).any(axis=1)]
   
    #get genes with rpkm > 0.1 in at least 10 samples. 
    data = expression_pcoding_mean.set_index('Unnamed: 0')
    indices = np.where(np.sum(data > 0.1, axis=1) >= 10)[0]
    filtered_data = data.iloc[indices,:]
    filtered_data.to_csv(output, header=True, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input expression file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput filtered expression file, only protein coding genes, rpkm > 0.1 in 10 samples")
    params = parser.parse_args()
    filter_genes(params.input, params.output)


if __name__ == "__main__":
    main()
