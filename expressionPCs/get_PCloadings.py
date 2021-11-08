import pandas as pd
import numpy as np
import sklearn.decomposition
import argparse


def annotate_ids(gene_ids):
    path = "/storage/cynthiawu/trans_eQTL/GTex.v8/gencode.v26.GRCh38.genes.annotations"
    gene_dict = {}
    with open(path) as f:
        for line in f:
            values = line.strip().split('\t')
            gene_dict[values[0]] = values[1]

    annotate = []
    for value in gene_ids:
        if value in gene_dict:
            annotate.append(gene_dict[value])
        else:
            annotate.append(value)
    return annotate


def get_PCloadings(input, output):
    expression = pd.read_csv(input, sep='\t', index_col=0)
    samples = list(expression.columns)
    expression_trans = expression.values.transpose()
    pca = sklearn.decomposition.PCA()
    pca.fit(expression_trans)
    expression_trans_pca = pca.transform(expression_trans)
    
    #annotate = annotate_ids(list(expression.index))
    col = ['PC' + str(i) for i in range(len(pca.explained_variance_))]
    #loadings = pd.DataFrame(pca.components_.T* np.sqrt(pca.explained_variance_), columns=col, index=annotate)
    loadings = pd.DataFrame(pca.components_.T* np.sqrt(pca.explained_variance_), columns=col, index=expression.index)
    loadings.to_csv(output, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input", required=True, help="Input expression file")
    parser.add_argument("-o", "--output", required=True, help="Output expressionPCs file")
    params = parser.parse_args()
    get_PCloadings(params.input, params.output)


if __name__ == "__main__":
    main()
