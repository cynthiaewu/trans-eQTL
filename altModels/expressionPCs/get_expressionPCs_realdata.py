import pandas as pd
import sklearn.decomposition
import argparse


def get_expressionPCs(input, output):
    #expression = pd.read_csv('/gymreklab-tscc/cynthiawu/Test_nullsnps/simulate_eqtls_only/FastMultivariate/Single_eqtl/SampleSize100/SingleParameter/numTarget_1000/Beta_1/Simulation_0/expression.csv', sep='\t')
    expression = pd.read_csv(input, sep='\t', index_col=0)
    #expression.set_index('Unnamed: 0', inplace=True)
    samples = list(expression.columns)
    expression_trans = expression.values.transpose()
    pca = sklearn.decomposition.PCA()
    pca.fit(expression_trans)
    expression_trans_pca = pca.transform(expression_trans)
    PCs = ['PC' + str(i) for i in range(len(expression_trans_pca))]
    PC_trans_df = pd.DataFrame(expression_trans_pca.T, columns=samples, index=PCs)
    PC_trans_df.to_csv(output, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input", required=True, help="Input expression file")
    parser.add_argument("-o", "--output", required=True, help="Output expressionPCs file")
    params = parser.parse_args()
    get_expressionPCs(params.input, params.output)


if __name__ == "__main__":
    main()
