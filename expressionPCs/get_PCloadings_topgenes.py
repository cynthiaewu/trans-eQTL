import pandas as pd
import argparse


def get_topgenes(input, output):
    loadings = pd.read_csv(input, sep='\t', index_col=0)
    normalized_loadings = (loadings-loadings.mean())/loadings.std()
    pcs = []
    positive = []
    negative = []
    for column in normalized_loadings:
        pcs.append(column)
        positive.append(list(normalized_loadings.loc[normalized_loadings[column] > 3].index))
        negative.append(list(normalized_loadings.loc[normalized_loadings[column] < -3].index))

    f = open(output, "w")
    f.write(f'PC\tPositive_corr\tNegative_corr\n')
    for i in range(len(pcs)):
        f.write(f'{pcs[i]}\t{positive[i]}\t{negative[i]}\n')       
    f.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input", required=True, help="Input expression file")
    parser.add_argument("-o", "--output", required=True, help="Output expressionPCs file")
    params = parser.parse_args()
    get_topgenes(params.input, params.output)


if __name__ == "__main__":
    main()
