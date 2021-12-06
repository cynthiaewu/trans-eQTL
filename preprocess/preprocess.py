import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats


def preprocess(tissue, output):
    EXPRFILE = "/storage/polo/GTEx_v8_European_only/%s_Euro/cleaned_tpm.csv"%tissue
    expr = pd.read_csv(EXPRFILE, index_col="Name")
    
    # Filter outlier genes
    means = expr.mean(axis=1)
    expr = expr[means<6000]

    # Require genes to have TPM >0.1 in >10 samples 
    expression_threshold=0.1
    min_samples = 10
    mask = ((np.sum(expr>expression_threshold,axis=1)>=min_samples)).values
    # apply normalization
    # Sample-level normalization: quantile normalize all samples to the average empirical distribution observed across samples
    M = normalize_quantiles(expr.loc[mask].values, inplace=False)

    R = inverse_quantile_normalization(M)
    dtype = np.float32
    quant_std_df = pd.DataFrame(data=R, columns=expr.columns, index=expr.loc[mask].index, dtype = dtype)

    expr_filt = quant_std_df

    expr_filt.to_csv(f'{output}/{tissue}-normfilt.csv')
    print('writing preprocessed expression file') 

def normalize_quantiles(M, inplace=False):
    """
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")

    Reference:
        [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003

    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    if not inplace:
        M = M.copy()

    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n

    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1

        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1

    if not inplace:
        return M


# Gene-level normalization: quantile normalize to a standard normal distribution.
def inverse_quantile_normalization(M):
    """
    After quantile normalization of samples, standardize expression of each gene
    """
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q       


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tissue", required=True, help="Tissue of expression data")
    parser.add_argument("-o", "--output", required=True, help="Output folder to write filtered expression file")
    params = parser.parse_args()
    preprocess(params.tissue, params.output)


if __name__ == "__main__":
    main()
