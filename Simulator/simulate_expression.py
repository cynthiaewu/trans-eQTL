import numpy as np
import pandas as pd
import yaml
from scipy.stats import bernoulli
import argparse


def generate_genotypes(sample_size, allele_freq):
    allele1 = bernoulli.rvs(allele_freq, size=sample_size)
    allele2 = bernoulli.rvs(allele_freq, size=sample_size)
    genotype = allele1 + allele2
    return genotype


def get_noise(num_genes, cov, sample_size):
    return [np.random.multivariate_normal(np.zeros(num_genes), cov)
           for i in range(sample_size)]


def write_xfile(array, output):
    try:
        shape = (array.shape[0], array.shape[1])
    except IndexError:
        shape = (array.shape[0], 1)
    col = ['Sample' + str(i) for i in range(shape[0])]
    row = ['SNP' + str(i) for i in range(shape[1])]
    df = pd.DataFrame(array.reshape(1, -1), index=row, columns=col)
    df.to_csv(f'{output}/genotypen.csv', index=True, header=True, sep='\t')


def write_yfile(array, output):
    col = ['Sample' + str(i) for i in range(len(array[0]))]
    row = ['Gene' + str(i) for i in range(len(array))]
    df = pd.DataFrame(array, index=row, columns=col)
    df.to_csv(f'{output}/expression.csv', index=True, header=True, sep='\t')


def model(config, beta_file, cov_file, seed, output):
    print(f'Seed = {seed}')
    np.random.seed(seed)

    with open(config) as f:
        params = yaml.load(f)
        print(params)
    beta = np.loadtxt(beta_file)
    cov = np.loadtxt(cov_file)
    num_genes = params['num_genes']
    allele_freq = params['allele_freq']
    sample_size = params['sample_size']

    X = generate_genotypes(sample_size=sample_size,
                                   allele_freq=allele_freq)
    noise = get_noise(num_genes=num_genes,
                      cov=cov,
                      sample_size=sample_size)

    Y = np.outer(beta, X)  + np.array(noise).T
    write_xfile(X, output)
    write_yfile(Y, output)
    print(f'Var: {np.var(Y)}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-b", "--beta_file", required=True, help="Input file with effect size distribution")
    parser.add_argument("-v", "--cov_file", required=True, help="Input file with gene covariance matrix")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Seed for random generator")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    model(config=params.config,
          beta_file=params.beta_file,
          cov_file=params.cov_file,
          seed=params.seed,
          output=params.output)


if __name__ == "__main__":
    main()
