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
    df.to_csv(f'{output}/genotype.csv', index=True, header=True, sep='\t')


def write_yfile(array, output):
    col = ['Sample' + str(i) for i in range(len(array[0]))]
    row = ['Gene' + str(i) for i in range(len(array))]
    df = pd.DataFrame(array, index=row, columns=col)
    df.to_csv(f'{output}/expression.csv', index=True, header=True, sep='\t')


def model(num_genes, allele_freq, sample_size, beta_file, cov, output):
    beta = np.loadtxt(beta_file)

    X = generate_genotypes(sample_size=sample_size,
                                   allele_freq=allele_freq)
    noise = get_noise(num_genes=num_genes,
                      cov=cov,
                      sample_size=sample_size)

    Y = np.outer(beta, X)  + np.array(noise).T
    write_xfile(X, output)
    write_yfile(Y, output)
    print(f'Var: {np.var(Y)}')


def iter_model(config, seed, iterations, output):
    print(f'Seed = {seed}')
    np.random.seed(seed)
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    num_genes = params['num_genes']
    allele_freq = params['allele_freq']
    sample_size = params['sample_size']
    identity = params['identity']
    if identity:
        cov = np.identity(num_genes)
        #cov_file = f'{output}cov.txt'
        #cov = np.loadtxt(cov_file)
        print('Identiy covariance matrix created')
    sim_prefix = 'Simulation'
    for i in range(iterations):
        folder = f'{output}{sim_prefix}_{i}/'
        if not identity:
            cov_file = f'{folder}cov.txt'
            cov = np.loadtxt(cov_file)
        model(num_genes, allele_freq, sample_size, f'{folder}beta.txt', cov, f'{folder}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    #parser.add_argument("-p", "--sim_prefix", required=True, help="Prefix for folders of simulated files (cov matrix and beta files)")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Seed for random generator")
    parser.add_argument("-i", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    iter_model(config=params.config,
          #sim_prefix=params.sim_prefix,
          seed=params.seed,
          iterations=params.iterations,
          output=params.output)


if __name__ == "__main__":
    main()
