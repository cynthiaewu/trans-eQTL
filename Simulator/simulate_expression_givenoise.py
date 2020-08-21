import numpy as np
import pandas as pd
import random
import yaml
from scipy.stats import bernoulli
#from scipy.stats import multivariate_normal
import argparse


def generate_genotypes(sample_size, allele_freq, num_snps):
    genotype = []
    for i in range(num_snps):
        allele1 = bernoulli.rvs(allele_freq, size=sample_size)
        allele2 = bernoulli.rvs(allele_freq, size=sample_size)
        genotype.append(allele1 + allele2)
        if i % 500 == 0:
            print(i)
    return np.array(genotype)


def get_noise(num_genes, noise_matrix, sample_size):
    choice = random.randint(0, len(noise_matrix)-1)
    noise = []
    for i in range(sample_size):
        noise.append(noise_matrix[choice])
    return noise


def write_xfile(array, num_snps, output):
    if num_snps == 1:
        array = array.reshape(1, -1)
    # try:
    # shape = array.shape
    # except IndexError:
    #     shape = (array.shape[0], 1)
    #print(shape)
    sample_size = array.shape[1]
    col = ['Sample' + str(i) for i in range(sample_size)]
    row = ['SNP' + str(i) for i in range(num_snps)]
    #print(col)
    #print(array)
    # if num_snps == 1:
    #     df = pd.DataFrame(array.reshape(1, -1), index=row, columns=col)
    # else:
    df = pd.DataFrame(array, index=row, columns=col)
    df.to_csv(f'{output}/genotype.csv', index=True, header=True, sep='\t')


def write_yfile(array, output):
    col = ['Sample' + str(i) for i in range(len(array[0]))]
    row = ['Gene' + str(i) for i in range(len(array))]
    df = pd.DataFrame(array, index=row, columns=col)
    df.to_csv(f'{output}/expression.csv', index=True, header=True, sep='\t')


def model(num_genes, allele_freq, sample_size, num_snps, beta_file, noise_matrix, output):
    
    beta = np.loadtxt(beta_file)
    print('starting generating genotypes')
    X = generate_genotypes(sample_size=sample_size, allele_freq=allele_freq, num_snps=num_snps)
    print('finished generating genotypes')
   #print(X)
    sum_X = 0
    print('started computing summation of beta and genotype')
    for i in range(num_snps):
        sum_X += (np.outer(beta[i], X[i]))
        if i % 500 == 0:
            print(i)
        #print(beta[i])
        #print(sum_X)
    #print(sum_X.shape)
    print('finished computing summation of beta and genotype')
    print('starting getting noise')
    
    noise = get_noise(num_genes, noise_matrix, sample_size)
    print('finished getting noise')
    #print(np.array(noise).shape)
    Y = sum_X + np.array(noise).T
    print('finished computing expression')
    #print(X.shape)
    #Y = np.outer(beta, X)  + np.array(noise).T
    write_xfile(X, num_snps, output)
    write_yfile(Y, output)
    #print(f'Var: {np.var(Y)}')


def iter_model(config, noise_file, seed, iterations, output):
    print(f'Seed = {seed}')
    np.random.seed(seed)
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    num_genes = params['num_genes']
    allele_freq = params['allele_freq']
    sample_size = params['sample_size']
    num_snps = params['num_snps']
    identity = params['identity']
    noise_matrix = np.loadtxt(noise_file)
    if identity:
        cov = np.identity(num_genes)
        #cov_file = f'{output}cov.txt'
        #cov = np.loadtxt(cov_file)
        print('Identiy covariance matrix created')
    sim_prefix = 'Simulation'
    for i in range(iterations):
        folder = f'{output}/{sim_prefix}_{i}/'
        if not identity:
            cov_file = f'{folder}cov.txt'
            cov = np.loadtxt(cov_file)
        model(num_genes, allele_freq, sample_size, num_snps,  f'{folder}beta.txt', noise_matrix, f'{folder}')
        print(f'Simulation {i}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-n", "--noise", required=True, help="Input noise matrix to randomly draw noise from")
    #parser.add_argument("-p", "--sim_prefix", required=True, help="Prefix for folders of simulated files (cov matrix and beta files)")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Seed for random generator")
    parser.add_argument("-i", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    iter_model(config=params.config,
          noise_file=params.noise,
          seed=params.seed,
          iterations=params.iterations,
          output=params.output)


if __name__ == "__main__":
    main()
