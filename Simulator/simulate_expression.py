import numpy as np
import pandas as pd
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
    return np.array(genotype)


def generate_peer(num_factors, sample_size, num_genes):
    factor_level = np.random.normal(0, 0.6, (sample_size, num_factors))
    factor_weight = []
    for i in range(num_factors):
        var_k = 0.8*((np.random.gamma(2.5, 0.6))**2)
        factor_weight.append(np.random.normal(0, var_k, num_genes))
    factor_weight = np.array(factor_weight)
    peer = np.dot(factor_level, factor_weight)
    print(f'Simulated {num_factors} PEER factor')
    return peer


def get_noise(num_genes, cov, sample_size, identity):
    if identity:
        print('Simulating noise (identity)')
        return np.random.normal(0, 1, (sample_size, num_genes))
    else:
        print('Simulating noise (gene correlation)')
        return np.random.multivariate_normal(np.zeros(num_genes), cov, size=sample_size)


def write_xfile(array, num_snps, output):
    if num_snps == 1:
        array = array.reshape(1, -1)
    
    sample_size = array.shape[1]
    col = ['Sample' + str(i) for i in range(sample_size)]
    row = ['SNP' + str(i) for i in range(num_snps)]
    df = pd.DataFrame(array, index=row, columns=col)
    df.to_csv(f'{output}/genotype.csv', index=True, header=True, sep='\t')


def write_yfile(array, output):
    col = ['Sample' + str(i) for i in range(len(array[0]))]
    row = ['Gene' + str(i) for i in range(len(array))]
    df = pd.DataFrame(array, index=row, columns=col)
    df.to_csv(f'{output}/expression.csv', index=True, header=True, sep='\t')


def model(num_genes, allele_freq, sample_size, num_snps, beta_file, num_factors, noise, output):
    
    beta = np.loadtxt(beta_file)
    
    if num_snps == 1:
        beta = beta.reshape(1, -1)
    
    #print('starting generating genotypes')
    X = generate_genotypes(sample_size=sample_size, allele_freq=allele_freq, num_snps=num_snps)
    #print('finished generating genotypes')
    sum_X = 0
    #print('started computing summation of beta and genotype')
    for i in range(num_snps):
        sum_X += (np.outer(beta[i], X[i]))
    
    #sum_X = np.dot(X.T, beta).T
    if num_factors != 0:
        peer = generate_peer(num_factors, sample_size, num_genes).T
        Y = sum_X + peer + np.array(noise).T
    else:
        Y = sum_X + np.array(noise).T
    #print('finished computing summation of beta and genotype')
    #print('starting getting noise')
    #noise = get_noise(num_genes=num_genes, cov=cov, sample_size=sample_size)
    #print('finished getting noise')
    #print('finished computing expression')
    #Y = np.outer(beta, X)  + np.array(noise).T
    write_xfile(X, num_snps, output)
    write_yfile(Y, output)


def iter_model(config, seed, iterations, output):
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
    num_factors = params['num_factors']
    beta = params['beta']
    if identity:
        cov = np.identity(num_genes)

        # simulate all noise for all iterations at once (numpy multivariate normal only has to run SD once)
        # Can run into memory error for large sample sizes (samplesize=500)
        #noise = get_noise(num_genes=num_genes, cov=cov, sample_size=sample_size*iterations)
        #print('Finished multivariate normal')
    sim_prefix = 'Simulation'
    for i in range(iterations):
        folder = f'{output}/{sim_prefix}_{i}/'
        if not identity:
            cov_file = f'{folder}cov.txt'
            cov = np.loadtxt(cov_file)
        if beta == 'sd':
            beta_location = f'{folder}/beta.txt'
        if beta == 'value':
            beta_location = f'{output}/beta.txt'
        
        noise_iter = get_noise(num_genes=num_genes, cov=cov, sample_size=sample_size, identity=identity)
        #noise_iter = noise[i*sample_size:(i*sample_size+sample_size),]
        #model(num_genes, allele_freq, sample_size, num_snps,  f'{folder}/beta.txt', f'{folder}')
        model(num_genes, allele_freq, sample_size, num_snps, beta_location, num_factors, noise_iter, f'{folder}')
        print(f'Simulation {i}, folder {output}')
    print(f'Finished simulations for {output}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Seed for random generator")
    parser.add_argument("-i", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    iter_model(config=params.config,
          seed=params.seed,
          iterations=params.iterations,
          output=params.output)


if __name__ == "__main__":
    main()
