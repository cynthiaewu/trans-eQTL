import numpy as np
import random
from scipy.stats import truncnorm
import os
import yaml
import time
import argparse

    
def generate_cov(num_genes):
    matrix = np.random.uniform(0, 1, (num_genes, num_genes))
    symm_matrix = (matrix + matrix.T)/2
    np.fill_diagonal(symm_matrix, 1)
    return symm_matrix


def generate_identity(num_genes):
    return np.identity(num_genes)


def generate_beta(num_genes, targets, fixed_betas):
    
    all_betas = []
    #t3 = time.time()
    #print(f'Beta zeros created: {t3-t2}')
    
    #targets = np.random.normal(0, 0.1, num_targets)
    #targets = [20, 40, 50, 60, 80, 100, 150, 200]
    #targets = [400, 600, 800, 1000, 1250, 1500, 1750, 2000]
    #targets = [250, 300, 350]
    #fixed_betas = [0.1, 0.2, 0.3, 0.4, 0.5]
    for num_t in targets:
        for beta_value in fixed_betas: 
            beta = np.zeros(num_genes)
            values = np.full(num_t, beta_value) 
            beta[:num_t] = values
            # iterations
            for i in range(100):
                all_betas.append(beta)

    return all_betas


def generator(num_genes, num_targets, identity, num_snps, num_nullsnps, beta_sd, beta_value, output):
    #t0 = time.time()
    if not identity:
        cov_matrix = generate_cov(num_genes)
        np.savetxt(f'{output}/cov.txt', cov_matrix)
        #t1 = time.time()
        #print(f'Cov matrix created: {t1-t0}')
        #t0 = time.time()
        #print(f'Cov matrix saved: {t0-t1}')
    beta = generate_beta(num_genes, num_targets, beta_value)
    np.savetxt(f'{output}/beta.txt', beta)
    #t6 = time.time()


def iter_generator(config, seed, iterations, output):
    print(f'Seed = {seed}')
    np.random.seed(seed)
    with open(config) as f:
        params = yaml.load(f)
    num_genes = params['num_genes']
    num_targets = params['num_targets']
    num_snps = params['num_snps']
    num_nullsnps = params['num_nullsnps']
    identity = params['identity']
    beta_sd = params['beta_sd']
    beta_value = params['beta_value']
    #if identity:
    #    cov_matrix = generate_identity(num_genes)
    #    np.savetxt(f'{output}/cov.txt', cov_matrix)
        

    for i in range(iterations):
        output_path = os.path.join(output, f'Simulation_{i}')
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        generator(num_genes, num_targets, identity, num_snps, num_nullsnps,  beta_sd, beta_value, output_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    #parser.add_argument("-g", "--genes", type=int, required=True, help="Number genes")
    #parser.add_argument("-t", "--num_targets", type=int, required=True, help="Number target genes for the trans-eqtl")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Seed for random generator")
    parser.add_argument("-i", "--iterations", type=int, default=0, help="# iterations for generating cov matrix and beta files")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    iter_generator(config=params.config,
          seed=params.seed,
          iterations=params.iterations,
          output=params.output)


if __name__ == "__main__":
    main()
