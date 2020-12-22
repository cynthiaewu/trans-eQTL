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


def generate_beta(num_genes, targets, num_nullsnps, beta_type, fixed_betas, rep):
    
    all_betas = []

    # simulate beta values drawn from normal distribution given standard deviation as input
    if beta_type == 'sd':
        for num_t in targets:
            for beta_sd in fixed_betas: 
                values = np.random.normal(0, beta_sd, num_t)
                beta = np.zeros(num_genes)
                beta[:num_t] = values
            # iterations
                for i in range(rep):
                    all_betas.append(beta)
    # simulate same fixed beta values for all target genes given beta value as input 
    if beta_type == 'value':
        for num_t in targets:
            for beta_value in fixed_betas: 
                beta = np.zeros(num_genes)
                values = np.full(num_t, beta_value) 
                beta[:num_t] = values
            # iterations
                for i in range(rep):
                    all_betas.append(beta)

    for i in range(num_nullsnps):
        all_betas.append((np.zeros(num_genes)))
    return all_betas


def generator(num_genes, num_targets, identity, num_snps, num_nullsnps, beta, beta_sd, beta_value, output_path, output, rep):
    if not identity:
        cov_matrix = generate_cov(num_genes)
        np.savetxt(f'{output}/cov.txt', cov_matrix)
    if beta == 'sd':
        beta = generate_beta(num_genes, num_targets,num_nullsnps,  beta, beta_sd, rep)
        np.savetxt(f'{output_path}/beta.txt', beta)
    if beta == 'value':
        beta = generate_beta(num_genes, num_targets, num_nullsnps, beta, beta_value, rep)
        np.savetxt(f'{output}/beta.txt', beta)


def iter_generator(config, seed, iterations, output):
    #print(f'Seed = {seed}')
    np.random.seed(seed)
    with open(config) as f:
        params = yaml.load(f)
    num_genes = params['num_genes']
    num_targets = params['num_targets']
    num_snps = params['num_snps']
    num_nullsnps = params['num_nullsnps']
    identity = params['identity']
    beta = params['beta']
    beta_sd = params['beta_sd']
    beta_value = params['beta_value']
    rep = params['rep']
    # don't actually have to create and save identity cov matrix, can just call numpy to create a identity matrix during calculations
    # faster than reading and writing the files for identity cov matrix

    #if identity:
    #    cov_matrix = generate_identity(num_genes)
    #    np.savetxt(f'{output}/cov.txt', cov_matrix)

    #print(f'Starting config generator for #targets: {num_targets} and beta: {beta_value}')
        

    for i in range(iterations):
        output_path = os.path.join(output, f'Simulation_{i}')
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        generator(num_genes, num_targets, identity, num_snps, num_nullsnps, beta, beta_sd, beta_value, output_path, output, rep)
    print(f'Finished config generator for #targets: {num_targets} and beta: {beta_sd}')
   

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
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
