import numpy as np
import random
from scipy.stats import truncnorm
import os
import argparse

    
def generate_cov(genes):
    matrix = np.random.uniform(0, 1, (genes, genes))
    symm_matrix = (matrix + matrix.T)/2
    np.fill_diagonal(symm_matrix, 1)
    return symm_matrix

def generate_beta(genes, num_targets):
    targets = truncnorm.rvs(0.5, 1, size=num_targets)
    sign = np.random.choice([-1, 1], num_targets)
    beta = np.concatenate([targets*sign, np.zeros(genes-num_targets)])
    return beta


def generator(genes, num_targets, output):
    cov_matrix = generate_cov(genes)
    np.savetxt(f'{output}/cov.txt', cov_matrix)
    beta = generate_beta(genes, num_targets)
    np.savetxt(f'{output}/beta.txt', beta)


def iter_generator(genes, num_targets, seed, iterations, output):
    print(f'Seed = {seed}')
    np.random.seed(seed)
    for i in range(iterations):
        output_path = os.path.join(output, f'Simulation_{i}')
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        generator(genes, num_targets, output_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genes", type=int, required=True, help="Number genes")
    parser.add_argument("-t", "--num_targets", type=int, required=True, help="Number target genes for the trans-eqtl")
    parser.add_argument("-s", "--seed", type=int, default=0, help="Seed for random generator")
    parser.add_argument("-i", "--iterations", type=int, default=0, help="# iterations for generating cov matrix and beta files")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    iter_generator(genes=params.genes,
          num_targets=params.num_targets,
          seed=params.seed,
          iterations=params.iterations,
          output=params.output)


if __name__ == "__main__":
    main()
