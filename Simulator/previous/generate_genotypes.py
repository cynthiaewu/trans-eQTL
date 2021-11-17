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


def write_xfile(array, num_snps, output):
    if num_snps == 1:
        array = array.reshape(1, -1)
    
    sample_size = array.shape[1]
    col = ['Sample' + str(i) for i in range(sample_size)]
    row = ['SNP' + str(i) for i in range(num_snps)]
    df = pd.DataFrame(array, index=row, columns=col)
    df.to_csv(f'{output}/genotype_null.csv', index=True, header=True, sep='\t')


def iter_model(config, seed, iterations, output):
    print(f'Seed = {seed}')
    np.random.seed(seed)
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    allele_freq = params['allele_freq']
    sample_size = params['sample_size']
    num_snps = params['num_snps']
    sim_prefix = 'Simulation'
    for i in range(iterations):
        folder  = f'{output}/{sim_prefix}_{i}'
        cur_eqtl = np.loadtxt(f'{folder}/genotype.csv', dtype=int, skiprows=1, usecols=range(1,sample_size+1))
        cur_eqtl = cur_eqtl.reshape(1, -1)
        X = generate_genotypes(sample_size=sample_size, allele_freq=allele_freq, num_snps=num_snps-1)
        all_snps = np.concatenate((cur_eqtl, X))
        write_xfile(all_snps, num_snps, folder)
        print(f'Simulation {i}, folder {output}')
    print(f'Finished simulations for {output}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    #parser.add_argument("-n", "--noise", required=True, help="Input noise matrix to randomly draw noise from")
    #parser.add_argument("-p", "--sim_prefix", required=True, help="Prefix for folders of simulated files (cov matrix and beta files)")
    parser.add_argument("-s", "--seed", type=int, default=1, help="Seed for random generator")
    parser.add_argument("-i", "--iterations", type=int, required=True, help="# iterations to simulate genotype and expression files")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    iter_model(config=params.config,
          #noise_file=params.noise,
          seed=params.seed,
          iterations=params.iterations,
          output=params.output)


if __name__ == "__main__":
    main()
