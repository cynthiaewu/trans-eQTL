import numpy as np
import yaml
from scipy.stats import bernoulli
import argparse


def generate_genotype(allele_freq):
    allele1 = bernoulli.rvs(allele_freq, size=1)
    allele2 = bernoulli.rvs(allele_freq, size=1)
    genotype = allele1 + allele2
    return int(genotype)


def model(config, beta_file, cov_file, geno_file, seed, output):
    with open(config) as f:
        params = yaml.load(f)
        print(params)
    beta = np.loadtxt(beta_file)
    cov = np.loadtxt(cov_file)
    genotype = np.loadtxt(geno_file)
    num_genes = params['num_genes']
    allele_freq = params['allele_freq']
    sample_size = params['sample_size']
    genotypes = []
    noise = []
    for i in range(sample_size):
        genotypes.append(generate_genotype(allele_freq))
        noise.append(np.random.multivariate_normal(np.zeros(num_genes), cov))
    genotypes = np.array(genotypes).reshape(1, -1)
    beta = beta.reshape(-1, 1)
    #print(beta.shape)
    #print(genotypes.shape)
    X = np.dot(beta,genotypes)
    print(X.shape)
    print((np.array(noise)).shape)
    Y = X + np.array(noise).T
    print(np.var(Y))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-b", "--beta_file", required=True, help="Input file with effect size distribution")
    parser.add_argument("-v", "--cov_file", required=True, help="Input file with gene covariance matrix")
    parser.add_argument("-g", "--geno_file", required=True, help="Input file with genotype matrix")
    parser.add_argument("-s", "--seed", required=True, help="Seed for random generator")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with simulated expression matrix")
    params = parser.parse_args()
    model(params.config, params.beta_file, params.cov_file, params.geno_file, params.seed, params.output)


if __name__ == "__main__":
    main()
