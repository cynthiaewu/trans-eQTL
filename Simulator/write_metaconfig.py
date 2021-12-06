import yaml
import argparse


def write_metaconfigs(input, num_targets, num_genes, num_snps, num_nullsnps, num_factors, samplesize,
        allele_freq, identity, beta, beta_value, sig_threshold):
    targets = [int(x) for x in num_targets.split(',')]
    beta_values = [float(x) for x in beta_value.split(',')]
    #beta_values = beta_values_str.split(',')
    #targets = [ 0, 20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500, 700, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]
    #beta_values = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]
    for tar in targets:
        for b in beta_values:
            parameters = {'num_genes': num_genes,
                    'num_snps': num_snps,
                    'num_nullsnps': num_nullsnps,
                    'num_factors': num_factors,
                    'sample_size': samplesize,
                    'allele_freq': allele_freq,
                    'num_targets': [tar],
                    'identity': identity,
		            'beta': beta,
                    #'beta_sd': 'NA',
                    'beta_value': [b],
                    'sig_threshold': sig_threshold}
            value = str(b).replace(".","")
            if (b.is_integer()):
                value = str(int(b))
            else:
                value = str(b).replace(".","")

            filename = f'{input}/numTarget_{tar}/Beta_{value}/metaconfig.yaml'
            #filename = f'/storage/cynthiawu/trans_eQTL/Scripts/Test_nullsnps/simulate_eqtls_only/FastMultivariate/Single_eqtl/SampleSize100/SingleParameter/numTarget_{tar}/Beta_{value}/metaconfig.yaml'
            with open(filename, 'w') as file:
                documents = yaml.dump(parameters, file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input folder with folders of targets and beta values")
    parser.add_argument("-t", "--num_targets", required=True, help="Input # target genes")
    parser.add_argument("-g", "--num_genes", type=int, default=15000, help="# of total genes")
    parser.add_argument("-x", "--num_snps", type=int, default=1, help="# of total snps")
    parser.add_argument("-n", "--num_nullsnps", type=int, default=0, help="# of null snps with no beta effect size on all genes")
    parser.add_argument("-f", "--num_factors", type=int, default=0, help="# of PEER factors")
    parser.add_argument("-s", "--samplesize", type=int, default=500, help="Sample size")
    parser.add_argument("-a", "--allele_freq", type=float, default=0.5, help="allele frequency of snps")
    parser.add_argument("-c", "--correlation", action='store_false', help="Include gene correlation (covariance matrix will be the gene correlation matrix of Nerve Tibial GTeX data")
    parser.add_argument("-b", "--beta", type=str, default='value', help="'value' or 'sd'")
    #parser.add_argument("-d", "--beta_sd", type=float, default=-1, help="standard deviation value of normal distribution with mean 0 to draw beta effect size values for the target genes")
    parser.add_argument("-v", "--beta_value", required=True, help="beta effect size for all target genes, either fixed value or standard deviation value")
    parser.add_argument("-p", "--sig_threshold", type=float, default=0.05, help="significance threshold for calculating the power of the eqtl methods")
    params = parser.parse_args()
    
    write_metaconfigs(input=params.input,
          num_targets=params.num_targets,
          num_genes=params.num_genes,
          num_snps=params.num_snps,
          num_nullsnps=params.num_nullsnps,
          num_factors=params.num_factors,
          samplesize=params.samplesize,
          allele_freq=params.allele_freq,
          identity=params.correlation,
          beta=params.beta,
          #beta_sd=params.beta_sd,
          beta_value=params.beta_value,
          sig_threshold=params.sig_threshold)


if __name__ == "__main__":
    main()
