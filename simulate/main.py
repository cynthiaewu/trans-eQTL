#!/usr/bin/env python3

"""
Tool to simulate expression data

Example command:
xQTL-simulate --out test --num-samples 10 \
     --num-genes 20 --num-nullsnps 2 \
     --beta 0.2 --beta-sd 0 --num-targets 2 \
     --run-matrix-eqtl
 """

import argparse
import numpy as np
import os
from rpy2 import robjects
from scipy.stats import bernoulli
import sys
from simulate import __version__

###########
# Default values

DEFAULT_NUM_GENES = 15000
DEFAULT_NUM_SAMPLES = 500
DEFAULT_NUM_SNPS = 1
DEFAULT_NUM_NULLSNPS = 0
DEFAULT_NUM_PEER = 0
DEFAULT_MAF = 0.5
DEFAULT_NUM_TARGETS = 0
DEFAULT_BETA = 0
DEFAULT_BETA_SD = 0

###########

def ERROR(msg):
	sys.stderr.write("[ERROR]: " + msg.strip() + "\n")
	sys.exit(1)

def MSG(msg):
	sys.stderr.write("[PROGRESS]: " + msg.strip() + "\n")

class ExprSimulator:
	def __init__(self, cov="identity", num_genes=DEFAULT_NUM_GENES, \
				num_samples=DEFAULT_NUM_SAMPLES, \
				num_snps=DEFAULT_NUM_SNPS, \
				num_nullsnps=DEFAULT_NUM_NULLSNPS, \
				num_peer=DEFAULT_NUM_PEER, \
				maf=DEFAULT_MAF, \
				num_targets=DEFAULT_NUM_TARGETS, \
				beta=DEFAULT_BETA, \
				beta_sd=DEFAULT_BETA_SD):
		if cov == "identity":
			cov = np.identity(num_genes)
		self.cov = cov
		self.num_genes = num_genes
		self.num_samples = num_samples
		self.num_snps = num_snps
		self.num_nullsnps = num_nullsnps
		self.num_peer = num_peer
		self.maf = maf
		self.num_targets = num_targets
		self.beta = beta
		self.beta_sd = beta_sd

	def Simulate(self, output_folder):
		# Get noise (numsamp x numgene matrix)
		noise = self.generate_noise()

		# Get genotypes
		X = self.generate_genotypes()

		# Get PEER
		peer = self.generate_peer()

		# Get matrix of effect sizes
		betas = self.generate_effects()

		# Get genetic effect. Assume effects of
		# each SNP are additive and independent
		sum_X = 0
		for i in range(self.num_snps):
			sum_X += (np.outer(betas[i], X[i]))
 
		# Put it together in one model
		Y = sum_X + peer + np.array(noise).T

		# Write to a file
		self.write_sim(X, Y, betas, output_folder)

	def generate_noise(self):
		return np.random.multivariate_normal(np.zeros(self.num_genes), \
        	self.cov, size=self.num_samples)

	def generate_effects(self):
		num_total_snps = self.num_snps+self.num_nullsnps
		betas = np.zeros((num_total_snps, self.num_genes))
		for i in range(self.num_snps):
			for j in range(self.num_targets):
				betas[i,j] = np.random.normal(self.beta, self.beta_sd)
		return betas

	def write_sim(self, X, Y, betas, output_folder):
		# Write genotypes
		outx = open(os.path.join(output_folder, "genotypes.csv"), "w")
		header = ["SNP"] + ["Sample%s"%i for i in range(self.num_samples)]
		outx.write(",".join(header)+"\n")
		for i in range(self.num_snps):
			outitems = ["SNP%s"%i]
			for j in range(self.num_samples):
				outitems.append(str(X[i, j]))
			outx.write(",".join(outitems)+"\n")
		outx.close()

		# Write phenotypes
		outy = open(os.path.join(output_folder, "expression.csv"), "w")
		header = ["Gene"] + ["Sample%s"%i for i in range(self.num_samples)]
		outy.write(",".join(header)+"\n")
		for i in range(self.num_genes):
			outitems = ["Gene%s"%i]
			for j in range(self.num_samples):
				outitems.append(str(Y[i,j]))
			outy.write(",".join(outitems)+"\n")
		outy.close()

		# Write betas
		outb = open(os.path.join(output_folder, "betas.csv"), "w")
		header = ["SNP"] + ["Gene%s"%i for i in range(self.num_genes)]
		outb.write(",".join(header)+"\n")
		for i in range(self.num_snps+self.num_nullsnps):
			outitems = ["SNP%s"%i]
			for j in range(self.num_genes):
				outitems.append(str(betas[i,j]))
			outb.write(",".join(outitems)+"\n")
		outb.close()

	def generate_genotypes(self):
		genotype = []
		for i in range(self.num_snps+self.num_nullsnps):
			allele1 = bernoulli.rvs(self.maf, size=self.num_samples)
			allele2 = bernoulli.rvs(self.maf, size=self.num_samples)
			genotype.append(allele1 + allele2)
		return np.array(genotype)

	def generate_peer(self):
		factor_level = np.random.normal(0, 0.6, (self.num_samples, self.num_peer))
		factor_weight = []
		for i in range(self.num_peer):
			var_k = 0.8*((np.random.gamma(2.5, 0.6))**2)
			factor_weight.append(np.random.normal(0, var_k, self.num_genes))
		factor_weight = np.array(factor_weight)
		peer = np.dot(factor_level, factor_weight).T
		return peer

def RunMatrixEQTL(output_folder):
	genotype = os.path.join(output_folder, "genotypes.csv")
	expression = os.path.join(output_folder, "expression.csv")
	eqtl_file = os.path.join(output_folder, "matrixEQTL.out")
	cmd = '''
library(MatrixEQTL);
useModel = modelLINEAR;
pvOutputThreshold = 1;
errorCovariance = numeric();
snps = SlicedData$new();
snps$fileDelimiter = ",";      # the TAB character
snps$fileOmitCharacters = "None"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile("%s");
gene = SlicedData$new();
gene$fileDelimiter = ",";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile("%s");
me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  output_file_name = "%s",
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
		'''%(genotype, expression, eqtl_file)
	print(cmd)
	robjects.r(cmd)

def main(args):
	np.random.seed(args.seed)

	# Set up covariance matrix
	if args.gxg_corr is None:
		cov = np.identity(args.num_genes)
	else:
		if not os.path.exists(args.gxg_corr):
			ERROR("Could not find file %s"%args.gxg_corr)
		cov = np.loadtxt(args.gxg_corr)

	# Set up output directory
	if not os.path.exists(args.out):
		os.mkdir(args.out)

	simulator = ExprSimulator(cov=cov, \
		num_genes=args.num_genes, \
		num_samples=args.num_samples, \
		num_snps=args.num_snps, \
		num_nullsnps=args.num_nullsnps, \
		num_peer=args.num_peer, \
		maf=args.maf, \
		num_targets=args.num_targets, \
		beta=args.beta, \
		beta_sd=args.beta_sd)

	for numsim in range(args.num_sim):
		MSG("Running simulation %s"%numsim)
		output_folder = os.path.join(args.out, "sim-%s"%numsim)
		if not os.path.exists (output_folder):
			os.mkdir(output_folder)
		simulator.Simulate(output_folder)
		if args.run_matrix_eqtl:
			RunMatrixEQTL(output_folder)

	MSG("Done!")
	return 0

def getargs(): 
    parser = argparse.ArgumentParser(__doc__)
    sim_group = parser.add_argument_group("Simulation parameters")
    sim_group.add_argument("--num-genes", help="Number of genes to simulate", type=int, default=DEFAULT_NUM_GENES)
    sim_group.add_argument("--num-snps", help="Number of SNPs to simulate", type=int, default=DEFAULT_NUM_SNPS)
    sim_group.add_argument("--num-nullsnps", help="Number of null SNPs to simulate", type=int, default=DEFAULT_NUM_NULLSNPS)
    sim_group.add_argument("--num-peer", help="Number of PEER factors to simulate", type=int, default=DEFAULT_NUM_PEER)
    sim_group.add_argument("--num-samples", help="Number of individuals", type=int, default=DEFAULT_NUM_SAMPLES)
    sim_group.add_argument("--maf", help="Minor allele frequency for simulated SNPs", type=float, default=DEFAULT_MAF)
    sim_group.add_argument("--num-targets", help="Number of target genes of the trans-eqtl", type=int, default=DEFAULT_NUM_TARGETS)
    sim_group.add_argument("--gxg-corr", help="File with gene-gene correlation matrix", type=str, default=None)
    sim_group.add_argument("--beta", help="Mean effect size", type=float, default=DEFAULT_BETA)
    sim_group.add_argument("--beta-sd", help="Std.dev of effect sizes", type=float, default=DEFAULT_BETA_SD)
    sim_group.add_argument("--num-sim", help="Number of iterations of the simulation to run", type=int, default=1)
    sim_group.add_argument("-s", "--seed", help="Seed for random generator", type=int, default=0, )
    other_group = parser.add_argument_group("Other options")
    other_group.add_argument("--run-matrix-eqtl", help="Run matrix eQTL on output", action="store_true")
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--out", help="Output prefix", type=str, required=True)
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def run():
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__":
    run()