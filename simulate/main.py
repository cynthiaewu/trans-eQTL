#!/usr/bin/env python3

"""
Tool to simulate expression data

Example command:
TODO

"""

import argparse
import numpy as np
import os
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
	def __init__(self, cov, num_genes=DEFAULT_NUM_GENES, \
				num_samples=DEFAULT_NUM_SAMPLES, \
				num_snps=DEFAULT_NUM_SNPS, \
				num_nullsnps=DEFAULT_NUM_NULLSNPS, \
				num_peer=DEFAULT_NUM_PEER, \
				maf=DEFAULT_MAF, \
				num_targets=DEFAULT_NUM_TARGETS, \
				beta=DEFAULT_BETA, \
				beta_sd=DEFAULT_BETA_SD):
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
		ERROR("Not implemented")

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

	simulator = ExprSimulator(cov, \
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