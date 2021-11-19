#!/usr/bin/env python3

"""
Tool to simulate expression data

Example command:
TODO

"""

import argparse
import os
import sys
from simulate import __version__

def main(args):
	print("Not implemented")
	return 1

def getargs(): 
    parser = argparse.ArgumentParser(__doc__)
    sim_group = parser.add_argument_group("Simulation parameters")
    sim_group.add_argument("--num-genes", help="Number of genes to simulate", type=int, default=15000)
    sim_group.add_argument("--num-snps", help="Number of SNPs to simulate", type=int, default=1)
    sim_group.add_argument("--num-nullsnps", help="Number of null SNPs to simulate", type=int, default=0)
    sim_group.add_argument("--num-peer", help="Number of PEER factors to simulate", type=int, default=0)
    sim_group.add_argument("--num-samples", help="Number of individuals", type=int, default=500)
    sim_group.add_argument("--maf", help="Minor allele frequency for simulated SNPs", type=float, default=0.5)
    sim_group.add_argument("--num-targets", help="Number of target genes of the trans-eqtl", type=int, default=0)
    sim_group.add_argument("--gxg-corr", help="File with gene-gene correlation matrix", type=str, default=None)
    sim_group.add_argument("--beta", help="Mean effect size", type=float, default=0)
    sim_group.add_argument("--beta-sd", help="Std.dev of effect sizes", type=float, default=0)
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