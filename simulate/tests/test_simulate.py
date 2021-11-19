import argparse
import os,sys
import pytest

from ..main import *

@pytest.fixture
def args(tmpdir):
	args = argparse.ArgumentParser()
	args.num_genes = DEFAULT_NUM_GENES
	args.num_snps = DEFAULT_NUM_SNPS
	args.num_nullsnps = DEFAULT_NUM_NULLSNPS
	args.num_peer = DEFAULT_NUM_PEER
	args.num_samples = DEFAULT_NUM_SAMPLES
	args.maf = DEFAULT_MAF
	args.num_targets = DEFAULT_NUM_TARGETS
	args.gxg_corr = None
	args.beta = DEFAULT_BETA
	args.beta_sd = DEFAULT_BETA_SD
	args.num_sim = 1
	args.seed = 10
	args.out = ""
	return args

# TODO this is just an example test
def test_GenerateNoise(args):
	sim = ExprSimulator(num_genes=10, num_samples=5)
	noise = sim.generate_noise()
	assert noise.shape[0] == 5
	assert noise.shape[1] == 10
