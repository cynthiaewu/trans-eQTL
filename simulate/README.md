# xQTL-simulate

xQTL-simulate is a tool for generated simulated gene expression and genotype datasets with trans-eQTL effects

## Basic usage

The basic command to run xQTL-simulate is:

```
xQTL-simulate [options] --out outfolder
```

This will generate simulation output files in folders within `outfolder`. By default, one simulation will be performed to generate data for 500 samples at 1 SNP and 15,000 genes with no trans-eQTL effects. See below for detailed usage info and a description of output files.

## Detailed usage

Type `xQTL-simulate --help` for a full list of available options. Basic options are described below:

### Simulation options
* `--num-genes <int>`: The number of genes to simulate expression for
* `--num-samples <int>`: The sample size of the simulated dataset
* `--num-snps <int>`: The number of SNPs with trans-eqtl effects to simulate
* `--num-nullsnps <int>`: The number of SNPs with no effects to simulate
* `--maf <float>`: The minor allele frequency of SNPs
* `--num-targets <int>`: The number of target genes of each trans-eqtl
* `--num-peer <int>`: Number of PEER-like technical covariate factors to simulate.
* `--gxg-corr <file>`: Path to gene-by-gene correlation matrix to use. If not specified, the identity matrix is used, indicating no gene by gene correlation.
* `--beta <float>`: Mean effect size of each QTL
* `--beta-sd <float>`: Standard deviation of effect sizes

### Output options
* `--out <string>`: path to the output folder to store results

### Other options
* `--num-sim <int>`: number of simulation rounds to perform
* `--seed <int>`: Seed to use for random number generation
* `--version`: print the xQTL version

## Output files

xQT-simulate will output a separate folder `$OUT/sim-$i` for each simulation round, where `$OUT` is the argument to `--out` and `$i` is the simulation number, starting at 0.

Each folder contains:
* `genotypes.csv`: A csv file with one row per SNP and one column per sample. SNP genotypes are 0, 1, or 2.
* `expression.csv`: A csv file with one row per gene and one column per sample. Values give the expression level of each sample/gene pair.
* `betas.csv`: A csv file wtih one row per SNP and one column per gene. Values give the effect size for each SNP/gene pair.

