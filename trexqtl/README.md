# trexqtl-run

trexqtl-run is a tool for to run Trex-QTL and/or CPMA on sorted snp-gene pairs association statistics to detect trans-eQTLs 

## Basic usage

The basic command to run trexqtl-run is:

```
trexqtl-run [options] --input input_file --out outfolder
```

This will run Trex-QTL and/or the CPMA method on sorted snp-gene pair association statistics file. 

## Detailed usage

Type `trexqtl-run --help` for a full list of available options. Basic options are described below:

### Run options
* `--trexqtl`: Run Trex-QTL
* `--cpma`: Run CPMA
* `--null_method` <string> : How to get the null distribution for test. Following statsOptions: ['chi2', 'eigen']
* `--precomputed_null` <string> : Precomputed file with null values for test stats. If user already has a precomputed (from previous runs or other methods) file with null Trex-QTL and/or CPMA null values to use, they can choose to use their file instead of recomputing the null values distribution. File must have the columns ['trexQTL'] and/or ['CPMA'] depending on which run options are chosen.
* `--threads` <int>: Number of threads to use


### Input options
* `--input <string>`: Matrix eQTL file. Must be sorted by SNP. Recommended tool to obtain association statistics file is Matrix-eQTL. Otherwise can also use other snp-gene association statistics file from other tools, but must have the columns ['SNP', 'gene', 't-stat'] Association statistics must be t-statistics. 
  
### Output options
* `--out <string>`: path to the output folder to store results

### Other options
* `--seed <int>`: Seed to use for random number generation
* `--version`: print the Trex-QTL version

## Output files

trexqtl-run will output a Trex-QTL and/or CPMA statistics file and a null values file if 'eigen' is chosen for `--null_method` in `$OUT` folder, where `$OUT` is the argument to `--out` 

* `results.tab`: Trex-QTL and/or CPMA statistics file with CPMA statistic, pvalue for CPMA statistic, Trex-QTL statistic, pvalue for Trex-QTL statistic, Trex-QTL predicted _t_, Trex-QTL predicted likelihood
* `null_sim.tab`: Computed file with null values for test stats if 'eigen' for `--null_method`is chosen. Saved null value file can be used for future runs using the `--precomputed_null` to save on compute resources. If a file is passed to `--precomputed_null`, no new null values file will be computed and saved.
  

