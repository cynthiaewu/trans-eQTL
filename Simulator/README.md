# Simulator
1. Generate config files with a given meta config yaml file

Iterations are the number of desired simulations. A folder will be created for each iteration/simulation and will contain "beta.txt" which contains the beta effect sizes each snp has on each gene and also a "cov.txt" file if the covariance/gene_correlation file is not the identity matrix.
```
python config_generator.py -c metaconfig_yaml_file -i iterations -o output_folder
```

Used for cpma method to have a simulation file with all trans-eQTLs with varying #target genes in one simulation file. The 'num_targets' in the meta_config_yaml file contains a list of #target genes to test and the 'beta_value' contains a list of the beta values to test. 100 iterations (to create necessary config files) for each pair of #target genes and beta value will be performed. 
```
config_generator_specific.py -c metaconfig_yaml_file -i iterations -o output_folder
```
2. Simulate genotype and expression files given that config generator has already been run to generate necessary config files. The metaconfig_yaml_file, iterations, and output_folder should be the same as (1.)
```
python simulate_expression.py -c metaconfig_yaml_file -i iterations -o output_folder
```
3. Run eQTL methods with the simulated genotype and expression files.
4. Calculate power of the results from the eQTL methods.

For Matrix eQTL method, calculate power across the Matrix eQTL results of simulations. Simulation folder contains all the Matrix eQTL results.
```
python calculate_power_across_sim.py -c metaconfig_yaml_file -f simulation_folder -i iterations -o output_folder
```
For CPMA method which tests one trans eQTL among many null snps (snps with beta=0 across all genes), calculate power across cpma results of simulations.
```
python calculate_power_across_sim_cpma.py -c metaconfig_yaml_file -f simulation_folder -i iterations -o output_folder
```
For CPMA method which tests many trans eQTLs in one simulation with no null snps, calculate power across cpma results of each pair of #target genes and beta value parameter. 

Before calculating power, obtain the pvalue from each cpma value for each snp by comparing the cpma value to the chi square distribution. pvalues are adjusted for the number of tests performed per pair of #target genes and beta value parameter. 
```
python compute_pvalue_chidist.py -i input_file -n num_test -o output_file
```
Input file of 'calculate_power_manysiminone_cpma.py' should be the output file from 'compute_pvalue_chidist.py'
```
python calculate_power_manysiminone_cpma.py -c metaconfig_yaml_file -i input_file -o output_folder
```
