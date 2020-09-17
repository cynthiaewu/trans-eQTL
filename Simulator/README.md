# To run all steps below in one script

For CPMA method, run_simulate_cpma_pipeline.py Input folder is where all the simulations and cpma files are to be written. Scripts folder is the path to all the trans eqtl scripts. Sample size is the desired sample size to be used in the simulations. 
 ```
 python run_simulate_cpma_pipeline.py -i input_folder -p scripts_folder -s samplesize
 ```
 For CPMA_x method, **NOT FINISHED**
 ```
 run_simulate_cpma_topx_pipeline.py
 ```
# Simulator
1. Generate metaconfig yaml file to be used for simulations. Input folder must contain the numTarget_x and Beta_x folders. Need to manually edit the parameters desired for simulations in this script. The parameter beta needs to be either "value" or "sd". If it is set to be "value," the parameter beta_value needs to be set to the desired beta effect size value which will be used for all target genes. If the parameter beta is set to be "sd," the beta effect size vallue for the target genes will be drawn from a distribution with mean 0 and SD set in the parameter beta_sd. The parameter num_targets is the #target genes the eqtl is to target. The parameter num_snps is the # of eqtls desired in the simulation. The parameter num_nullsnps is the # of null snps with no beta effect size desired in the simulation. The parameter sample_size is the sample size for the genotype and expression files in the simulations. The parameter allele_freq is used for creating the genotype files for the snps. The parameter identity is set to be True for using the identity matrix to simulate the noise matrix used for simulations or False if a gene covariance matrix is given to simulate the noise matrix. The parameter rep is used for "manysiminone" simulations where one simulation contains multiple eqtls with different parameters. The parameter sig_threshold is used to set the significance for calculating the power of the eqtl methods. **NEED MORE WORK**
   ```
   python write_metaconfig.py -i input_folder -s samplesize
   ```
   
2. Generate config files with a given metaconfig yaml file

   Iterations are the number of desired simulations. A folder will be created for each iteration/simulation and will contain "beta.txt" which contains the beta effect sizes each snp has on each gene and also a "cov.txt" file if the covariance/gene_correlation file is not the identity matrix.
   ```
   python config_generator.py -c metaconfig_yaml_file -i iterations -o output_folder
   ```

   Used for cpma method to have a simulation file with all trans-eQTLs with varying #target genes in one simulation file. The 'num_targets' in the meta_config_yaml file contains a list of #target genes to test and the 'beta_value' contains a list of the beta values to test. 100 iterations (to create necessary config files) for each pair of #target genes and beta value will be performed. 
   ```
   config_generator_specific.py -c metaconfig_yaml_file -i iterations -o output_folder
   ```
3. Simulate genotype and expression files given that config generator has already been run to generate necessary config files. The metaconfig_yaml_file, iterations, and output_folder should be the same as (1.)
   ```
   python simulate_expression.py -c metaconfig_yaml_file -i iterations -o output_folder
   ```
4. Run eQTL methods with the simulated genotype and expression files.
5. Calculate power of the results from the eQTL methods.

   For Matrix eQTL method, calculate power across the Matrix eQTL results of simulations. Simulation folder contains all the Matrix eQTL results.
   ```
   python calculate_power_across_sim.py -c metaconfig_yaml_file -f simulation_folder -i iterations -o output_folder
   ```
   For CPMA method which tests one trans eQTL among many null snps (snps with beta=0 across all genes), calculate power across cpma results of simulations.
   ```
   python calculate_power_across_sim_cpma.py -c metaconfig_yaml_file -f simulation_folder -i iterations -o output_folder
   ```
   For CPMA method which tests many trans eQTLs in one simulation with no null snps, calculate power across cpma results of each pair of #target genes and beta value parameter. (The CPMAx method currently only use the pvalues of the top 10% genes for each snp for its calculation.)

   Before calculating power, obtain the pvalue from each cpma value for each snp by comparing the cpma value to the chi square distribution. pvalues are adjusted for the number of tests performed per pair of #target genes and beta value parameter. For cpma_type parameter, 0 is for cpma, 1 is for cpmax. 
   ```
   python compute_pvalue_chidist.py -i input_file -c cpma_type -n num_test -o output_file
   ```
   Input file of 'calculate_power_manysiminone_cpma.py' should be the output file from 'compute_pvalue_chidist.py' The script is for calculating power for many simulations done in a single file. 
   ```
   python calculate_power_manysiminone_cpma.py -c metaconfig_yaml_file -i input_file -o output_folder
   ```

   calculate_power_singleqtl_cpma.py is for simulations that contain a single eqtl per simulation (no null snp or other eqtl is simulated in the same file). For cpma_type parameter, 0 is for cpma, 1 is for cpmax. 
   ```
   python calculate_power_singleqtl_cpma.py -c metaconfig_yaml_file -t cpma_type -f input_folder -i iterations
   ```
