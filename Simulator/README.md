# Simulator
1. Create a metaconfig yaml file with desired settings for simulations. User can copy over the example metaconfig.yaml and manually change parameters or they can use write_metaconfig.py script to generate a metaconfig.yaml file. The script can generate multiple metaconfig.yaml files for various numTarget genes and betas. If using the script, the input folder must contain the numTarget_x and Beta_y folders in the folder format numTarget_x/Beta_y with x=number of target genes and y=beta value in the format of the value without the decimal point (e.g. 0.05 would be 005).
   
#### metaconfig.yaml parameters:
- num_genes = # of total genes, default = 15,000 genes
- num_snps = # of total snps, default = 1 snp
- num_nullsnps = # of null snps with no beta effect size on all genes, default = 0 null snp
- num_factors = # of PEER factors, default = 0 factors
- sample_size = # samples/individuals, default = 500 samples
- allele_freq = allele frequency of snps, default = 0.5
- num_targets = # target genes of eqtl. If using script to generate metaconfig.yaml for multiple target genes, set this to be a list of target gene values separated by , with no spaces in between (e.g. 100,200,300). This can also be one value if only testing one # target genes 
- identity = True if no gene correlation, False if gene correlation (covariance matrix will be the gene correlationn matrix of Nerve Tibial GTeX data), default = True
- beta = 'value' or 'sd', default = "value" 
  - If set to "value," the parameter beta_value needs to be set to the desired beta effect size value which will be used for all target genes. 
  - If set to "sd," the beta effect size value for the target genes will be drawn from a normal distribution with mean 0 and SD set in the parameter beta_sd. 
- beta_sd = 'NA' if beta = 'value', else the standard deviation value of normal distribution with mean 0 to draw beta effect size values for the target genes.
- beta_value = 'NA' if beta = 'sd', else beta effect size for all target genes. If using script to generate metaconfig.yaml for multiple beta values, set this to be a list of beta values in float format separated by , with no spaces in between (e.g. 0.05,0.1,0.2). This can also be one value if only testing one # beta value 
- sig_threshold = significance threshold for calculating the power of the eqtl methods

To run the script to generate metaconfig.yaml files. 
   ```
   python write_metaconfig.py -x input_folder -t num_targets {-g num_genes} {-v num_snps} {-n num_nullsnps} {-f num_factors} {-s samplesize} {-a allele_freq} {-c identity} {-b beta} {-d beta_sd} {-v beta_value} {-h sig_threshold}
   ```
   
2. Generate config files with a given metaconfig yaml file

   Iterations are the number of desired simulations. A folder will be created for each iteration/simulation and will contain "beta.txt" which contains the beta effect sizes each snp has on each gene if beta_sd is chosen.
   ```
   python config_generator.py -c metaconfig_yaml_file -i iterations -o output_folder
   ```

3. Simulate genotype and expression files given that config generator has already been run to generate necessary config files. The metaconfig_yaml_file, iterations, and output_folder should be the same as (1.)
   ```
   python simulate_expression.py -c metaconfig_yaml_file -i iterations -o output_folder
   ```
4. Run eQTL methods with the simulated genotype and expression files.
5. Calculate power of the results from the eQTL methods.

   For Matrix eQTL method, calculate power across the Matrix eQTL results of simulations. Simulation folder contains all the Matrix eQTL results. num_snps_real is the number of snps/hypothesis tests (default = 10,000 snps) to correct for. An eqtl is successfully detected if at least one of the pvalues of its snp-gene pairs passes the significance threshold. If null snps were included in the simulations, use calculate_power_matrixeqtl_onesuccess_nullsnps.py, otherwise, use calculate_power_matrixeqtl_onesuccess.py
   ```
   python calculate_power_matrixeqtl_onesuccess.py -c metaconfig_yaml_file -f simulation_folder -n num_snps_real -i iterations 
   ```
   ```
   python calculate_power_matrixeqtl_onesuccess_nullsnps.py -c metaconfig_yaml_file -f simulation_folder -n num_snps_real -i iterations 
   ```
   
   For CPMA and xQTL: 
   Before calculating power, obtain the pvalue from each cpma value for each snp by comparing the cpma value to the chi square distribution. pvalues are adjusted for the number of tests performed per pair of #target genes and beta value parameter. For cpma_type parameter, 0 is for cpma, 1 is for cpmax. 
   ```
   python compute_pvalue_chidist.py -i input_file -c cpma_type -n num_test -o output_file
   ```

   calculate_power_singleqtl_cpma.py is for simulations that contain a single eqtl per simulation (no null snp or other eqtl is simulated in the same file). For cpma_type parameter, 0 is for cpma, 5 is for x-QTL. 
   ```
   python calculate_power_singleqtl_cpma.py -c metaconfig_yaml_file -t cpma_type -f input_folder -i iterations
   ```

# To run all steps above in one script (in progress)

For CPMA method, run_simulate_cpma_pipeline.py Input folder is where all the simulations and cpma files are to be written. Scripts folder is the path to all the trans eqtl scripts. Sample size is the desired sample size to be used in the simulations. 
 ```
 python run_simulate_cpma_pipeline.py -i input_folder -p scripts_folder -s samplesize
 ```
 For CPMA_x method, 
 ```
 run_simulate_cpma_topx_pipeline.py -i input_folder -p scripts_folder -s samplesize -x topx
 ```   

