import sys
import subprocess
import os
import argparse

def cpma_pipeline(input_folder, output, scripts_folder):
    if not os.path.isdir(output):
        os.mkdir(output)
    
    for i in range(1, 23):
        chr_folder = os.path.join(output, f'chr{i}')
        if not os.path.isdir(chr_folder):
            os.mkdir(chr_folder)

    for i in range(1, 23):
        PCs_folder = os.path.join(f'{output}/chr{i}', 'expressionPCs')
        if not os.path.isdir(PCs_folder):
            os.mkdir(PCs_folder)
    eqtl_file = f'expressionPCs/matrixeQTL_PCs'
    '''
    for i in range(1, 23):
        pca_cmd = (f'python {scripts_folder}/expressionPCs/get_expressionPCs_realdata.py -f {input_folder}/chr{i}_gene_intersect.tsv -o {output}/chr{i}/expression_PCs.csv'.split(' '))
        pca = subprocess.Popen(pca_cmd).wait()
    print('Finished computing PCA')
    '''
    for i in range(1, 23):
        #matrix_cmd = f'Rscript /storage/polo/GTEx_v8_matrix_eQTL/Scripts/Run_Matrix_eQTL_PC_PEER_quantile_norm.r {input_folder}/chr{i}_SNP_intersect.tsv {output}/chr{i}/expression_PCs.csv {input_folder}/chr{i}_covariance_intersect.tsv {output}/chr{i}/{eqtl_file}_chr{i}'.split(' ')
        matrix_cmd = f'Rscript /storage/polo/GTEx_v8_matrix_eQTL/Scripts/Run_Matrix_eQTL_PC_PEER_quantile_norm.r {input_folder}/chr{i}_SNP_intersect.tsv {output}/expressionPCs_noPEER_allchrgenes {input_folder}/chr{i}_covariance_intersect.tsv {output}/chr{i}/{eqtl_file}_chr{i}'.split(' ')
        subprocess.call(matrix_cmd)
    
 # Obtain zscores and pvalues in a snp by gene matrix format from matrix eQTL output
    pvalue_file = f'{eqtl_file}_pvalue_converted'
    zscore_file = f'{eqtl_file}_zscore'
    
    for i in range(1, 23):
        values_cmd = f'python {scripts_folder}/CPMA/get_values.py -i {output}/chr{i}/{eqtl_file}_chr{i} -z {output}/chr{i}/{zscore_file}_chr{i}'.split(' ')
        values = subprocess.Popen(values_cmd).wait()
        print(f'Finished getting pvalues and zscores for chr{i}')
    print('Finished getting pvalues and zscores for all chr')
     
    for i in range(1, 23):
        values_cmd = f'python {scripts_folder}/CPMA/tstat_to_pvalue.py -i {output}/chr{i}/{zscore_file}_chr{i} -p {output}/chr{i}/{pvalue_file}_chr{i}'.split(' ')
        values = subprocess.Popen(values_cmd).wait()
        print(f'Finished converting zscores to pvalues for chr{i}')
    print('Finished converting zscores to pvalues for all chr')
    
    #Calculate cpma values for all chr with matrix eqtl pvalues output
    cpma_file = f'{eqtl_file}_cpma_converted'
    
    for i in range(1, 23):
        cpma_cmd = f'python {scripts_folder}/CPMA/calculate_cpma_topx.py -i {output}/chr{i}/{pvalue_file}_chr{i} -x 1.0 -o {output}/chr{i}/{cpma_file}_chr{i}'.split(' ')
        cpma = subprocess.Popen(cpma_cmd).wait()
        print(f'Finished calculating cpma for chr{i}')
    print('Finished calculating cpma for all chr')
    
    
    for i in range(1, 23):
        chidist_cmd = (f'python {scripts_folder}/Simulator/compute_pvalue_chidist_realdata.py -i {output}/chr{i}/{cpma_file}_chr{i} -o {output}/chr{i}/{cpma_file}_pval_chr{i}'.split(' '))
        chidist = subprocess.Popen(chidist_cmd).wait()
    

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type=str, help='Input folder with expression files and matrix eqtl outputs')
    parser.add_argument('-o', '--output', type=str, help='Destination folder for output')
    parser.add_argument('-s', '--scripts_folder', type=str, help='Folder with scripts')
    input_values = parser.parse_args()

    cpma_pipeline(input_values.input, input_values.output, input_values.scripts_folder)


if __name__ == '__main__':
    main()
