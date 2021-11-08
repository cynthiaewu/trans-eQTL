import sys
import subprocess
import os
import pandas as pd
import argparse

def matrixeqtl_pipeline(tissue, factor, output, scripts_folder):
    if not os.path.isdir(output):
        os.mkdir(output)
    
    for i in range(1, 23):
        chr_folder = os.path.join(output, f'chr{i}')
        if not os.path.isdir(chr_folder):
            os.mkdir(chr_folder)

    for i in range(1, 23):
        eqtl_folder = os.path.join(f'{output}/chr{i}', 'matrixeqtl')
        if not os.path.isdir(eqtl_folder):
            os.mkdir(eqtl_folder)
    eqtl_file = f'matrixeqtl/matrixeQTL_results'
    
    factors = pd.read_csv(f'{output}/{tissue}_{factor}_peer_factors.tsv', sep='\t')
    residuals = pd.read_csv(f'{output}/{tissue}_{factor}_peer_residuals.tsv', sep='\t', nrows=1)
    factors.index = residuals.columns
    factors = factors.drop(columns='V1')
    factors.T.to_csv(f'{output}/{tissue}_{factor}_peer_factors_name.tsv', sep='\t')

    
    for i in range(1, 23):
        #matrix_cmd = f'Rscript /storage/polo/GTEx_v8_matrix_eQTL/Scripts/Run_Matrix_eQTL_PC_PEER_quantile_norm.r {input_folder}/chr{i}_SNP_intersect.tsv {output}/chr{i}/expression_PCs.csv {input_folder}/chr{i}_covariance_intersect.tsv {output}/chr{i}/{eqtl_file}_chr{i}'.split(' ')
        matrix_cmd = f'Rscript /storage/polo/GTEx_v8_matrix_eQTL/Scripts/Run_Matrix_eQTL_PC_PEER_quantile_norm.r {output}/chr{i}/chr{i}_SNP_intersect.tsv {output}/{tissue}_{factor}_peer_factors_name.tsv {output}/chr{i}/chr{i}_covariance_intersect.tsv {output}/chr{i}/{eqtl_file}_chr{i}_{factor}factors'.split(' ')
        subprocess.call(matrix_cmd)
        gzip_cmd = f'gzip {output}/chr{i}/{eqtl_file}_chr{i}_{factor}factors'.split(' ')
        subprocess.call(gzip_cmd)
    

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tissue', type=str, help='Tissue for gtex data')
    parser.add_argument('-f', '--factors', type=str, help='# peer factors')
    parser.add_argument('-o', '--output', type=str, help='Destination folder for output')
    parser.add_argument('-s', '--scripts_folder', type=str, help='Folder with scripts')
    input_values = parser.parse_args()

    matrixeqtl_pipeline(input_values.tissue, input_values.factors, input_values.output, input_values.scripts_folder)


if __name__ == '__main__':
    main()
