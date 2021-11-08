import sys
import subprocess
import os
import argparse

def matrixeqtl_pipeline(input_folder, output, scripts_folder):
    if not os.path.isdir(output):
        os.mkdir(output)
    '''
    for i in range(1, 23):
        chr_folder = os.path.join(output, f'chr{i}')
        if not os.path.isdir(chr_folder):
            os.mkdir(chr_folder)

    for i in range(1, 23):
        eqtl_folder = os.path.join(f'{output}/chr{i}', 'matrixeqtl')
        if not os.path.isdir(eqtl_folder):
            os.mkdir(eqtl_folder)
    eqtl_file = f'matrixeqtl/matrixeQTL_results'
    '''
    '''
    for i in range(1, 23):
        pca_cmd = (f'python {scripts_folder}/expressionPCs/get_expressionPCs_realdata.py -f {input_folder}/chr{i}_gene_intersect.tsv -o {output}/chr{i}/expression_PCs.csv'.split(' '))
        pca = subprocess.Popen(pca_cmd).wait()
    print('Finished computing PCA')
    '''
    for i in range(1, 23):
        #matrix_cmd = f'Rscript /storage/polo/GTEx_v8_matrix_eQTL/Scripts/Run_Matrix_eQTL_PC_PEER_quantile_norm.r {input_folder}/chr{i}_SNP_intersect.tsv {output}/chr{i}/expression_PCs.csv {input_folder}/chr{i}_covariance_intersect.tsv {output}/chr{i}/{eqtl_file}_chr{i}'.split(' ')
        #matrix_cmd = f'Rscript /storage/polo/GTEx_v8_matrix_eQTL/Scripts/Run_Matrix_eQTL_PC_PEER_quantile_norm.r {output}/chr{i}/chr{i}_SNP_intersect.tsv {output}/Nerve-Tibial_peer_residuals.tsv {output}/chr{i}/chr{i}_covariance_intersect.tsv {output}/chr{i}/{eqtl_file}_chr{i}'.split(' ')
        matrix_cmd = f'Rscript /storage/cynthiawu/trans_eQTL/Scripts/MatrixeQTL/Run_Matrix_eQTL_PC_PEER_quantile_norm.r {output}/chr{i}/chr{i}_SNP_intersect.tsv {output}/Thyroid_25_peer_residuals_format.tsv {output}/chr{i}/chr{i}_covariance_intersect.tsv {output}/chr{i}/matrixeqtl_results_ciseqtls_chr{i}'.split(' ')
        subprocess.call(matrix_cmd)
        #gzip_cmd = f'gzip {output}/chr{i}/{eqtl_file}_chr{i}'.split(' ')
        gzip_cmd = f'gzip {output}/chr{i}/matrixeqtl_results_ciseqtls_chr{i}'.split(' ')
        print(gzip_cmd)
        subprocess.call(gzip_cmd)
    

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type=str, help='Input folder with expression files and matrix eqtl outputs')
    parser.add_argument('-o', '--output', type=str, help='Destination folder for output')
    parser.add_argument('-s', '--scripts_folder', type=str, help='Folder with scripts')
    input_values = parser.parse_args()

    matrixeqtl_pipeline(input_values.input, input_values.output, input_values.scripts_folder)


if __name__ == '__main__':
    main()
