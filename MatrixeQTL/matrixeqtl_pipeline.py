import sys
import subprocess
import os
import argparse

def matrixeqtl_pipeline(input_folder, scripts_folder, topx):
    genotype = f'{input_folder}/genotype.csv'
    expression = f'{input_folder}/expression.csv'
    cpma_folder = os.path.join(input_folder, 'CPMA')
    cpmax_folder = os.path.join(input_folder, 'CPMAx')
    if not os.path.isdir(cpma_folder):
        os.mkdir(cpma_folder)
    if not os.path.isdir(cpmax_folder):
        os.mkdir(cpmax_folder)
    eqtl_file = f'{cpma_folder}/gene-snp-eqtl'
    
    #Perform matrix eQTL to get gene-snp pairs
    #matrix_cmd = f'Rscript /storage/cynthiawu/trans_eQTL/Scripts/MatrixeQTL/gene-SNP_pairs.R -g {genotype} -e {expression} -o {eqtl_file}'.split(' ')
    matrix_cmd = f'Rscript {scripts_folder}/MatrixeQTL/gene-SNP_pairs.R -g {genotype} -e {expression} -o {eqtl_file}'.split(' ')
    subprocess.call(matrix_cmd)
    print(f'Finished matrix eQTL, {input_folder}')


def iterate_folders(folder, scripts_folder, topx, iterations):
    for i in range(iterations):
        input_folder = f'{folder}/Simulation_{i}'
        print(f'Starting CPMA for Simulation {i}, {folder}')
        matrixeqtl_pipeline(input_folder, scripts_folder, topx)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='Input folder with simulated expression and genotype files')
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-x', '--topx', default=0.1, type=float, help='Top x of genes to be used in cpma calculation (0 to 1)')
    parser.add_argument('-i', '--iterations', type=int, help='# of simulated folders to run cpma pipeline on')
    params = parser.parse_args()

    iterate_folders(folder=params.folder, scripts_folder=params.scripts_folder, topx=params.topx, iterations=params.iterations)


if __name__ == '__main__':
    main()
