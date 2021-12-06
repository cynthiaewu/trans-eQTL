import sys
import subprocess
import os
import argparse

def gamma_pipeline(input_folder, scripts_folder):
    cpma_folder = os.path.join(input_folder, 'CPMA')
    gamma_folder = os.path.join(input_folder, 'gammaModel')
    if not os.path.isdir(gamma_folder):
        os.mkdir(gamma_folder)
    pvalue_file = f'{cpma_folder}/gene-snp-eqtl_pvalue'
    gamma_file = f'{gamma_folder}/gene-snp-eqtl_gammapvalue'

    gamma_cmd = (f'python {scripts_folder}/gammaModel/calculate_gammaModel.py -i {pvalue_file} -o {gamma_file}'.split(' '))
    subprocess.call(gamma_cmd)
    print(f'Finished calculating gamma, {input_folder}')
    
    


def iterate_folders(folder, scripts_folder, iterations):
    for i in range(iterations):
        input_folder = f'{folder}/Simulation_{i}'
        print(f'Starting gamma for Simulation {i}, {folder}')
        gamma_pipeline(input_folder, scripts_folder)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='Input folder with simulated expression and genotype files')
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-i', '--iterations', type=int, help='# of simulated folders to run cpma pipeline on')
    params = parser.parse_args()

    iterate_folders(folder=params.folder, scripts_folder=params.scripts_folder, iterations=params.iterations)


if __name__ == '__main__':
    main()
