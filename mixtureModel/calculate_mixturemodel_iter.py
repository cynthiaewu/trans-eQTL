import sys
import subprocess
import os
import argparse

def mixture_pipeline(input_folder, scripts_folder):
    cpma_folder = os.path.join(input_folder, 'CPMA')
    mixture_folder = os.path.join(input_folder, 'mixtureModel')
    if not os.path.isdir(mixture_folder):
        os.mkdir(mixture_folder)
    pvalue_file = f'{cpma_folder}/gene-snp-eqtl_pvalue'
    mixture_file = f'{mixture_folder}/gene-snp-eqtl_mixturepvalue'

    mixture_cmd = (f'python {scripts_folder}/mixtureModel/calculate_mixturemodel.py -i {pvalue_file} -o {mixture_file}'.split(' '))
    subprocess.call(mixture_cmd)
    print(f'Finished calculating mixture, {input_folder}')
    
    


def iterate_folders(folder, scripts_folder, iterations):
    for i in range(iterations):
        input_folder = f'{folder}/Simulation_{i}'
        print(f'Starting mixture for Simulation {i}, {folder}')
        mixture_pipeline(input_folder, scripts_folder)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='Input folder with simulated expression and genotype files')
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-i', '--iterations', type=int, help='# of simulated folders to run cpma pipeline on')
    params = parser.parse_args()

    iterate_folders(folder=params.folder, scripts_folder=params.scripts_folder, iterations=params.iterations)


if __name__ == '__main__':
    main()
