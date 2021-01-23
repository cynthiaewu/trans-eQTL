import sys
import subprocess
import argparse
import numpy as np

def run_alltissues(tissue_list):
    #tissues = np.loadtxt('/storage/cynthiawu/trans_eQTL/GTex.v8/tissue_list.txt', dtype=str)
    tissues = np.loadtxt(tissue_list, dtype=str)
    for tissue in tissues:        
        print(f'Started cpma pipeline for {tissue}')
        with open(f'/storage/cynthiawu/trans_eQTL/GTex.v8/log/{tissue}_2.log', 'w') as logfile:
            pipeline_cmd = f'taskset -c 1-21 python -u run_cpma_pipeline_realdata.py -i /storage/polo/GTEx_v8_matrix_eQTL/GTEx_v8_Euro_only_PEER/{tissue} -o /storage/cynthiawu/trans_eQTL/GTex.v8/{tissue}'.split(' ')
            proc = subprocess.Popen(pipeline_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
            for line in proc.stdout:
                sys.stdout.write(line)
                logfile.write(line)
                logfile.flush()
            proc.wait()
            print(f'Finished cpma pipeline for {tissue}')


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tissue_list', type=str, help='Input list of all tissues to run real data cpma pipeline')
    input_values = parser.parse_args()

    run_alltissues(input_values.tissue_list)


if __name__ == '__main__':
    main()
