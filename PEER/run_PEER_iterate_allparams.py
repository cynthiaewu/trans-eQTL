import sys
import subprocess
import argparse
import numpy as np

def run_peer(input_folder, scripts_folder, param_list):
    parameters = np.loadtxt(param_list, dtype=str)
    peer_cmd = []
    print(parameters)
    core = 1
    for param in parameters:
        peer_cmd.append(f'taskset -c {core} python {scripts_folder}/PEER/run_PEER_iterate.py -f {input_folder}/{param} -p {scripts_folder} -i 100'.split(' '))
        core += 1
    #proc = [subprocess.Popen(i, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True) for i in peer_cmd]
    proc = [subprocess.Popen(i) for i in peer_cmd]
    for p in proc:
        '''
        cur_file = f'{input_folder}/{param}/logfile'
        with open(cur_file, 'w') as logfile:
            for line in p.stdout:
                sys.stdout.write(line)
                logfile.write(line)
                logfile.flush()
            print(p)
        '''
        p.wait()
    print('Finished calculating peer')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='Input folder with simulated expression and genotype files')
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-l', '--param_list', type=str, help='Input list of all parameters to run peer pipeline')
    params = parser.parse_args()

    run_peer(input_folder=params.folder, scripts_folder=params.scripts_folder, param_list=params.param_list)


if __name__ == '__main__':
    main()
