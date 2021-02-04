import sys
import subprocess
import os
import argparse

def peer_pipeline(input_folder, scripts_folder):
    peer_cmd = f'Rscript {scripts_folder}/PEER/Run_PEER.r {input_folder}'.split(' ')
    subprocess.call(peer_cmd)
    print(f'Finished PEER, {input_folder}')


def iterate_folders(folder, scripts_folder, iterations):
    for i in range(iterations):
        input_folder = f'{folder}/Simulation_{i}'
        print(f'Starting PEER for Simulation {i}, {folder}')
        peer_pipeline(input_folder, scripts_folder)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='Input folder with simulated expression and genotype files')
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-i', '--iterations', type=int, help='# of simulated folders to run cpma pipeline on')
    params = parser.parse_args()

    iterate_folders(folder=params.folder, scripts_folder=params.scripts_folder, iterations=params.iterations)


if __name__ == '__main__':
    main()
