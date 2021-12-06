import sys
import subprocess
import os
import argparse

def pipeline(output, scripts_folder, tissue, factor):
    preprocess_cmd = f'python {scripts_folder}/preprocess/preprocess.py -t {tissue} -o {output}'.split(' ')
    #subprocess.call(preprocess_cmd)
    print(f'Finished preprocess, {tissue}')

    peer_cmd = f'Rscript {scripts_folder}/preprocess/runPeer.r {output} {tissue} {factor}'.split(' ')
    subprocess.call(peer_cmd)
    print(f'Finished PEER, {tissue}')



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=str, help='Output folder')
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument('-t', '--tissue', type=str, help='Tissue of expression data')
    parser.add_argument('-f', '--factor', type=int, help='# of PEER factors')
    params = parser.parse_args()

    pipeline(output=params.output, scripts_folder=params.scripts_folder, tissue=params.tissue, factor=params.factor)


if __name__ == '__main__':
    main()
