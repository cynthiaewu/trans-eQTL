import subprocess
import argparse


def run_matrixeqtl_pipeline(input_folder, scripts_folder, topx, samplesize):
    #targets = [ 0, 20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500, 700, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]
    targets = [ 1, 5, 10, 15, 20, 30, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500, 700, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]
    #targets = [5, 10, 15, 30]
    #targets = [0, 5, 10, 15, 20, 30, 40, 60, 80, 100, 200, 300, 400, 700, 1000, 5000, 10000, 15000]
    #targets = [1, 10, 100, 1000]
    #targets = [10000]
    #beta_values = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 1]
    beta_values = [0, 0.01, 0.05, 0.1, 1]
    #beta_values = [0.01]
    #beta_values = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]
    #beta_values = [0, 0.1, 1]

    targets_str = f'{targets}'.replace(' ', '').replace('[', '').replace(']', '')
    beta_values_str = f'{beta_values}'.replace(' ', '').replace('[', '').replace(']', '')

    print('Starting running matrixeqtl')    
    
    for tar in targets:
        cpma_cmd = []
        for beta in beta_values:
            value = str(beta).replace(".","")
            cpma_cmd.append(f'python {scripts_folder}/MatrixeQTL/matrixeqtl_pipeline.py -f {input_folder}/numTarget_{tar}/Beta_{value} -p {scripts_folder} -x {topx} -i 100'.split(' '))
        cpma_procs = [ subprocess.Popen(i) for i in cpma_cmd]
        for p in cpma_procs:
            p.wait()
    print('Finished matrixeqtl') 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_folder", required=True, help="Input folder with num_target and Beta folders")
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument("-x", "--topx", default=0.1, type=float, help="Top x percent of genes to be used for cpma calculation")
    parser.add_argument("-s", "--samplesize", type=int, default=0, help="Sample size")
    params = parser.parse_args()

    run_matrixeqtl_pipeline(input_folder=params.input_folder,
          scripts_folder=params.scripts_folder,
          topx=params.topx,
          samplesize=params.samplesize)


if __name__ == "__main__":
    main()
