import subprocess
import argparse


def mixture_pipeline(input_folder, scripts_folder):
    #targets = [ 0, 5, 10, 15, 20, 30, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500, 700, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]
    beta_values = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]
    #beta_values = [ 0.01, 0.02, 0.03, 0.04, 0.4]
    #beta_values = [0, 0.01, 0.1, 1]

    #targets = [0, 5, 10, 15, 20, 30, 40, 60, 80, 100, 200, 300, 400, 700, 1000, 5000, 10000, 15000]
    targets = [5, 10, 15, 30]
    #targets = [10, 100, 1000]
    #beta_values = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 1]
     
    print('Starting running mixture pipeline')

    
    for tar in targets:
        mixture_cmd = []
        for beta in beta_values:
            value = str(beta).replace(".","")
            mixture_cmd.append(f'python {scripts_folder}/mixtureModel/calculate_mixturemodel_iter.py -f {input_folder}/numTarget_{tar}/Beta_{value} -p {scripts_folder} -i 100'.split(' '))
        mixture_procs = [ subprocess.Popen(i) for i in mixture_cmd]
        for p in mixture_procs:
            p.wait()
    print('Finished calculating mixture') 
    
    
    print('Starting to calculate power')
    for tar in targets:
        power_cmd = []
        for beta in beta_values:
            value = str(beta).replace(".","")
            power_cmd.append(f'python {scripts_folder}/Simulator/calculate_power_singleqtl_cpma.py -c {input_folder}/numTarget_{tar}/Beta_{value}/metaconfig.yaml -m 5 -x 1.0 -f {input_folder}/numTarget_{tar}/Beta_{value} -i 100'.split(' '))
        power_procs = [ subprocess.Popen(i) for i in power_cmd]
        for p in power_procs:
            p.wait()
    print('Finished calculating power')
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_folder", required=True, help="Input folder with num_target and Beta folders")
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    params = parser.parse_args()

    mixture_pipeline(input_folder=params.input_folder,
          scripts_folder=params.scripts_folder)


if __name__ == "__main__":
    main()
