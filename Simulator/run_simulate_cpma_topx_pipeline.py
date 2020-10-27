import subprocess
import argparse


def sim_cpmax_pipeline(input_folder, scripts_folder, topx, samplesize):
    #targets = [ 0, 20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500, 700, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]
    targets = [ 0, 5, 10, 15, 20, 30, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500, 700, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]
    #targets = [5, 10, 15, 30]
    #targets = [0, 20, 40]
    beta_values = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]
    #beta_values = [0, 0.1, 1]

    targets_str = f'{targets}'.replace(' ', '').replace('[', '').replace(']', '')
    beta_values_str = f'{beta_values}'.replace(' ', '').replace('[', '').replace(']', '')

    '''
    metaconfig_cmd = f'python {scripts_folder}/Simulator/write_metaconfig.py -i {input_folder} -t {targets_str} -b {beta_values_str} -s {samplesize}'.split(' ')
    subprocess.call(metaconfig_cmd)
    print('Finished writing metaconfig files')

    print('Starting generating config files')
    #configen_cmd = []
    for tar in targets:
        configen_cmd = []
        for beta in beta_values:
            value = str(beta).replace(".","")
            configen_cmd.append(f'python {scripts_folder}/Simulator/config_generator_specific.py -c {input_folder}/numTarget_{tar}/Beta_{value}/metaconfig.yaml -i 100 -o {input_folder}/numTarget_{tar}/Beta_{value}'.split(' '))
        configen_procs = [ subprocess.Popen(i) for i in configen_cmd]
        #print(configen_procs)
        for p in configen_procs:
            #print(p)
            p.wait()
    print('Finished generating config files')
    '''

    print('Starting simulations')
    for tar in targets:
        simulate_cmd = []
        for beta in beta_values:
            value = str(beta).replace(".","")
            simulate_cmd.append(f'python {scripts_folder}/Simulator/simulate_expression_givenoise.py -c {input_folder}/numTarget_{tar}/Beta_{value}/metaconfig.yaml -i 100 -o {input_folder}/numTarget_{tar}/Beta_{value}'.split(' '))
        simulate_procs = [ subprocess.Popen(i) for i in simulate_cmd]
        for p in simulate_procs:
            p.wait()
    print('Finished simulating files')
    
    print('Starting running cpma pipeline')
    
    
    for tar in targets:
        cpma_cmd = []
        for beta in beta_values:
            value = str(beta).replace(".","")
            cpma_cmd.append(f'python {scripts_folder}/CPMA/run_cpmax_pipeline_sim.py -f {input_folder}/numTarget_{tar}/Beta_{value} -p {scripts_folder} -x {topx} -i 100'.split(' '))
        cpma_procs = [ subprocess.Popen(i) for i in cpma_cmd]
        for p in cpma_procs:
            p.wait()
    print('Finished calculating cpma') 
     
   
    print('Starting comparing to chi distribution')
    for tar in targets:
        chidist_cmd = []
        for beta in beta_values:
            value = str(beta).replace(".","")
            chidist_cmd.append(f'python {scripts_folder}/Simulator/compute_pvalue_chidist.py -i {input_folder}/numTarget_{tar}/Beta_{value} -m 1 -x {topx} -t 1 -n 100'.split(' '))
        chidist_procs = [ subprocess.Popen(i) for i in chidist_cmd]
        for p in chidist_procs:
            p.wait()
    print('Finished comparing to chi distribution')
   

    print('Starting to calculate power')
    for tar in targets:
        power_cmd = []
        for beta in beta_values:
            value = str(beta).replace(".","")
            power_cmd.append(f'python {scripts_folder}/Simulator/calculate_power_singleqtl_cpma.py -c {input_folder}/numTarget_{tar}/Beta_{value}/metaconfig.yaml -m 1 -x {topx} -f {input_folder}/numTarget_{tar}/Beta_{value} -i 100'.split(' '))
        power_procs = [ subprocess.Popen(i) for i in power_cmd]
        for p in power_procs:
            p.wait()
    print('Finished calculating power')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_folder", required=True, help="Input folder with num_target and Beta folders")
    parser.add_argument('-p', '--scripts_folder', type=str, help='Input folder with scripts')
    parser.add_argument("-x", "--topx", default=0.1, type=float, help="Top x percent of genes to be used for cpma calculation")
    parser.add_argument("-s", "--samplesize", type=int, default=0, help="Sample size")
    params = parser.parse_args()

    sim_cpmax_pipeline(input_folder=params.input_folder,
          scripts_folder=params.scripts_folder,
          topx=params.topx,
          samplesize=params.samplesize)


if __name__ == "__main__":
    main()
