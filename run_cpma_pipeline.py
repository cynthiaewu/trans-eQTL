import sys
import subprocess
import os
import argparse

def cpma_pipeline(genotype, expression, tissue, covariate, european, output):
    genotype_file1 = genotype + '/GTExNormalizedSNPGenotypes_chr1.table.gz'
    intersect_file = f'{output}/intersect_samples.txt'
    expression_inter = f'{output}/Clean_expression_{tissue}_inter'
    covariates_inter = f'{output}/gtex650_{tissue}_intersect.pca'
    intersect_cmd = f'python find_sample_intersect.py -g {genotype_file1} -e {expression} -c {covariate} -i {intersect_file} -p {expression_inter} -q {covariates_inter}'.split(' ')
    #intersect = subprocess.Popen(intersect_cmd).wait()
    print('Intersecting Files Done')

    for i in range(1, 23):
        output_path = os.path.join(output, f'chr{i}')
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
    # subprocess.call(["mkdir", sample])

    genotype_inter = f'GTExNormalizedSNPGenotypes_inter'
    for i in range(1, 3):
        preprocess_cmd = f'python preprocess_genotypefile.py -g {genotype}/GTExNormalizedSNPGenotypes_chr{i}.table.gz -i {intersect_file} -o {output}/chr{i}/{genotype_inter}_chr{i}'.split(' ')
        #print(preprocess_cmd)
        #preprocess = subprocess.Popen(preprocess_cmd).wait()
        print(f'Finished preprocessing chr{i}')
    print('Preprocessing Genotype Files Done')

    genotype_filter = f'{genotype_inter}_filter'
    filtered_snp_list = '/storage/cynthiawu/trans_eQTL/GTex_filteredsnps_pos.txt'
    for i in range(1, 3):
        filtersnp_cmd = f'python filter_snps.py -i {output}/chr{i}/{genotype_inter}_chr{i} -s {filtered_snp_list} -o {output}/chr{i}/{genotype_filter}_chr{i} -c {i}'.split(' ')
        #filtersnp = subprocess.Popen(filtersnp_cmd).wait()
        print(f'Finished Filtering chr{i}')
    print('Filtering Snps Done')
    
    expression_filter = f'{expression_inter}_filter'
    filtergene_cmd = f'python filter_genes.py -i {expression_inter} -o {expression_filter}'.split(' ')
    #filtergene = subprocess.Popen(filtergene_cmd).wait()
    print('Filtering Genes Done')
    
    finalgenofile = genotype_filter
    finalexpfile = expression_filter
    if european == 1:
        eurosamp = '/storage/cynthiawu/trans_eQTL/European_samples.txt'
        finalgenofile = f'{genotype_filter}_european'
        for i in range(1, 3):
            eurogeno_cmd = f'python get_european_samples.py -i {output}/chr{i}/{genotype_filter}_chr{i} -e {eurosamp} -t 0 -o {output}/chr{i}/{finalgenofile}_chr{i}'.split(' ')
            print(eurogeno_cmd)
            eurogeno = subprocess.Popen(eurogeno_cmd).wait()
            print(f'Obtained European samples for chr{i}')
        finalexpfile = f'{expression_filter}_european'
        euroexp_cmd = f'python get_european_samples.py -i {expression_filter} -e {eurosamp} -t 1 -o {finalexpfile}'.split(' ')
        euroexp = subprocess.Popen(euroexp_cmd).wait()
        print('Obtained European samples for expression data')
        eurocovar_cmd = f'python get_european_samples_covariates.py -i {covariates_inter} -e {eurosamp} -o {output}/gtex650_eurosamp_covariates.pca'.split(' ')
        eurocovar = subprocess.Popen(eurocovar_cmd).wait()
        print('Obtained European samples for covariates')

    '''
    readcounts = os.path.join(output_path, sample + "_readcounts.txt")
    subprocess.call(["python", "prepare_readcounts.py", "-i", readcounts_path, "-o", readcounts])
    '''

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-g', '--genotype', type=str, help='Input folder with genotype files')
    parser.add_argument('-e', '--expression', type=str, help='Input expression file')
    parser.add_argument('-t', '--tissue', type=str, help='Tissue type')
    parser.add_argument('-c', '--covariate', type=str, help='Input covariate file')
    parser.add_argument('-s', '--european', type=int, help='0: Include all samples, 1: Include only European samples')
    parser.add_argument('-o', '--output', type=str, help='Destination folder for output')
    input_values = parser.parse_args()

    cpma_pipeline(input_values.genotype, input_values.expression, input_values.tissue, input_values.covariate, input_values.european, input_values.output)


if __name__ == '__main__':
    main()
