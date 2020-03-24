import pandas as pd
import argparse


def filter_snps(input, output):
    info_formatted = pd.read_csv('/storage/cynthiawu/trans_eQTL/GTEx_snp_AC_AF_HWE_formatted.txt', sep='\t')
    #snps = pd.read_csv('/storage/cynthiawu/trans_eQTL/Nerve-Tibial/chr2/GTExNormalizedSNPGenotypes_inter_coding_chr2.table', sep='\t')
    snps = pd.read_csv(input, sep='\t')
    snps_info = snps.join(info_formatted.set_index('Position'), on='chrom_start', how='left')
    snps_info = snps_info.dropna()
    keep = []
    for index, row in snps_info.iterrows():
        ac = False
        af = False
        for count in row['AC'].split(','):
            if int(count) >= 3:
                ac = True
        for freq in row['AF'].split(','):
            if float(freq) >= 0.01:
                af = True
        if ac & af & (row['HWE'] >= 0.01):
            keep.append(True)
        else:
            keep.append(False)
    filtered_snps_info = snps_info[keep]
    final_snps = filtered_snps_info.drop(columns = ['AC', 'AF', 'HWE'])
    final_snps.to_csv(output, index=False, header=True, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input genotype file")
    parser.add_argument("-o", "--output", required=True, help="Ouptput filtered genotype file, snps with AC >= 3, AF >= 0.01, HWE >= 0.01")
    params = parser.parse_args()
    filter_snps(params.input, params.output)


if __name__ == "__main__":
    main()
