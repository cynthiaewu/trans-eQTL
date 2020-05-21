import pandas as pd
import argparse


def format_snps(input, output):
    # output of bcftools with desired filters
    #data = pd.read_csv('/storage/cynthiawu/trans_eQTL/GTex_filteredsnps.txt', sep='\t', names=['chr', 'pos', 'AF', 'AC', 'HWE'])
    data = pd.read_csv(input, sep='\t', names=['chr', 'pos', 'AF', 'AC', 'HWE'])
    data = data[data.chr != 'X']
    data['chr'] = 'chr' + data['chr'].astype(str)
    data['pos'] = data['pos'].astype(str)
    data['Position'] = data[['chr', 'pos']].agg('_'.join, axis=1)
    snps = list(data['Position'])
    #with open('/storage/cynthiawu/trans_eQTL/GTex_filteredsnps_pos.txt', 'w') as f:
    with open(output, 'w') as f:
        for sample in snps:
            f.write(sample + '\n') 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file of the saved output of bcftools with applied filtering parameters")
    parser.add_argument("-o", "--output", required=True, help="Ouptput list of postions of the snps that passed filtering")
    params = parser.parse_args()
    format_snps(params.input, params.output)


if __name__ == "__main__":
    main()
