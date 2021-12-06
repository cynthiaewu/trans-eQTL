import pandas as pd
import argparse


def normalize(input, output):
    df = pd.read_csv(input, sep='\t', index_col=0)
    normalized_df=(df-df.mean())/df.std()
    normalized_df.to_csv(output, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Matrix eqtl results output file as input file")
    parser.add_argument("-o", "--output", required=True, help="Output folder with simulated files")
    params = parser.parse_args()

    normalize(input=params.input,
          output=params.output)


if __name__ == "__main__":
    main()
