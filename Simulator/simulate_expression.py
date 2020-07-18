import numpy as np
import yaml
import argparse


def model(config, seed, output):
    with open(config) as f:
        params = yaml.load(f)
        print(params)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", required=True, help="Input config file with parameters")
    parser.add_argument("-s", "--seed", required=True, help="Seed for random generator")
    parser.add_argument("-o", "--output", required=True, help="Ouptput file with simulated expression matrix")
    params = parser.parse_args()
    model(params.config, params.seed, params.output)


if __name__ == "__main__":
    main()
