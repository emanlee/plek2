import numpy as np
import pandas as pd
from functions import filter_fasta, get_kmer, get_ORF, contact, prediction, init_data, output_acc
import argparse
import sys

def parse_arguments():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Instruction of PLEK2(more details in README.txt):')
    parser.add_argument('-i', '--input_file', help='Input file in fasta format', required=True)
    parser.add_argument('-m', '--model', help='ve: vertebrate, pl: plant', required=True)
    # parser.add_argument('-o', '--output_file', help='Output file', required=True)
    return parser.parse_args()

def main():
    args = parse_arguments()

    try:
        # Read fasta file
        file = filter_fasta(args.input_file)
    except FileNotFoundError:
        print("Error: File not found.")
        sys.exit(1)

    # Extract k-mer sequences
    get_kmer('kmer_seqs')

    # Extract ORF sequences
    get_ORF('seq_lines.fasta')

    # Get associations between k-mers and ORFs
    contact('kmer_6.txt', 'orf_length.txt')

    print("loading data ...")

    # Initialize data
    nums = init_data("features.txt")

    # Prepare data for prediction
    dataframe = pd.DataFrame(nums)
    dataset = dataframe.values

    X = np.zeros(shape=(len(dataset), 1, 5461, 1))

    for i, x in enumerate(dataset):
        x = np.asarray(x)
        x = np.resize(x, (1, 5461, 1))
        X[i] = x

    print("Loading Model and predicting ...")

    md = args.model

    # Perform prediction
    out = prediction(X, md)

    print(out)

    # Output prediction results
    prediction_result = output_acc(out)

    # out_file = open(args.output_file, 'w+')

if __name__ == "__main__":
    main()
