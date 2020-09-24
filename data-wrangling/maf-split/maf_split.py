#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def argParser():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--outputfile','-o',type=str,help='specify filename of your outputfile')
    parser.add_argument('--inputfile', '-i', type=str, help='specify filename of your inputfile')
    return parser.parse_args()

args = argParser()
outfile = args.outputfile
infile = args.inputfile


def clean_data(infile):
    data = pd.read_csv(infile,sep='\t')
    df = pd.DataFrame(data)

    rslt_df = df.loc[df['Project_Code'] == 'Ovary-AdenoCA']
    print(rslt_df)


def main():
    clean_data(infile)
    #calculate_entropy(new_dict)


if __name__ == '__main__':
    main()
