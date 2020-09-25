#!/usr/bin/env python

import argparse
import pandas as pd
from itertools import islice

def argParser():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--tumortype', '-t', type=str, help='specify tumor tissue type')
    parser.add_argument('--vcffile','-v',type=str,help='specify filename of vcf')
    parser.add_argument('--tsvfile', '-f', type=str, help='specify filename of tsv file with tumor types mapped')

    return parser.parse_args()

def read_vcf(vcffile):

    with open(vcffile, "r") as myfile:
        head = list(islice(myfile, 255))

    with open('testing.vcf', "w") as f2:
        for item in head:
            f2.write(item)


def read_mapper(tsvfile, tumortype):
    tcga_ids_bytissue = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    #tcga_ids_bytissue = []
    mapperfile = pd.read_csv(tsvfile,sep='\t')
    for i in range(0,len(mapperfile.index)):
        tissue = mapperfile.iloc[i][1]
        if tissue == tumortype:
            tcga_ids_bytissue.append(mapperfile.iloc[i][3])
    print(tcga_ids_bytissue)
    with open('tissue_ids.txt', 'w') as f:
    
        for i in tcga_ids_bytissue:
            f.write(i+'\n')
    return tcga_ids_bytissue

def get_columns(vcffile, tcga_ids_bytissue):
    columns = [0,1,2,3,4,5,6,7,8]
    with open(vcffile,'r') as f:
        for l in f:
            if l.startswith('#CHROM'):
                header_list = l.split()
                break
    for i,item in enumerate(header_list):
        if item in tcga_ids_bytissue:
            columns.append(i)
    print(columns)
    return columns
def write_header(columns,vcffile):
    header_final = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    with open(vcffile,'r') as f:
        header = islice(f, 255, 256)
        for item in header:
            if item in columns:
                header_final.append(item)    
    print(header_final)
    return header_final         

def filter_columns(columns,vcffile):     
    with open('testing-ovary.vcf','a') as f:
        with open(vcffile,'r') as v:
            for l in v:
                if not l.startswith(('##','#')):
                    line_list = l.split()
                    line_to_add = []
                    for i,item in enumerate(line_list):
                        if i in columns:
                            line_to_add.append(item)
                        
                    f.write(('\t'.join(line_to_add) + '\n'))

def main():
    args = argParser()
    tumortype = args.tumortype
    vcffile = args.vcffile
    tsvfile = args.tsvfile
    read_vcf(vcffile)
    tcga_ids_bytissue = read_mapper(tsvfile,tumortype)
    
    columns = get_columns(vcffile,tcga_ids_bytissue)
    write_header(columns,vcffile)
    vcffile = filter_columns(columns,vcffile)

if __name__ == '__main__':
    main()
