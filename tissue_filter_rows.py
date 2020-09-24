#!/usr/bin/env python
import requests
import json
import argparse
import linecache
import pandas as pd
import io
from itertools import islice

def argParser():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--tumortype', '-t', type=str, help='specify tumor tissue type')
    parser.add_argument('--vcffile','-v',type=str,help='specify filename of vcf')
    parser.add_argument('--tsvfile', '-f', type=str, help='specify filename of tsv file with tumor types mapped')
    parser.add_argument('--pathway', '-p', type=str, help='specify pathway of genes in all lowercase, either hrr or mmr')

    return parser.parse_args()

def read_vcf(vcffile):
     
    with open(vcffile, "r") as myfile:
        head = list(islice(myfile, 253))

    with open('fileName.csv', "w") as f2:
        for item in head:
            f2.write(item)

    return vcffile

def read_tsv(tsvfile, tumortype):
    tcga_ids_bytissue = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    mapperfile = pd.read_csv(tsvfile,sep='\t')
    for i in range(0,len(mapperfile.index)):
        tissue = mapperfile.iloc[i][1]
        if tissue == tumortype:
            tcga_ids_bytissue.append(mapperfile.iloc[i][3])
    return tcga_ids_bytissue

def hrr(tcga_ids_bytissue,vcffile,tumortype):
    with open(vcffile,'r') as f:
        for l in f:
            if not l.startswith(('##','#')):
                line_list = l.split()
                if line_list[0] != 'X' or 'Y':
                    chrom = int(line_list[0])
                    pos = int(line_list[1])
                    print(chrom)

#SSBP1
                if chrom == 7:
                    if 141438121 <= pos <= 141450288:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#RAD50
                if chrom == 5:
                    if 131892616 <= pos <= 131980313:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#MRE11
                if chrom == 11:
                    if 94150466 <= pos <= 94227040:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#NBS1 (NBN, NBS)
                if chrom == 8:
                    if 90945564 <= pos <= 90996952:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#ATM
                if chrom == 11:
                    if 108093559 <= pos <= 108239829:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#CHEK2
                if chrom == 22:
                    if 29083731 <= pos <= 29137822:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#TP53
                if chrom == 17:
                    if 7571720 <= pos <= 7590868:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#BARD1
                if chrom == 2:
                    if 215590370 <= pos <= 215674428:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#CtIP / RBBP8

                if chrom == 18:
                    if 20513295 <= pos <= 20606451:
                        with open('fileName.csv','a') as o:
                            o.write(l)

#BRIP1
                if chrom == 17:
                    if 59756547 <= pos <= 59940920:
                        with open('fileName.csv','a') as o:
                            o.write(l)

#BRCA1
                if chrom == 17:
                    if 41196312 <= pos <= 41277500:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#PALB2
                if chrom == 16:
                    if 23614481 <= pos <= 23652678:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#BRCA2
                if chrom == 13:
                    if 32889617 <= pos <= 32973809:
                        with open('fileName.csv','a') as o:
                            o.write(l)
#RAD51
                if chrom == 13:
                    if 40987327 <= pos <= 41024356:
                        with open('fileName.csv','a') as o:
                            o.write(l)

def mmr(tcga_ids_bytissue, vcf_dataframe,tumortype):
    pos_list = []
    select_columns = vcf_dataframe[tcga_ids_bytissue]
    for i in range(0, len(select_columns.index)):
        pos = select_columns.iloc[i][1]
        if 37034841 <= pos <= 37092337:
            pos_list.append(pos)
        if 47630206 <= pos <= 47739716:
            pos_list.append(pos)
        if 48010221 <= pos <= 48034092:
            pos_list.append(pos)
        if 6012870 <= pos <= 6048737:
            pos_list.append(pos)
    select_rows = select_columns.loc[select_columns['POS'].isin(pos_list)]

    select_rows.to_csv("%s.vcf" % tumortype,sep='\t',header=True, index=False)



def main():
    args = argParser()
    tumortype = args.tumortype
    vcffile = args.vcffile
    tsvfile = args.tsvfile
    pathway = args.pathway
    vcf_dataframe = read_vcf(vcffile)
    tcga_ids_bytissue = read_tsv(tsvfile,tumortype)
    if pathway == "hrr":
        hrr(tcga_ids_bytissue,vcf_dataframe,tumortype)
    if pathway == "mmr":
        mmr(tcga_ids_bytissue,vcf_dataframe,tumortype)



if __name__ == '__main__':
    main()
