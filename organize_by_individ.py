#!/usr/bin/env python
import argparse
from itertools import islice

def argParser():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--vcffile','-v',type=str,help='specify filename of vcf')

    return parser.parse_args()

def write_header(vcffile):
    with open(vcffile, "r") as myfile:
        head = list(islice(myfile, 257))

    with open('individuals.tsv', "w") as f2:
        for item in head:
            f2.write(item)
    with open('individuals.tsv','a') as f3:
        f3.write("ID	CHROM	POS	REF	ALT	INFO"+'\n')

def reorganize(vcffile):
    with open(vcffile,'r') as f:
        for l in f:
            if l.startswith('#CHROM'):
                header_list = l.split()
                print(header_list)
                break

        for l in f:
            if not l.startswith('#'):
                line_list = l.split()
                chrom = int(line_list[0])
                pos = int(line_list[1])
                ref = line_list[3]
                alt = line_list[4]
                for i,genotype in enumerate(line_list):
                    reorg_linelist = []
                    if '.:.' not in genotype: #!= './.:.:.:.:.:.:.:.:.:.:.:.:.:.' or './.:.:.:.:.:.:.:.:.:.:.:.:.:.:.':
                        if header_list[i].startswith('TCGA'):
                            tcga_id = header_list[i]
                            print(tcga_id,i)
                            reorg_linelist=[tcga_id,str(chrom),str(pos),ref,alt,genotype]
        #                print(reorg_linelist)
                            with open("individuals.tsv",'a') as indiv_file:
                                indiv_file.write('\t'.join(reorg_linelist) + '\n')

                

def main():
    args = argParser()
    vcffile = args.vcffile
    write_header(vcffile)   
    reorganize(vcffile)

if __name__ == '__main__':
    main()
