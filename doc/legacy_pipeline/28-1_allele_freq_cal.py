#!/usr/bin/python

import sys
import os
import glob
import argparse

import matplotlib
import gzip
matplotlib.use('Agg')
from numpy.random import randn
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def to_percent(y, position):
# Ignore the passed in position. This has the effect of scaling the default
# tick locations.
    s = str(10 * y)

# The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'



def histogram_plot(list_allele_frq, output_png_file):
# Make a normed histogram. It'll be multiplied by 100 later.
    plt.hist(list_allele_frq, bins=50, normed=True)
    plt.xlim(0,1)
    formatter = FuncFormatter(to_percent)
    plt.gca().yaxis.set_major_formatter(formatter)

#    plt.show()
#    save(output_png_file, ext="png", close=False, verbose=True)
    plt.savefig(output_png_file, bbox_inches='tight', format='png')
    return


def run_program(input_vcf, output_prefix, sample):

    list_allele_frq = []

    output_file = '%s.xls' %output_prefix
    output_png_file = '%s.png' %output_prefix

    caller_type = "gatk"

    if input_vcf.endswith(".vcf.gz"):
        input_fp = gzip.open(input_vcf, 'r')
    else:
        input_fp = open(input_vcf, 'r')
    output_fp = open(output_file, 'w')
    print >> output_fp, 'chr_id\tbp\trs\tref\talt\tdp\tref_dp\talt_dp\tratio'

    for line in input_fp:
        if line.startswith("##INFO=<ID=DP4"):
            caller_type = 'samtools'
        elif line.startswith("##INFO=<ID=AO,"):
            caller_type = 'proton'
        elif line.startswith("#CHROM"):
            print "#CHROM"
            header_index = {}
            units = line.strip().split()
            for unit in units:
                header_index[unit] = units.index(unit)
        elif line.startswith("#"):
            continue
        else:
            units = line.strip().split()
            chr_id = units[0]
            bp = units[1]
            rs = units[2]
            ref = units[3]
            alt = units[4]
            info = units[7]
            list_info = info.split(';')
            genotype_info = units[header_index[sample]]
            gt = genotype_info.split(':')[0]

            if ',' in alt:
                continue

            if caller_type == 'samtools':
                for item in list_info:
                    if item.startswith('DP='):
                        dp = int(item.split('=')[1])
                    if item.startswith('DP4='):
                        foward_ref  = int(item.split('=')[1].split(',')[0])
                        reverse_ref = int(item.split('=')[1].split(',')[1])
                        foward_alt  = int(item.split('=')[1].split(',')[2])
                        reverse_alt = int(item.split('=')[1].split(',')[3])
                        ref_dp = foward_ref + reverse_ref
                        alt_dp = foward_alt + reverse_alt
                        dp = ref_dp + alt_dp
                        if dp == 0:
                            continue
                        ratio = round(float(alt_dp) / (ref_dp + alt_dp), 3)

            elif caller_type == 'proton':
                for item in list_info:
                    if item.startswith('AF='):
                        ratio = float(item.split('=')[1])
                    if item.startswith('DP='):
                        dp = int(item.split('=')[1])
                    if item.startswith('AO='):
                        alt_dp = int(item.split('=')[1])

                ref_dp = dp - alt_dp

            else:
                gt = genotype_info.split(':')[0]
                ad = genotype_info.split(':')[1]
                ref_dp = int(ad.split(',')[0])
                alt_dp = int(ad.split(',')[1])
                dp = ref_dp + alt_dp
                if dp == 0:
                    continue
                ratio = round(float(alt_dp) / (ref_dp + alt_dp), 3)

            if len(ref) == 1 and len(alt) == 1 and dp >= 10 and gt == "0/1":
                list_allele_frq.append(ratio)

                list_output = []
                list_output.append(chr_id)
                list_output.append(bp)
                list_output.append(rs)
                list_output.append(ref)
                list_output.append(alt)
                list_output.append(str(dp))
                list_output.append(str(ref_dp))
                list_output.append(str(alt_dp))
                list_output.append(str(ratio))
                print >> output_fp, '\t'.join(list_output)

    histogram_plot(list_allele_frq, output_png_file)
    input_fp.close()
    output_fp.close()
    return

def usage():
    message='''
python %s --input "*.vcf"

-i, --input  : input vcf file
-o, --output : output prefix
-s, --sample : sample name

##optional
    ''' %sys.argv[0]
    print message

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-s', '--sample')
    args = parser.parse_args()
    try:
        len(args.input) > 0

    except:
        usage()
        sys.exit(2)

    run_program(args.input, args.output, args.sample)

if __name__ == '__main__':
    main()

