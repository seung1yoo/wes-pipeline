#!/usr/bin/python

import sys
import os
import glob
import gzip
import argparse


def run_program(input_file, filt_depth, wes_config_file):
    input_vcf_list = glob.glob(input_file)
    Dic_Sample_snv = {}
    for input_vcf in input_vcf_list:
        Sample_ID = os.path.basename(input_vcf).split(".")[0]
        Dic_Sample_snv[Sample_ID] = {}
        with gzip.open(input_vcf, 'r') as vcf_reader:
            for line in vcf_reader:
                if line.startswith("#"):
                    continue
                else:
                    units = line.strip().split()
                    if len(units) != 10:
                        continue
                    chr_ID = units[0]
                    bp = units[1]
                    ref = units[3]
                    alt = units[4]
                    info = units[7].split(';')
                    for item in info:
                        if item.startswith('DP='):
                            dp = item.split('DP=')[1]
                    if len(ref) == 1 and len(alt) == 1:
                        gt = units[9].split(':')[0]
                        if int(dp) >= int(filt_depth):
                            key = '%s_%s' % (chr_ID, bp)
                            Dic_Sample_snv[Sample_ID][key] = gt

    if wes_config_file:
        sampleID_deliveryID_pairInfo = {}
        with file(wes_config_file) as wes_config_reader:
            for line in wes_config_reader:
                if line.startswith("delivery_tbi_id"):
                    items = line.split("=")[1].strip()
                    for pairInfo in items.split(","):
                        deliveryID, sampleID = pairInfo.split(":")[:2]
                        sampleID_deliveryID_pairInfo[sampleID] = deliveryID

    sample_list = Dic_Sample_snv.keys()
    sample_list.sort()
    set_pair = set()

    print "\t".join(["Sample ID", "Sample ID", "Total SNVs", "Same SNVs", "Diff SNVs", "Concordance ratio"])
    for Idv_1_key in sample_list:
        for Idv_2_key in sample_list:
            if Idv_1_key != Idv_2_key:
                pair_1_id = '%s:%s' % (Idv_1_key, Idv_2_key)
                pair_2_id = '%s:%s' % (Idv_2_key, Idv_1_key)
                if pair_1_id in set_pair:
                    continue
                elif pair_2_id in set_pair:
                    continue
                else:
                    total_gt = 0
                    exact_gt = 0
                    nonexact_gt = 0
                    for bp_key in Dic_Sample_snv[Idv_2_key].keys():
                        if Dic_Sample_snv[Idv_1_key].has_key(bp_key):
                            Sample_1_gt = Dic_Sample_snv[Idv_1_key][bp_key]
                            Sample_2_gt = Dic_Sample_snv[Idv_2_key][bp_key]
                            total_gt += 1
                            if Sample_1_gt == Sample_2_gt:
                                exact_gt += 1
                            else:
                                nonexact_gt += 1
                    set_pair.add(pair_1_id)
                    set_pair.add(pair_2_id)
                    pct = round(exact_gt / float(total_gt) * 100, 2)
                    if wes_config_file:
                        print '%s\t%s\t%s\t%s\t%s\t%s' % (sampleID_deliveryID_pairInfo[Idv_1_key], sampleID_deliveryID_pairInfo[Idv_2_key], total_gt, exact_gt, nonexact_gt, pct)
                    else:
                        print '%s\t%s\t%s\t%s\t%s\t%s' % (Idv_1_key, Idv_2_key, total_gt, exact_gt, nonexact_gt, pct)

    return

def usage():
    message = '''
python %s --input "*.vcf"

-i, --input: input vcf file ex) "*.vcf"

##optional
    ''' % sys.argv[0]
    print message


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-d', '--depth', default=10)
    parser.add_argument('-c', '--wesconfig', default="")

    #    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    try:
        len(args.input) > 0

    except:
        usage()
        sys.exit(2)

    run_program(args.input, args.depth, args.wesconfig)


#    try:
#        run_program(args.input)
#    except:
#        print 'ERROR'

if __name__ == '__main__':
    main()

