#!/usr/bin/python

import sys
import glob
import os
import locale
import time
import argparse


def MakeDirectory(Dir):
    if not os.path.exists(Dir):
        os.system('mkdir -p %s' %Dir)

def fastq_statistics_run(first_in_fastq, second_in_fastq, first_filt_fastq, second_filt_fastq, sample, output_file):

    fiFirst = os.popen("zcat %s" % first_in_fastq)
    fiSecond = os.popen("zcat %s" % second_in_fastq)
    if first_filt_fastq[-3:] == '.gz':
        filtFirst = os.popen("zcat %s" % first_filt_fastq)
    else:
        filtFirst = os.popen("cat %s" % first_filt_fastq)
    if second_filt_fastq[-3:] == '.gz':
        filtSecond = os.popen("zcat %s" % second_filt_fastq)
    else:
        filtSecond = os.popen("cat %s" % second_filt_fastq)

    fiFirst_ReadCnt = 0
    fiSecond_ReadCnt = 0
    fiFirst_Base = 0
    fiSecond_Base = 0

    filt_First_ReadCnt = 0
    filt_Second_ReadCnt = 0
    filt_First_Base = 0
    filt_Second_Base = 0

    while True:
        line1_id = fiFirst.readline().rstrip()
        line1_seq = fiFirst.readline().rstrip()
        line1_strand = fiFirst.readline().rstrip()
        line1_quality = fiFirst.readline().rstrip()
        line2_id = fiSecond.readline().rstrip()
        line2_seq = fiSecond.readline().rstrip()
        line2_strand = fiSecond.readline().rstrip()
        line2_quality = fiSecond.readline().rstrip()
        
        if not line1_id or not line2_id:
            break

        fiFirst_ReadCnt += 1
        fiSecond_ReadCnt += 1
        fiFirst_Base += len(line1_seq)
        fiSecond_Base += len(line2_seq)

    while True:
        line1_id = filtFirst.readline().rstrip()
        line1_seq = filtFirst.readline().rstrip()
        line1_strand = filtFirst.readline().rstrip()
        line1_quality = filtFirst.readline().rstrip()
        line2_id = filtSecond.readline().rstrip()
        line2_seq = filtSecond.readline().rstrip()
        line2_strand = filtSecond.readline().rstrip()
        line2_quality = filtSecond.readline().rstrip()
        
        if not line1_id or not line2_id:
            break

        filt_First_ReadCnt += 1
        filt_Second_ReadCnt += 1
        filt_First_Base += len(line1_seq)
        filt_Second_Base += len(line2_seq)


    Total_read_cnt = fiFirst_ReadCnt + fiSecond_ReadCnt
    Total_base_cnt = fiFirst_Base + fiSecond_Base
    Total_filt_read_cnt = filt_First_ReadCnt + filt_Second_ReadCnt
    Total_filt_base_cnt = filt_First_Base + filt_Second_Base

    Total_filt_read_pct = round(float(Total_filt_read_cnt)/Total_read_cnt * 100, 2)
    Total_filt_base_pct = round(float(Total_filt_base_cnt)/Total_base_cnt * 100, 2)

    list_output = []
    list_output.append(sample)
    list_output.append(str(Total_read_cnt))
    list_output.append(str(Total_filt_read_cnt))
    list_output.append(str(Total_filt_read_pct))
    list_output.append(str(Total_base_cnt))
    list_output.append(str(Total_filt_base_cnt))
    list_output.append(str(Total_filt_base_pct))
    
    with open(output_file, 'w') as output_fp:
        print >> output_fp, '\t'.join(list_output)
    return

def usage():
    message='''
python %s --project_path ~/project_path --input_path ~/input_data_path --sample_id A0001 --log_path ~/log_path -o ~/output_path -r option_parameter
    
-f, --first_in_fastq     : rawdata,  first in fastq file
-s, --second_in_fastq    : rawdata,  second in fastq file
-a, --first_filt_fastq   : filter,  first in fastq file
-l, --second_filt_fastq  : filter,  second in fastq file
-o, --output
-i, --sample_id

''' %sys.argv[0]
    print message

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--first_in_fastq')
    parser.add_argument('-s', '--second_in_fastq')
    parser.add_argument('-a', '--first_filt_fastq')
    parser.add_argument('-l', '--second_filt_fastq')
    parser.add_argument('-o', '--output')
    parser.add_argument('-i', '--sample_id')
    args = parser.parse_args()
    try:
        len(args.first_in_fastq) > 0

    except:
        usage()
        sys.exit(2)

    fastq_statistics_run(args.first_in_fastq, args.second_in_fastq, args.first_filt_fastq, args.second_filt_fastq, args.sample_id, args.output)

if __name__ == '__main__':
    main()

