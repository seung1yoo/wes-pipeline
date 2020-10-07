#!/usr/bin/python

import sys
import glob
import os
import locale
import time
import argparse
import numpy
import subprocess
from pipeline import *

dic_ts_tv = {
        'A:G':'Transition',
        'G:A':'Transition',
        'C:T':'Transition',
        'T:C':'Transition',
        'A:T':'Transversion',
        'T:A':'Transversion',
        'A:C':'Transversion',
        'C:A':'Transversion',
        'G:C':'Transversion',
        'C:G':'Transversion',
        'G:T':'Transversion',
        'T:G':'Transversion',
        }

def comma_number(n):
    locale.setlocale(locale.LC_ALL, '')
    if str(n).find('.') != -1:
        decimal_cnt = len(str(n)) - str(n).find('.')
        format_str = '%' + '.%sf' % str(decimal_cnt)
    else:
        format_str = '%d'
    return locale.format(format_str, n, 1)

def MakeDirectory(Dir):
    if not os.path.exists(Dir):
        os.system('mkdir -p %s' %Dir)

def write_head(output_fp, seq_type):
    output_head = []
    output_head.append('Sample_ID')
    output_head.append('sequence_read')
    output_head.append('clean_read')
    output_head.append('clean_read_rate(%)')
    output_head.append('sequence_base')
    output_head.append('clean_base')
    output_head.append('clean_base_rate(%)')
    output_head.append('deduplication_read')
    output_head.append('deduplication_rate(%)')
    output_head.append('mapping_read')
    output_head.append('mapping_rate(%)')
    output_head.append('insert_size(mean)')
#    output_head.append('insert_size(median)')
    output_head.append('insert_size(std)')
    output_head.append('unique_read')
    output_head.append('unique_rate(%)')

    if seq_type == "WES":
        output_head.append('on-target_read')
        output_head.append('on-target_rate(%)')
        output_head.append('on-target_depth(mean)')
#        output_head.append('on-target_depth(median)')
        output_head.append('on-target_depth(std)')
        output_head.append('target_bp')
    else:
        output_head.append('depth(mean)')
#        output_head.append('depth(median)')
        output_head.append('depth(std)')
        output_head.append('reference_bp')

    output_head.append('coverage_1X_rate(%)')
    output_head.append('coverage_5X_rate(%)')
    output_head.append('coverage_10X_rate(%)')
    output_head.append('coverage_20X_rate(%)')
    output_head.append('coverage_50X_rate(%)')
    if seq_type == "WES":
        output_head.append('chrX_depth(mean)')
        output_head.append('chrY_depth(mean)')
        output_head.append('chrX/chrY_ratio')
    output_head.append('total_variants')
    output_head.append('total_snv')
    output_head.append('total_insert')
    output_head.append('total_deletion')
    output_head.append('novel_variants(dbSNP)')
    output_head.append('novel_variants_ratio')
    output_head.append('hetero_variants')
    output_head.append('homo_variants')
    output_head.append('hetero_homo_ratio')
    output_head.append('transition_variants')
    output_head.append('transversion_variants')
    output_head.append('ts_tv_ratio')
    print >> output_fp, '\t'.join(output_head)
    return

def fasta_statistics(reference):

    reference_size = 0

    with open(reference) as input_fp:
        for line in input_fp:
            if line.startswith('>'):
                continue
            else:
                line_sequence = line.strip()
                for allele in line_sequence:
                    if allele != "N" and allele != "n":
                        reference_size += 1

    print 'referecne size : %s' %reference_size
    return reference_size

def statistics_run(statistics_fastq_path, merge_flagstat, dedup_flagstat, uniqread_flagstat, coverage_file, depth_file, output_file, sample_id, reference, seq_type, samtools, snpeff_file, qualimap_file):

    Dic_statistics = {}

    list_statistics_fastq_file = glob.glob('%s/*' %statistics_fastq_path)

    if seq_type == "WGS":
        reference_size = fasta_statistics(reference)
    elif seq_type == "WES":
        reference_size = 0
    else:
        print 'ERROR : seq type : %s' %seq_type
        sys.exit()

    sequence_read_cnt       = 0 
    sequence_clean_read_cnt = 0
    sequence_base_cnt       = 0
    sequence_clean_base_cnt = 0

    for statistics_fastq_file in list_statistics_fastq_file:
        with open(statistics_fastq_file) as input_fp:
            for line in input_fp:
                units = line.strip().split('\t')
#            sample_id               = units[1]
                sequence_read_cnt       += int(units[1])
                sequence_clean_read_cnt += int(units[2])
#                cln_fastq_pct           = units[3]
                sequence_base_cnt       += int(units[4])
                sequence_clean_base_cnt += int(units[5])
#                cln_fastq_base_pct      = int(units[6])
   
    cln_fastq_pct      = round((sequence_clean_read_cnt / float(sequence_read_cnt)) * 100, 2)
    cln_fastq_base_pct = round((sequence_clean_base_cnt / float(sequence_base_cnt)) * 100, 2)

    Dic_statistics[sample_id] = []
    Dic_statistics[sample_id].append(sample_id)
    Dic_statistics[sample_id].append(str(sequence_read_cnt))
    Dic_statistics[sample_id].append(str(sequence_clean_read_cnt))
    Dic_statistics[sample_id].append(str(cln_fastq_pct))
    Dic_statistics[sample_id].append(str(sequence_base_cnt))
    Dic_statistics[sample_id].append(str(sequence_clean_base_cnt))
    Dic_statistics[sample_id].append(str(cln_fastq_base_pct))

# deduplication statistics
    fp = open(dedup_flagstat, 'r')
#    for line in fp.xreadlines():
#        units = line.strip().split()
#        if units[3] == 'in':
#            dedup_read_cnt = int(units[0])
#        if units[3] == "mapped":
#            mapped_read_cnt = units[0]
#            break
    for line in fp.xreadlines():
        units = line.strip().split('\t')
        if units[0] == 'SN':
            if units[1] == "raw total sequences:":
                dedup_read_cnt = int(units[2])
            elif units[1] == "reads mapped:":
                mapped_read_cnt = int(units[2])
            elif units[1] == "mismatches:":
                mismatches_cnt = int(units[2])
            elif units[1] == "error rate:":
                error_rate = float(units[2])
            elif units[1] == "bases mapped:":
                total_bases = float(units[2])
            elif units[1] == "insert size average:":
                insert_size_mean = round(float(units[2]), 2)
            elif units[1] == "insert size standard deviation:":
                insert_size_sd = float(units[2])
    fp.close()

    cln_read_cnt = int(Dic_statistics[sample_id][2])
    dedup_pct = round(dedup_read_cnt/float(cln_read_cnt) *100, 2)
    mapped_rate = round(float(mapped_read_cnt) / float(dedup_read_cnt) * 100, 2)

    Dic_statistics[sample_id].append(str(dedup_read_cnt))
    Dic_statistics[sample_id].append(str(dedup_pct))
    Dic_statistics[sample_id].append(str(mapped_read_cnt))
    Dic_statistics[sample_id].append(str(mapped_rate))
    Dic_statistics[sample_id].append(str(insert_size_mean))
    Dic_statistics[sample_id].append(str(insert_size_sd))
#    Dic_statistics[sample_id].append(str(total_bases))

# insertion size statistics

    set_id = set()
    dic_id = {} 
    list_insert_size = []

#    cmd = '%s view %s' % (samtools, dedup_bam)
#    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#    for line in proc.stdout:
#        units = line.strip().split('\t')
#        chr_ID = units[2]
#        bp = units[3]
#        insert_size = int(units[8]) 
#        read_id = units[0] 
#        if units[6] == '=':
#            if insert_size > 0 and insert_size < 1000: 
#                list_insert_size.append(insert_size)

#    avg_insert =  numpy.mean(list_insert_size) 
#    median_insert = numpy.median(list_insert_size)
#    sd_insert = numpy.std(list_insert_size)
#    min_insert = numpy.min(list_insert_size)
#    max_insert = numpy.max(list_insert_size)
#    insert_1q = numpy.percentile(list_insert_size, 25)
#    insert_3q = numpy.percentile(list_insert_size, 75)
#    Dic_statistics[sample_id].append(str(avg_insert))
#    Dic_statistics[sample_id].append(str(median_insert))
#    Dic_statistics[sample_id].append(str(sd_insert))

# UniqRead
    fp = open(uniqread_flagstat, 'r')
#    for line in fp.xreadlines():
#        units = line.strip().split()
#        if units[3] == 'in':
#            uniq_read = int(units[0])
    for line in fp.xreadlines():
        units = line.strip().split('\t')
        if units[0] == 'SN':
            if units[1] == "raw total sequences:":
                uniq_read_cnt = int(units[2])
    fp.close()

    uniq_rate = round(float(uniq_read_cnt) / float(dedup_read_cnt) * 100, 2)

    Dic_statistics[sample_id].append(str(uniq_read_cnt))
    Dic_statistics[sample_id].append(str(uniq_rate))

# Average Depth
    total_depth = 0
    depth_1X = 0
    depth_2X = 0
    depth_5X = 0
    depth_10X = 0
    depth_20X = 0
    depth_30X = 0
    depth_50X = 0

    if seq_type == "WES":
        list_depth = []
        list_chrX_depth = []
        list_chrY_depth = []

        fp = open(depth_file, 'r')
        for line in fp.xreadlines():
            units = line.strip().split()
            chr_ID = units[0]
            start_bp = units[1]
            end_bp = units[2]
            depth = int(float(units[-1]))
            list_depth.append(depth)
            if depth >= 1:
                depth_1X += 1
            if depth >= 2:
                depth_2X += 1
            if depth >= 5:
                depth_5X += 1
            if depth >= 10:
                depth_10X += 1
            if depth >= 20:
                depth_20X += 1
            if depth >= 30:
                depth_30X += 1
            if depth >= 50:
                depth_50X += 1
            total_depth += int(depth)

            if chr_ID == "chrX":
                list_chrX_depth.append(depth)
            if chr_ID == "chrY":
                list_chrY_depth.append(depth)

            reference_size += 1
#        print total_depth, reference_size
        fp.close()
#    print total_depth, reference_size
        avg_depth    = round(numpy.mean(list_depth), 2)
        median_depth = round(numpy.median(list_depth), 2)
        std_depth    = round(numpy.std(list_depth), 2)

#    avg_depth     = round(float(total_depth) / reference_size, 2)
        chrX_avg_depth    = round(numpy.mean(list_chrX_depth), 2)
        chrY_avg_depth    = round(numpy.mean(list_chrY_depth), 2)
        chrX_chrY_ratio    = round(chrX_avg_depth / chrY_avg_depth, 2)

        depth_1X_avg  = round(depth_1X / float(reference_size) * 100, 2)
        depth_5X_avg  = round(depth_5X / float(reference_size) * 100, 2)
        depth_10X_avg = round(depth_10X / float(reference_size) * 100, 2)
        depth_20X_avg = round(depth_20X / float(reference_size) * 100, 2)
        depth_50X_avg = round(depth_50X / float(reference_size) * 100, 2)

        Ind_target_read_cnt = 0

        fp = open(coverage_file, 'r')
        for line in fp.xreadlines():
            units = line.strip().split()
            chr_ID = units[0]
            start_bp = units[1]
            end_bp = units[2]
            read_cnt = units[-4]
            on_bp_cnt = units[-3]
            target_bp_cnt = units[-2]

            Ind_target_read_cnt += int(read_cnt)
        fp.close()

        ontarget_read_pct = str(round(Ind_target_read_cnt/float(dedup_read_cnt) *100, 2))
        Dic_statistics[sample_id].append(str(Ind_target_read_cnt))
        Dic_statistics[sample_id].append(str(ontarget_read_pct))

    else:
        with open(qualimap_file) as input_fp:
            for line in input_fp:
                if line.strip().startswith('mean coverageData'):
                    units = line.strip().split(' = ')
                    avg_depth = units[1].replace('X', '')
                if line.strip().startswith('std coverageData'):
                    units = line.strip().split(' = ')
                    std_depth = units[1].replace('X', '')
                if line.strip().startswith('number of bases'):
                    units = line.strip().split(' = ')
                    reference_size_include_N = units[1].replace(',', '').split()[0]
                if line.strip().startswith('There'):
                    units = line.strip().split(' >= ')
                    coverage_depth = units[1]
                    if coverage_depth == "1X":
                        depth_1X = line.strip().split()[3].replace('%', '')
                        depth_1X_base = int(float(depth_1X) * int(reference_size_include_N) / 100)
                        depth_1X_avg  = round(depth_1X_base / float(reference_size) * 100, 2)
#                        print >> sys.stderr, 'depth_1X      : %s' %depth_1X
#                        print >> sys.stderr, 'depth_1X_base : %s' %depth_1X_base
#                        print >> sys.stderr, 'depth_1X_avg  : %s' %depth_1X_avg
                    if coverage_depth == "5X":
                        depth_5X = line.strip().split()[3].replace('%', '')
                        depth_5X_base = int(float(depth_5X) * int(reference_size_include_N) / 100)
                        depth_5X_avg  = round(depth_5X_base / float(reference_size) * 100, 2)
                    if coverage_depth == "10X":
                        depth_10X = line.strip().split()[3].replace('%', '')
                        depth_10X_base = int(float(depth_10X) * int(reference_size_include_N) / 100)
                        depth_10X_avg  = round(depth_10X_base / float(reference_size) * 100, 2)
                    if coverage_depth == "20X":
                        depth_20X = line.strip().split()[3].replace('%', '')
                        depth_20X_base = int(float(depth_20X) * int(reference_size_include_N) / 100)
                        depth_20X_avg  = round(depth_20X_base / float(reference_size) * 100, 2)
                    if coverage_depth == "50X":
                        depth_50X = line.strip().split()[3].replace('%', '')
                        depth_50X_base = int(float(depth_50X) * int(reference_size_include_N) / 100)
                        depth_50X_avg  = round(depth_50X_base / float(reference_size) * 100, 2)
                

    Dic_statistics[sample_id].append(str(avg_depth))
#    Dic_statistics[sample_id].append(str(median_depth))
    Dic_statistics[sample_id].append(str(std_depth))
    Dic_statistics[sample_id].append(str(reference_size))
    Dic_statistics[sample_id].append(str(depth_1X_avg))
    Dic_statistics[sample_id].append(str(depth_5X_avg))
    Dic_statistics[sample_id].append(str(depth_10X_avg))
    Dic_statistics[sample_id].append(str(depth_20X_avg))
    Dic_statistics[sample_id].append(str(depth_50X_avg))
    if seq_type == "WES":
        Dic_statistics[sample_id].append(str(chrX_avg_depth))
        Dic_statistics[sample_id].append(str(chrY_avg_depth))
        Dic_statistics[sample_id].append(str(chrX_chrY_ratio))

# SNV 
    total_var    = 0
    count_novel  = 0
    count_snp    = 0
    count_ins    = 0
    count_del    = 0
    count_hetero = 0
    count_homo   = 0
    count_ts     = 0
    count_tv     = 0
    dic_snpeff_head = {}
    set_snp = set()

    with open(snpeff_file) as input_fp:
        for line in input_fp:
            units = line.strip().split('\t')
            if line.startswith("CHROM"):
                head = units
                for idx in range(len(head)):
                    col_name = head[idx]
                    dic_snpeff_head[col_name] = idx
                    dic_snpeff_head[idx] = col_name
                    if 'GT' in col_name:
                        genotype_idx = idx

            else:
                chr_idx = dic_snpeff_head['CHROM']
                pos_idx = dic_snpeff_head['POS']
                rs_idx = dic_snpeff_head['ID']
                ref_idx = dic_snpeff_head['REF']
                alt_idx = dic_snpeff_head['ALT']
                var_type_idx = dic_snpeff_head['VARTYPE']

                chr_id  = units[chr_idx]
                bp      = units[pos_idx]
                rs      = units[rs_idx]
                ref      = units[ref_idx]
                alt      = units[alt_idx]
                vartype = units[var_type_idx]

######### 
                bp_key = '%s:%s' %(chr_id, bp)
                if bp_key in set_snp:
                    continue
                else:
                    set_snp.add(bp_key)

######### 
                if vartype == "SNP":
                    ts_tv_key = '%s:%s' %(ref, alt)
                    if dic_ts_tv[ts_tv_key] == "Transition":
                        count_ts += 1
                    if dic_ts_tv[ts_tv_key] == "Transversion":
                        count_tv += 1

                hetero_homo_type = ''

                genotype = units[genotype_idx]

                if '/' in genotype:
                    allele_1 = genotype.split('/')[0]
                    allele_2 = genotype.split('/')[1]
                if '|' in genotype:
                    allele_1 = genotype.split('|')[0]
                    allele_2 = genotype.split('|')[1]

                if allele_1 == allele_2:
                    count_homo += 1
                else:
                    count_hetero += 1

                total_var += 1
                if rs == ".":
                    count_novel += 1
                if vartype == "SNP":
                    count_snp +=1 
                if vartype == "INS":
                    count_ins +=1 
                if vartype == "DEL":
                    count_del +=1 
    
    novel_rate       = round(int(count_novel) / float(total_var), 2)
    hetero_homo_rate = round(int(count_hetero) / float(count_homo), 2)
    ts_tv_rate = round(int(count_ts) / float(count_tv), 2)

    Dic_statistics[sample_id].append(str(total_var))
    Dic_statistics[sample_id].append(str(count_snp))
    Dic_statistics[sample_id].append(str(count_ins))
    Dic_statistics[sample_id].append(str(count_del))
    Dic_statistics[sample_id].append(str(count_novel))
    Dic_statistics[sample_id].append(str(novel_rate))

    Dic_statistics[sample_id].append(str(count_hetero))
    Dic_statistics[sample_id].append(str(count_homo))
    Dic_statistics[sample_id].append(str(hetero_homo_rate))

    Dic_statistics[sample_id].append(str(count_ts))
    Dic_statistics[sample_id].append(str(count_tv))
    Dic_statistics[sample_id].append(str(ts_tv_rate))

    output_fp = open(output_file, 'w')
    write_head(output_fp, seq_type)
    print >> output_fp, '\t'.join(Dic_statistics[sample_id])
    print >> sys.stderr, '%s  Complete : %s' %(Present_time(), sample_id)

    End_time = Present_time()
    print 'Preccess End : %s' %End_time
    output_fp.close()
    return

def usage():
    message='''
python %s --project_path ~/project_path --input_path ~/input_data_path --sample_id A0001 --log_path ~/log_path -o ~/output_path -r option_parameter
    
-f, --fastq_statistics     : rawdata,  first in fastq file
-m, --merge_flagstat
-b, --dedup_bam
-d, --dedup_flagstat
-u, --uniqread_flagstat
-k, --coverage
-q, --qualimap
-c, --depth
-o, --output
-s, --sample_id
-r, --reference
-p, --program

''' %sys.argv[0]
    print message

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastq_statistics_path')
    parser.add_argument('-m', '--merge_flagstat')
#    parser.add_argument('-b', '--dedup_bam')
    parser.add_argument('-d', '--dedup_flagstat')
    parser.add_argument('-u', '--uniqread_flagstat')
    parser.add_argument('-k', '--coverage')
    parser.add_argument('-q', '--qualimap')
    parser.add_argument('-c', '--depth')
    parser.add_argument('-o', '--output')
    parser.add_argument('-s', '--sample_id')
    parser.add_argument('-n', '--seq_type')
    parser.add_argument('-r', '--reference')
    parser.add_argument('-p', '--program')
    parser.add_argument('-g', '--snpeff')
    args = parser.parse_args()
    try:
        len(args.fastq_statistics_path) > 0
        len(args.merge_flagstat) > 0
        len(args.dedup_flagstat) > 0
        len(args.uniqread_flagstat) > 0

    except:
        usage()
        sys.exit(2)

    statistics_run(args.fastq_statistics_path, args.merge_flagstat, args.dedup_flagstat, args.uniqread_flagstat, args.coverage, args.depth, args.output, args.sample_id, args.reference, args.seq_type, args.program, args.snpeff, args.qualimap)
    return

if __name__ == '__main__':
    main()

