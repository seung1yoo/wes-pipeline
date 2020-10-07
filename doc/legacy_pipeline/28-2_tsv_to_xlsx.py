#!/usr/bin/python

import sys
import os
import argparse
import openpyxl

def read_tsv_to_worksheet(input_file, worksheet):
    with open(input_file) as input_fp:
        for line in input_fp:
            units = line.strip().split('\t')
            list_output = []
            for unit in units:
                list_output.append(str(unit))
            worksheet.append(list_output)

    return worksheet

def worksheet_set_bold_row(worksheet, row_number_list):
    for row_number in row_number_list:
        for cell in worksheet["%d:%d" % (row_number,row_number)]:
            cell.font = openpyxl.styles.Font(bold=True)

def worksheet_auto_colum_width(worksheet, set_width_list = []):
    for i, col in enumerate(worksheet.columns):
        max_length = 0
        column = col[0].column
        if i < len(set_width_list):
            worksheet.column_dimensions[column].width = set_width_list[i]
        else:
            for cell in col:
                if cell.coordinate in worksheet.merged_cells:
                    continue
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(str(cell.value))
                except:
                    pass
            adjusted_width = max_length * 1.1
            worksheet.column_dimensions[column].width = adjusted_width

def xlsx_read(input_file, output_file ,bold_header, sheetname):
    workbook = openpyxl.Workbook()
    workbook_sheet1 = workbook.active
    workbook_sheet1.title = sheetname

    workbook_sheet1 = read_tsv_to_worksheet(input_file, workbook_sheet1)
    if bold_header:
        worksheet_set_bold_row(workbook_sheet1, [1])
    # worksheet_auto_colum_width(workbook_sheet1)

    workbook.save("%s" %output_file)
    return

def usage():
    message='''
python %s

-i, --input      : excel xlsx file

##optional
-o, --output : default(output.xlsx)
''' %sys.argv[0]
    print message


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tsv', help='Input tsv')
    parser.add_argument('--xlsx', default="output.xlsx", help='Output xlsx name')
    parser.add_argument('--header', default=False, action='store_true', help='Bold output xlsx header')
    parser.add_argument('--sheetname', default="Sheet1", help='Output xlsx sheet name')
    args = parser.parse_args()
    try:
        len(args.tsv) > 0

    except:
        usage()
        sys.exit(2)

    xlsx_read(args.tsv, args.xlsx, args.header, args.sheetname)

if __name__ == '__main__':
    main()
