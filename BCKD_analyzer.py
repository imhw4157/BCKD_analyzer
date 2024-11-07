import os
import regex
from collections import Counter
import pandas as pd

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows


ind = os.getcwd()
f_list = open('target_list.txt').readlines()
direc = f_list[0].strip()

WT_ref_seq = "TCTTTGTCGCAAGGTGGTTAGAAGCAGATAGCGTGTTGGCCTAGCGGTGCTGCTATTTGAATAAGGCATCTGGTGGGACAGCCGTATGTCTATCCTGTGAGAAAGGAGGAGTGGTTATTAGCCAAGATGGCGATAGTGATGTCTGCGGCCAAGATCTGGAGGCCGAGCCGTGGCCTGCGCCAGGCTGCTCTTCTCCTGTTG"
left_indicator = "TTTGTCGCAAGGTGG"
right_indicator = "GCTGCTCTTCTCCTG"
right_indicator2 = "GGCTGCTCTTCTCCT"
attB = "GGCTTGTCGACGACGGCGGTCTCCGTCGTCAGGATCAT"
WT_marker1 = "TGGTTATTAG"
integ_indicator = "GCGGTAGCGATCGC"


def minimum_frequency(records, n=2):
    seq_form_fastq = [record[1] for record in records] 
    counts = Counter(seq_form_fastq)  
    return [record for record in records if counts[record[1]] >= n] 


def pick_indicator_with_mismatch(seq_input):
    indicator_records = []
    for seq in seq_input:
        is_WT_L = regex.search(rf"({left_indicator}){{s<={2}}}", seq[1]) 
        is_WT_R = regex.search(rf"({right_indicator}){{s<={2}}}", seq[1])
        is_WT_R2 = regex.search(rf"({right_indicator2}){{s<={2}}}", seq[1])        
        if is_WT_L and (is_WT_R or is_WT_R2):
            indicator_records.append(seq) 
    return indicator_records 


def count_integ(seq_input): 
    cnt = {'dimer': 0, 'attB': 0, 'WT': 0, 'integration': 0, 'trash': 0}
    for seq in seq_input: 
        if integ_indicator in seq[1]:
            cnt['integration'] += 1
        elif regex.search(rf"({attB}){{s<={2}}}", seq[1]):
            cnt['attB'] += 1
        elif WT_marker1 in seq[1]:
            cnt['WT'] += 1
        else:
            cnt['trash'] += 1
    return cnt


wb = Workbook()
sheet = wb.active
sheet.title = "Analysis Result"

sheet_headers = ['Index', 'Min.Freq', 'WT', 'attB', 'Integration', 'Trash']
sheet.append(sheet_headers)

for t in f_list[1:]:
    each_line_list = t.split('\t')
    if len(each_line_list) <= 1:
        continue
    
    file = each_line_list[0]
    file_name = f"{file}.fastqjoin"
    os.chdir(os.path.join(ind, direc))

    records = []
    with open(file_name, "r") as handle:
        for record in FastqGeneralIterator(handle):
            records.append(record)

    record_indicator = pick_indicator_with_mismatch(records)
    record_minimum = minimum_frequency(record_indicator, n=2)
    cnt = count_integ(record_minimum)

    result_row = [file, len(record_indicator), cnt['WT'], cnt['attB'], cnt['integration'], cnt['trash']]
    sheet.append(result_row)
    os.chdir(ind)

    print(f"\u25B6 #{file_name} analysis completed.")
    print("--------------------------------------------------------------------------")

wb.save('sub.xlsx')
print('Jobs done! made by HW')
