'''
Dan Adler, Farhan Damani
Computational Genomics Final

Get read metrics for each trial.
Used to see if reads were mapped to the right spot, consecutively (i.e trace a read through a sequence)

COMMAND:
python Read_Metrics.py PATH_TO_FILE/reads_trial_x.txt > output.txt

Data organized by:
true_offset,contig,offset,likelihood
'''

'''
Organize reads data by contig.

Essentailly places in subsequent dictionaries to organize by dict[contig][offset] = [[read, likelihood]...]

@ param the input file thas has been read as a giant string.
@ return the data structure
'''

# Imports
import sys
import numpy as np
import pandas as pd

def organize_data(fileIn):
    fileIn = fileIn.split('\n')
    reads_data = {}
    i = 0
    while len(fileIn[i]) > 0:
        r,s,o,l = fileIn[i].split(',')
        if int(s) not in reads_data:
            reads_data[int(s)] = {}
        if int(o) not in reads_data[int(s)]:
            reads_data[int(s)][int(o)] = []
        reads_data[int(s)][int(o)].append([int(r),float(l)])
        i += 1

    return reads_data

'''
Now take data and organize them incrementally

so we get data[contig] = [[offsets][reads at each offset][frequency info]].

@ param reads_data is the output from organize_data
@ return the structure
'''
def organize_data_2(reads_data):
    data_pairs = {}
    contigs = reads_data.keys()
    contigs.sort()
    for c in contigs:
        offsets = range(0,max(reads_data[c]))
        freq = [0] * max(reads_data[c])
        reads = [''] * max(reads_data[c])
        for o in offsets:
            # Go through each read in the offset
            if o in reads_data[c]:
                o_reads = []
                num = 0
                # Count how many reads exist in that offset
                for r in reads_data[c][o]:
                    o_reads.append(r[0])
                    num += 1
                # Update the frequency for that offset + 100
                curr = 0
                while (o + curr) < len(freq) and curr < 100:
                    freq[o + curr] += num
                    curr += 1
                reads[o] = o_reads
        data_pairs[c] = [offsets, freq, reads]
    return data_pairs

'''
Format data into a stem-leaf type chart.

@ param data_pairs is the output from 
'''
def format_data(data_pairs):
    table_data = {}
    for c in data_pairs:
        table_data[c] = []
        m = max(data_pairs[c][0])
        d = range(100,m+100,100)
        r = []
        curr = 0
        i = 0
        while i < len(data_pairs[c][0]):
            s = str(d[curr] - 100) + '-' + str(d[curr]) + ': '
            temp = []
            while i < d[curr] and i < len(data_pairs[c][2]):
                for j in data_pairs[c][2][i]:
                    temp.append(str(j))
                i += 1
            rm = s + ' '.join(temp)
            if len(rm) > 0:
                table_data[c].append(rm)
            curr += 1
    return table_data

'''
Print the output.
@ param table_data is the organized output from above
'''
def print_output(table_data):
    curr = 0
    for c in table_data:
        print('Contig: ' + str(curr))
        curr += 1
        for s in table_data[c]:
            print(s)
        print('\n')

def read_metrics():
    reads_data = organize_data(open(sys.argv[1], 'r').read());
    data_pairs = organize_data_2(reads_data)
    table_data = format_data(data_pairs)
    print_output(table_data)

read_metrics()