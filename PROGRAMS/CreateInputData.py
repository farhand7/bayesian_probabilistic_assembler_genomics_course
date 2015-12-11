'''
Dan Adler, Farhan Damani
Computational Genomics Final

CreateTestOutput.py

COMMAND:
python CreateInputData.py k input.fa output.fa k

Takes a genome (all caps no spaces) and
samples all k-mers from that genome from that genome.
'''

# Imports
import sys

'''
Opens file into a string.
@param infile is the file to be opened
@return the output string (string from the infile)
'''
def openFile(infile):
    data = open(infile, 'r').read()
    return data

'''
Samples n random reads of length k from that genome (all kmers).
@param data is the input genome.
@param k is the length of reads to sample.
@return reads of that data.
'''
def sampleReads(data, k):
    reads = []
    for i in range(len(data)-k+1):
        reads.append('>' + str(i+1))
        reads.append(data[i:i+100])
    return '\n'.join(reads)

'''
Output string we want to file
@param outfile name is the name of the output file
@param output_string is the otuput string
'''
def outputFile(outfile, output_string):
    of = open(outfile, 'w')
    of.write(output_string)
    of.close()

'''
Main function
'''
def main():
    # Open file
    data = openFile(sys.argv[2])
    # Get random reads
    reads = sampleReads(data,int(sys.argv[1]))
    # Output file
    outputFile(sys.argv[3],reads)
main()