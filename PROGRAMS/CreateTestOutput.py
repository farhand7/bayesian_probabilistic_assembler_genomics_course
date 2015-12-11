'''
Dan Adler, Farhan Damani
Computational Genomics Final

CreateTestOutput.py

COMMAND:
python CreateTestOutput.py input.fa output.fa

Given a file of data with A's,T's,C's,G's representing a genome, this file
creates a continuous sequence of reads by taking away extraneous whitespace.

The input file is a sequence of a,t,c,g of any case with any spacing.
The output file is a sequence without any spaces and call capitalized.
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
Turn file string into output string we want:
Ex: atcg gcta = ATCGGTCA
@param data is the inputted string not formatted currectly.
@return the string formatted what we want.
'''
def formatString(data):
    data_no_space = ''.join(data.split())
    return data_no_space.upper()



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
    data = openFile(sys.argv[1])
    # Format string
    output_string = formatString(data)
    # Output string
    outputFile(sys.argv[2], output_string)
main()