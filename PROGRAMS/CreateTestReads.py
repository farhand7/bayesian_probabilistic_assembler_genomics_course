'''
Dan Adler, Farhan Damani
Computational Genomics Final

CreateTestReads.py

COMMAND:
python CreateTestReads.py error_percent input.txt output.txt

Given a file containing kmer file names, we take in each file of kmers, 
shuffle them, and then create errors based upon an inputted error_percent
number between 0 and 100.
'''
# Imports
import sys
import random
import math

'''
Opens file into a string.
@param infile is the file to be opened
@return the output string (string from the infile)
'''
def openFile(infile):
    data = open(infile, 'r').read()
    return data

'''
Read genome data
@param genome is the genome to read.
@return a list of the genomes
'''
def readGenomeData(genome):
    g = genome.split('\n')
    output_genome = []
    for i in g:
        if i[0] != ">":
            output_genome.append([i,curr])
        else:
            curr = i[1:]
    return output_genome

'''
Add error to a genome.
@param p is the percent error we want
@param genomes is the genome kmers
@return the genomes with error
'''
def addError(p, genomes):
    number_error = int(math.floor(p/100 * len(genomes)* len(genomes[0])))
    for i in range(number_error):
        nts = ['A', 'T', 'C', 'G']
        # Choose random genome
        g_number = random.randint(0,len(genomes) - 1)
        # Choose random offset
        nt_number = random.randint(0,len(genomes[0]) - 1)
        nt = genomes[g_number][0][nt_number]
        nts.remove(nt)
        replacement_nt = random.choice(nts)
        g = list(genomes[g_number][0])
        g[nt_number] = replacement_nt
        genomes[g_number][0] = ''.join(g)
    return genomes

'''
Format output
@param genome is an error induced genome.
@return the formatted output string
'''
def formatOutput(genome):
    g = ''
    for i in range(len(genome)):
        g += '>' + str(genome[i][1]) + '\n'
        g += genome[i][0] + '\n'
    return g

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
    filenames = openFile(sys.argv[2]).split('\n')
    genomes = []
    for i in range(len(filenames)):
        # Open k-mer file
        genome = openFile(filenames[i])
        # Filter out read numbers
        genomes.extend(readGenomeData(genome))
    # Shuffle genomes
    random.shuffle(genomes)
    # Create error
    error_genomes = addError(float(sys.argv[1]),genomes)
    # Format output
    error_genomes = formatOutput(error_genomes)
    # Output reads
    outputFile(sys.argv[3],error_genomes)
main()