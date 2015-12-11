'''
Dan Adler, Farhan Damani
Computational Genomics Final

Driver Program

Runs the Bayesian Alignment algorithm.  Simply:
1) Initalizes reads in a map and randomly distributes reads to contigs
2) Maps the reads to the most likely alignment
3) Performs a new consensus sequence
4) Checks whether reads should merge and merges if needed

COMMAND:
python Driver.py PATH_TO_INPUT/input.txt PATH_TO_OUTPUT_FOLDER iterations

Ex:
python Driver.py ../DATA/1_wuhan_input_ordered.txt ../OUTPUT 20

Note:
If you want to create the figure in the qm likelihood file you MUST comment in the following lines
#likelihood_list = []
#likelihood_list.append(likelihood) (many lines do this)
#qm._create_likelihood_figure(likelihood_list)
'''

# IMPORTS
import Initialization as init
import Consensus_Sequence as cs
import Read_Mapping as rm
import Merge_Contigs as mc
import Likelihood as ll
import Quality_Metrics as qm
import sys

# NUMBER OF ITERATIONS
NUM_ITERS = int(sys.argv[3])
def _main():

    # Open files and process reads dictionary
    f = open(sys.argv[1], 'r') 
    reads_dict = init._process(f)

    # Get contigs from first consensus sequence
    contigs = cs.run_consensus(reads_dict)
    contig_file = open(sys.argv[2] + '/contig.txt', 'w+')
    ll_file = open(sys.argv[2] + '/likelihood.txt', 'w+')

    # Set initial parameters 
    likelihood = 0
    likelihood_new = 0
    #likelihood_list = []

    for i in range(NUM_ITERS):

        '''FILE WRITES'''
        # Contigs file write data
        contig_file.write('%s\tstart\t' %(str(i)))
        for c in contigs:
            contig_file.write('%s\t' %(str(c)))
        contig_file.write('\n')
        contig_file.flush()
        # Likelihood file write data
        ll_file.write('%s\t%s\t%s\n' %(str(i), str(likelihood), str(len(contigs)))), ll_file.flush()
        #likelihood_list.append(float(likelihood))
        # Reads file write data
        reads_file = open(sys.argv[2] + '/reads_trial_' + str(i) + '.txt','w')
        for r in reads_dict:
            for l in reads_dict[r]:
                reads_file.write(str(l[3])+','+str(l[0])+','+str(l[1])+str(',')+str(l[3])+'\n')
        reads_file.close()
        '''COMPUTATION OF ALGORITHM'''
        # Update likelihood
        likelihood = likelihood_new
        # Map reads
        reads_dict = rm.run(reads_dict, contigs)
        # Run Consensus Sequence
        contigs = cs.run_consensus(reads_dict)
        # Print data to file
        contig_file.write('%s\tmerge\t' %(str(i)))
        for c in contigs:
            contig_file.write('%s\t' %(str(c)))
        contig_file.write('\n')
        # Run merge
        contigs, reads_dict = mc.run_merge(contigs,reads_dict) # how do we know if a merge has happened..do we need to know?
        # Get new likelihood
        likelihood_new = ll._likelihood(reads_dict,contigs)


    '''FILE WRITES'''
    # Reads file write data
    reads_file = open(sys.argv[2] + '/reads_trial_' + str(i+1) + '.txt','w')
    for r in reads_dict:
        for l in reads_dict[r]:
            reads_file.write(str(l[3])+','+str(l[0])+','+str(l[1])+str(',')+str(l[3])+'\n')
    reads_file.close()
    # Print data to file
    for c in contigs:
        contig_file.write('1000\tend\t%s\n' %(str(c)))
    ll_file.write('%s\t%s\t%s\n' %(str(NUM_ITERS), str(likelihood), str(len(contigs)))), ll_file.flush()
    #likelihood_list.append(likelihood)
    #qm._create_likelihood_figure(likelihood_list)

_main()