Dan Adler
Farhan Damani

Computational Genomics 
600.439/639
Final Project

FOLDER ORGANIZATION
________________________

1) DATA:

Contains the original input utilized as inputs to our algorithm.  Actual data in .fa (or .txt) format is found in the main DATA folder.  Each line of these files contains a read indicator (filler line >x) and then the following line contains a read.  For our results, most of the data utilized the 1_wuhan_input.txt and 1_wuhan_input_ordered.txt, which are random orderings of the first 612 k-mers of the Wuhan aphid virus, taken from the NCBI reference sequence (NC_028388).

More will be explained on the creation of these input files in the PROGRAMS section.

Subfolders contain:

a) Original_Data: i.e. .data files which are original NCBI genome data files from different viruses, and .fa files (though not truly .fa files) which contains these genomes with extraneous information removed.  All .data files are sourced to their original downloadable locations, but as a reference:
m13mp18.data: https://www.neb.com/products/n4040-m13mp18-single-stranded-dna
wuhanaphivirus.data: http://www.ncbi.nlm.nih.gov/nuccore/952002249
papamelieravirus.data: http://www.ncbi.nlm.nih.gov/nuccore/NC_028378.1

b) Genomes: The Original_Data/_.fa files placed in a format with spaces removed and nucleotides capitalized for coding conveniences.

b) Kmers: The k-mers from genomes in the Genomes/_.fa folder.  Each line is either an identifier to the offset where that k-mer starts, or the k-mer itself.  For example, the line following ‘>6’ is the k-mer of length 100 that starts at offset 6 in the genome specified by the filename.  Note that all k-mers taken in this assignment were of length 100 for convenience.  Files specified by _x.fa are k-mer filers with the first x% reads of the genome were taken.

________________________


2) PROGRAMS

In this file is all the code utilized to a) generate data to input into our algorithm, b) run the algorithm itself, c) test pieces of the algorithm, and d) analyze data from our algorithm during each step and after termination.

a) Data to generate input.  The following files were used to generate input data:

i) CreateTestOutput.py can be run with the following command:

COMMAND:
python CreateTestOutput.py input.fa output.fa

It takes a _.fa file in the DATA/Original_Data folder (with just a genome that might have extraneous space, not all caps, etc.) and creates the .fa files in the DATA/Genomes folder by simply removing all spaces and capitalizing letters.  It is called ‘CreateTestOutput.py’ because the genome files are the ideal assembly of the reads.

ii) CreateTestReads.py can be run with the following command:

COMMAND:
python CreateInputData.py input.fa output.fa k

It takes in a _.fa file from the DATA/Genomes folder and creates a k-mer file located in the DATA/Kmers folder.  The parameter, k, specifies the length of the k-mer you will take, and each k-mer is sampled in order.

iii) CreateTestReads.py can be run with the following command:

COMMAND:
python CreateTestReads.py e input.txt output.txt

This file takes the k-mer data from DATA/Kmers and does two things.  It first randomly selects each k-mer and shuffles them.  Second, it introduces e% error into the data set by randomly choosing e% of the nucleotides, and sampling from a uniform (A,T,C,G) distribution to replace it.  This counts as creating our ‘sequencing’ error.  After running this program, the input data for the algorithm is generated.  

Note that the input.txt file is not the k-mers themselves.  An input.txt file for this program is a text file that specifics on individual lines k-mer file paths where data should be taken to make the randomized sample.  Thus, metagenomic sets can be created from individual k-mer files.

b) Algorithm Files:

The following files are used to run the algorithm.  Specifics about the algorithm can be found in our write-up, but brief explanations follow.

i) Driver.py.  Simply can be run by typing:

COMMAND:
python Driver.py PATH_TO_INPUT/input.txt PATH_TO_OUTPUT_FOLDER iterations

Ex:
python Driver.py ../DATA/1_wuhan_input_ordered.txt ../OUTPUT 20

This runs the algorithm using the data specified by your input.txt in the first line of the main.  It simply randomly creates contigs, and then under a certain number of iterations (specified by NUM_ITERS) maps the reads to their most likely alignment, computes new contigs, checks what needs to merged, and updates the likelihood.  It saves output data (described below) on the contigs, likelihood and read mappings at each iteration in the output directory specified by argument 2.  Lastly, the third argument specifies the number of iterations the algorithm will conduct, which should change depending on the size of your dataset.  The following files then are called by the Driver.

ii) Initialization.py

This file takes in k-mers and creates a dictionary of reads with those k-mers.  The parameters at the top define global parameters used throughout the algorithm.  CONTIG_LENGTH sets of a uniform(0,CONTIG_LENGTH) distribution for contains, NUM_CONTIGS is the guess of the ‘expected’ contains that we should have, and K_MER is the K_MER length.  DEFAULT_PROB is the initial probability, which is set to a small value.

The reads dictionary outputted will be utilized through the algorithm.  Each entry of the dictionary is a read.  The reads map to a list of lists.  Each individual list specifies a single contig (s), offset in the contig (o), negative log likelihood of it being in the current position (p), and its true position in the original genome, which will be used later for analysis purposes.

In addition, contigs are randomly created by mapping reads to random locations, and then computing a consensus sequence after initialization (more explained below).

iii) Read_Mapping.py are the functions used to initialize a read mapping step.  A full read mapping can be run using the method:

Read_Mapping.run(reads_dict, contigs) where reads_dict is a read_dict as specified above, and contigs is the current contigs list.  Essentially, each read is aligned to its optimal alignment in the contigs, and we compute the joint probability of all variables at that alingment position (i.e. the likelihood). If the likelihood at that position is greater than the likelihood for the same read at the previous time step, then we map the read to this new position. If the new likelihood is not greater, than we use an approximate inference procedure to accept/reject the new mapping (see write-up for details on this). Other methods define various distributions used to compute probabilities. Lastly, we compute the optimal alignment using the hamming distance.

The output is an updated dictionary with the new mappings to all the latent variables.

iv) Consensus_Sequence.py goes is used to update the contigs. It can be run using the method:

Consensus_Sequence.run_consensus(reads)

It essentially goes through the reads dictionary and finds the most likely nucleotides in every position of the contigs, and updates the contigs as follows.  The new contig list is outputted.

v) Merge_Contigs.py is used to check the current contigs list and see if two contigs overlap enough (exactly 15 or more nucleotides at the suffix of one contig matches with the prefix of another).  It can be run with using:

Merge_Contigs.run_merge(contigs, reads, OVERLAP_LENGTH = 15)

where contigs are your current contigs, reads is your reads dictionary, and OVERLAP_LENGTH (defaulted to 15) is the overlap that needs to occur.  If no merge can occur the function outputs an empty list, else it outputs a new contigs list and new reads dictionary updated by the contig merge.  Note that this file utilized code from Dan Adler Homework 4 assignment.  We also utilized Ben’s github suffixPrefixMatch method from the ‘Unpaired assembly challenge’, which can be found here:

http://nbviewer.ipython.org/github/BenLangmead/comp-genomics-class/blob/master/projects/UnpairedAsmChallenge.ipynb

vi) Likelihood.py can be run by utilizing the following:

Likelihood._likelihood(reads)

This function simply compute an updated likelihood of the reads data with their current mappings by summing through all the -log(p’s) for each read and mapping.  Note that we want the computed value to go down, since -log(p)’s summed have higher likelihood when they are smaller.  The sum is outputted.

c) Testing files

i) Test_Consensus_Sequence.py provides unittests that validate whether the Consensus_Sequence.py is working properly.  It can simply be run:

python Test_Consensus_Sequence.py

ii) Test_Merge_Contigs.py provides unittests that validate whether the Merge_Contigs.py file is working properly.  It can be simply run:

python Test_Merge_Contigs.py

d) Analysis files

i) QualityMetrics.py can be run during Driver.py to create plot figures to evaluate the percent change over contigs over time steps and to create the boxplot and other figures. The percent change function uses edit distance from Ben’s Github to do so from:

http://nbviewer.ipython.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_DP_EditDist.ipynb

One can simply run the file by typing:

python Quality_Metrics.py

It takes a change_in_contigs_data.txt folder from the SRC_OUTPUT/trial* directory (explained below), and computes a percent change over time figure for that trial contig data.  The data generated was taken over time during each step of the Driver.py file (see code).

ii) Read_Metrics.py

This file looks at data collected during run for each read dictionary.  The data collected can be found in the SRC_OUTPUT/Dan_trial*.  This data was taken from the wuhan virus (see above).  The code can be run using the command:

python Read_Metrics.py input.txt > output.txt

The input is a reads_trial_x.txt file (See Dan_trial_two), and the data in each of the x.txt files is on each line represents a read: true_offset,current contig, current offset, likelihood, on step x of the algorithm.  All of the reads_trial_x.txt files then represent the read mapping on each step.  Running the code produces a text file that bins what contigs/offsets reads were mapped to during a certain algorithmic step.  One can see an example of this output in the SRC_OUTPUT/Dan_trial_two/Dan_trial_two_reads_20.txt file.
 
4) SRC_OUTPUT

This is where the output data, taken during running of the algorithm, was collected.  

Essentially, folders specified by a trial contain a contig.txt file that shows each contig during the step of the algorithm. 

The other folders, Dan_trial_x.txt contains the contigs (in contigs.txt), likelihood (in likelihood.txt) and read mapping in each algorithmic step.  contigs.txt and likelihood.txt contain all steps in one file, while each read mapping is split up into its now file reads_trial_x.txt.  There are different amounts of read mapping files depending upon how many trials were run of the algorithm.  

The last file type is the Dan_trial_y_reads_x.txt file, which are outputs of Read_Metrics.py showing binned data of where each reads are mapped to.  Each section shows a contig, and the numbers i-j: represent offsets of that contig where the reads to their right were mapped to.

5) figures 

Contains various figures utilized within our presentation and write-up.


=============================================================================================
Group members did:

Dan Adler helped with the initial modeling of the algorithm to actually create the simple distributions to be used.  He came up with the idea of utilizing a geometric prior for the number of contigs distribution, and created/found all of the input data that was used for running the algorithm and the analysis.  Dan also was in charge of writing/testing the merge and consensus sequence steps of the algorithm, and also analyzed the read mapping data to see how reads were distributed within different portions of the algorithm.  Lastly, Dan created the README.txt and wrote various portions of the report.

Farhan Damani came up with the idea and found the Laserson paper. After setting up a meeting with Dr. Alexis Battle, Farhan incorporated Dr. Battle’s advice and later designed the graphical model, including assumptions, priors, and dependencies. Farhan explored approximate inference techniques to make probability computations tractable and identified the key algorithmic steps from the Laserson paper that should be used in the model. He wrote out the math that was necessary for implementation and coded the initialization, driver, and probabilistic algorithmic step (read mapping), requiring likelihood computations, local alignment, sampling through approximate inference techniques, and a scoring scheme to compute the joint probability of x and y. Lastly, he analyzed likelihood convergence, how the contigs evolved over time-steps, and wrote significant portions of the write-up and the slides.


