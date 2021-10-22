
######################

# Bianca De Sanctis
# October 18 2021

# Simulates three genomes named "true", "false" and "query", where query is assumed to be from a population closer to true than false, and calculates the probability of correctly assigning reads drawn from the query to the true or false sequence under different assignment protocols.

# Input: Parameters for simulation. Run with "-h" command to show all parameter options. Generally flexible but assumes diploidy, and that the reads are all the same length, but the latter can be relaxed easily with a wrapper script.

# Output: The number of query reads assigned correctly (to the true sequence), incorrectly (to the false sequence), or not assigned. Results can be printed or written to file.
# Optionally produces fasta files of the true sequence, false sequence and query reads for downstream input into classification engines such as Clark, Kraken etc, or for ancient DNA damage simulators such as gargammel.

# Requires msprime as the simulation engine.

####################### 
# Set up

import random
import msprime
import os
import sys
import argparse

# Read input
parser = argparse.ArgumentParser(description='Simulate three lineages and calculate assignment accuracy.')
parser.add_argument("-L", "--genome_length", help='length of reference sequences',default=15000,metavar='') 
parser.add_argument("-N", "--pop_size", default=10000, help='population size',metavar='')
parser.add_argument("-R", "--num_reads", default=1000, help='number of reads',metavar='')
parser.add_argument("-k", "--read_length", default=70, help='read length',metavar='')
parser.add_argument("-mu", "--mutation_rate", default=1e-7, help='mutation rate (per site per generation)', metavar='')
parser.add_argument("-Tqt", default=1, help='query-true divergence time (mya)', metavar='')
parser.add_argument("-Ttf", default=15, help='true-false divergence time (mya)', metavar='')
parser.add_argument("-gen", "--gen_length", default=20, help='generation length (years)', metavar='')
parser.add_argument("-qA", "--query_age", default=0, help='query sample age (generations)', metavar='')
parser.add_argument("-qT", "--true_age", default=0, help='true sample age (generations)', metavar='')
parser.add_argument("-qF", "--false_age", default=0, help='false sample age (generations)', metavar='')
parser.add_argument("-r", "--recomb_rate", default=1e-8, help='recombination rate (per site per generation)', metavar='')
parser.add_argument("-o", "--output_file", default="results.txt", help='output file name', metavar='')
parser.add_argument("--write_files", default=True, help='write output to a file (default True)', metavar='')
parser.add_argument("--print_results", default=True, help='print results to screen (default True)', metavar='')
parser.add_argument("--create_fastas", default=False, help='create reference sequence and query read fastas (default False)', metavar='')

args = parser.parse_args()

# Set parameters
seqlength = int(args.genome_length) # How long do you want your chromosome to be?
pop_size = int(args.pop_size)
num_reads = int(args.num_reads)
read_length= int(args.read_length)
mutation_rate = float(args.mutation_rate) # per site per generation
Tqt = float(args.Tqt) # in millions of years
Ttf = float(args.Ttf) # in millions of years
generation_length = float(args.gen_length) # in years
query_age = float(args.query_age) # ages are in generations
true_age = float(args.true_age)
false_age = float(args.false_age)
recomb_rate = float(args.recomb_rate)

# Where do you want the output to go, and do you want to create fastas?
write_files = args.write_files
file_path = args.output_file
print_results = args.print_results
create_fastas = args.create_fastas

####################
# Run

# Simulate an ancestral history for 3 diploid samples under the coalescent
# with recombination on a seqlength region with parameters defined above..
newick = "((true:" + str(Tqt) + ",query:" + str(Tqt) + "):" + str(Ttf-Tqt) + ",false:" + str(Ttf) + ")"
demography = msprime.Demography.from_species_tree(
    newick,
    time_units="myr",
    initial_size=pop_size,
    generation_time=generation_length)

ts = msprime.sim_ancestry(
	samples= [msprime.SampleSet(1, population="true", time=true_age),msprime.SampleSet(1, population="query", time=query_age), msprime.SampleSet(1, population="false", time=false_age)],
	demography=demography,
	recombination_rate=recomb_rate,
	ploidy=2,
	sequence_length=seqlength,
	random_seed=123456)

mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=5678)
# Default mutation model is msprime.JC69. 

# Create arbitrary sequence of same length, because msprime doesn't bother to simulate non-variable sites.
bases = ["A","C","T","G"]
bgseq = random.choices(bases,k=seqlength)
# Assumes each nucleotide occurs at 25% frequency.

# Now use msprime output to create true and false reference sequences, replacing the variable sites in bgseq. Recall python is 0-based but msprime is not.
trueseq = list(bgseq); falseseq = list(bgseq); queryseq1 = list(bgseq); queryseq2 = list(bgseq)
for var in mts.variants():
	pos = var.site.position
	# Here we are abritrarily picking one strand of the diploid true and false sequences to be our "reference", but sampling query reads from both strands of its sequence. 
	trueseq[int(var.site.position)-1] = var.alleles[var.genotypes[0]]
	falseseq[int(var.site.position)-1] = var.alleles[var.genotypes[4]]
	queryseq1[int(var.site.position)-1] = var.alleles[var.genotypes[2]]
	queryseq2[int(var.site.position)-1] = var.alleles[var.genotypes[3]]
num_true_false_diffs = sum( trueseq[j] != falseseq[j] for j in range(len(trueseq)))

# Now get "reads" (random substrings) from the query sequences.
readset = []
readstarts = []
for i in range(0,num_reads):
	startpos = random.randint(0,seqlength-read_length-1)
	strand = random.randint(0,2)
	if strand == 0:
		read = queryseq1[startpos:startpos+read_length]
	else:
		read = queryseq2[startpos:startpos+read_length]
	readset.append(read)
	readstarts.append(startpos)

# How many mismatches does each read have with the true and false sequences?
lm_correct = 0; lm_incorrect = 0; lm_nodec = 0
em_correct = 0; em_incorrect = 0; em_nodec = 0
for i in range(0,num_reads):
	query = readset[i]
	true = trueseq[readstarts[i]:readstarts[i]+read_length]
	false = falseseq[readstarts[i]:readstarts[i]+read_length]
	truemm = sum( query[j] != true[j] for j in range(len(query)))
	falsemm = sum( query[j] != false[j] for j in range(len(query)))
	if truemm < falsemm:
		lm_correct = lm_correct + 1
	elif falsemm < truemm:
		lm_incorrect = lm_incorrect + 1
	else:
		lm_nodec = lm_nodec + 1
	if truemm == 0 and falsemm > 0:
		em_correct = em_correct + 1
	elif truemm > 0 and falsemm == 0:
		em_incorrect = em_incorrect + 1
	else:
		em_nodec = em_nodec + 1

outputstring = "L " + str(seqlength) + " k " + str(read_length) + " N " + str(pop_size) + " R " + str(num_reads) + " mu " + str(mutation_rate) + " Tqt " + str(Tqt) + " Ttf " + str(Ttf) + " gen " + str(generation_length) + " qA " + str(query_age) + " tA " + str(true_age) + " fA " + str(false_age) + " r " + str(recomb_rate)
em_result_string = "Exact match: Correct: " + str(em_correct) + " Incorrect: " + str(em_incorrect) + " No_assignment: " + str(em_nodec)
lm_result_string = "Least mismatch: Correct: " + str(lm_correct) + " Incorrect: " + str(lm_incorrect) + " No_assignment: " + str(lm_nodec)
diff_string = "Your reference sequences differed at " + str(num_true_false_diffs) + " sites."

# Print results
if print_results:
	print(outputstring)
	print(diff_string)
	print(em_result_string)
	print(lm_result_string)
	
if write_files:
	f = open(file_path,"a")
	f.write(outputstring); f.write("\n")
	f.write(diff_string); f.write("\n")
	f.write(em_result_string); f.write("\n")
	f.write(lm_result_string); f.write("\n")
	f.close()

# Output true and false fasta files (for one strand each, as chosen above) and a fasta file containing all the query reads. 
if create_fastas:
	ft = open("true.fa","w")
	ft.write(">True"); ft.write("\n")
	ft.write(''.join(trueseq)); ft.close()

	ff = open("false.fa","w")
	ff.write(">False"); ff.write("\n")
	ff.write(''.join(falseseq)); ff.close()

	fq = open("query.fa","w")
	for i in range(0,len(readset)):
		fq.write(">Query"+str(i+1)); fq.write("\n")
		fq.write(''.join(readset[i])); fq.write("\n")
	fq.close()






