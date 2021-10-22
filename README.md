Oct 22, 2021

Bianca De Sanctis, bdd28@cam.ac.uk

Two code files to go with the paper "Accuracy of supervised binning for metagenomic sequences".

There is an R script theoretical-calculation.R which runs the theoretical accuracy calculation. It contains a main function to calculate the theoretical accuracy, with the assignment method and the type of calculation ("correct","incorrect" or "no") as input. Other input options can be seen at the top of the script in the R function. 
The script also contains the code to run the function over a parameter range and make the main figure in the paper.

There is also a python script simulate.py which runs the simulation calculation, with a wide range of parameter input options. It requires msprime (see https://pypi.org/project/msprime/). To see the input options type 

python simulate.py -h

The script calculates results for both assignment methods, and prints the results to screen or writes them to a file. It can also write fasta files for the true and false sequences and query reads for downstream analysis, such as adding deamination or running existing binning software.
