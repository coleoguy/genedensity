TODO: 
1) a clear and explicit list of the requirements for a genome to be analyzed
2) a clear and explicit list of the requirements for a contig in a genome to be included


for each genome we analyze at most the largest 2N contigs




A genome is analyzed if
1. the results of its analysis does not exist in results/vertebrates/individual_species_results
2. the chromosome number of the species is available in data/vertebrates/chromnums.csv

A contig is included if
1. it is among the 2n longest contigs in the genome
2. its size is at least 10,000,000 bp
3. it has at least 2 genes

In addition, the results of a genome analysis will only be recorded if the linear model of gene count vs contig size returns a p-value less than 0.05.


