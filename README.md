This is an attempt to make an automated workflow for phylogenetic tree generation that includes:

Blast search of a seed sequence within a fasta file of its protein family

Blast search of an outgroup sequence within a fasta file of its own family

Clusters the results to a chosen % identity

Combines both filkes together, and ensures consistent fasta formatting

Conducts a multiple sequence alignemtn with mafft

Creates a phylogenetic tree of the resulting alignment


The inputs for the code are accession numbers of the seed sequences and fasta files of the databases in which they will be blasted.

The resulting final code is very impressive and does include a generous amount of help from copilot, but I did attempt to code each step manually, and fell on copilot when i realised how my intended project was a bit more challenging than I first thought. spend many many hours trying to make trees exclusively using python and the bash shell, and the ramblings of an insane person lie in the files AutoTree trouble 1 and 2

The course was extremely insightful, and I can now confidently create automation scripts for when I need to investigate new phylogenies. Thanks very much for all your help :)
