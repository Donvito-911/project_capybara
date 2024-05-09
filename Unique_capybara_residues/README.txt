How to interpret the output of the 'extract_unique_residues_allow_oneX.py' script:

1) The file 'Unique_residues_index_PRANK.out' contains a table with the computed residue-concentration index:
number of unique-capybara residues / alignment length

Structure:
Cluster name (Orthofinder)
Mouse Ensembl Protein ID
Guinea Pig Ensembl Protein ID
Index value	

Ej.
OG0013987	ENSMUSP00000035437.8	ENSCPOP00000006950.2	0.001876172607879925

--> 'NA' in the index column means that length of the pruned alignment was <50% of the original alignment and was
excluded for further analyses.

2) The file 'Unique_residues_positions_PRANK.out' contains a detailed description of the positions of each
unique-capybara residue.

Structure:
Cluster name (Orthofinder)
Length of the pruned alignment
Number of single substitutions found in the alignment
Specific substitutions that occurred in the capybara lineage

Ej.
Alignment:    OG0013987
   Alignment length:  1066
   Single polymorphisms:  299
     S963A: 'Capybara'
     L142P: 'Capybara'
     
-->S963A = A change of a Serine (S) for an Alanine (A) at position 963
--> In this case the index would be 2/1066 = 0.0018 (see first example)