#Script to calculate GC content from a set of genes.
#Santiago Herrera-Alvarez (01-22-2019)

import os, sys
import numpy as np

from Bio.Align import MultipleSeqAlignment 
from Bio.SeqRecord import SeqRecord 
from Bio.SeqUtils import GC123
from Bio import SeqIO

fasta_file = "capybara_genes.fasta" # Input fasta file

gene_code=['Capybara_gene_code']
gc_content_tot=['GC_content_total']
gc_first=['GC_first_codon_pos']
gc_second=['GC_second_codon_pos']
gc_third=['GC_third_codon_pos']
gene_length=['Gene_length']
protein_length=['Prot_length']

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
for seq in fasta_sequences:
    gene_code.append(seq.id) #ID exists only for 'Record' objects
    gc_content_tot.append(GC123(seq)[0])
    gc_first.append(GC123(seq)[1])
    gc_second.append(GC123(seq)[2])
    gc_third.append(GC123(seq)[3])
    gene_length.append(len(seq))
    protein_length.append(len(seq)/3)


gene_code=np.array(gene_code)
gc_content_tot=np.array(gc_content_tot)
gc_first=np.array(gc_first)
gc_second=np.array(gc_second)
gc_third=np.array(gc_third)
gene_length=np.array(gene_length)
protein_length=np.array(protein_length)
tbl=np.column_stack((gene_code,gc_content_tot,gc_first,gc_second,gc_third,gene_length,protein_length))
tbl=tbl.tolist()

print('\n'.join(['\t'.join(['{:4}'.format(item) for item in row]) for row in tbl]))





