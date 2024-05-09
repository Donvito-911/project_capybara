from __future__ import division

from glob import glob
from cogent import LoadSeqs, DNA
from cogent.evolve.pairwise_distance import TN93Pair
from cogent.phylo import distance
from cogent.core.alignment import Alignment
import os, sys
import numpy as np

from Bio.Align import MultipleSeqAlignment 
from Bio.SeqRecord import SeqRecord 
from Bio.codonalign.codonalphabet import default_codon_table, default_codon_alphabet 
from Bio.codonalign.codonseq import _get_codon_list, CodonSeq, cal_dn_ds

try:
		dirseq=sys.argv[1]
		print os.path.split(dirseq)[-1] + "...OK!"
except:
		print "Please pass a directory with fasta sequences"

os.chdir(dirseq)
path=os.getcwd() + '/*'
sequences = glob(path)

file_1 = 'pairwise_dnds_Gpig-Capybara-Rat.out'
filepath1=os.path.join(dirseq.rsplit('/',1)[0],file_1)
f=open(filepath1,'w')

# Arreglo para guadar el nombre del cluster de ortologos
Cpor_seqs=[]
Cap_seqs=[]
Rat_seqs=[]

# Arreglo para guardar los valores de omega
CGP = [] #Cap-GuineaPig
GPR = [] #GuineaPig-Rat

just_nucs = lambda x: set(''.join(x)) <= set('ACGT')

#s = 0
for sequence in sequences:
		#ortologos_nombres.append(os.path.split(sequences[s])[-1].split('_')[0])
		aln1 = LoadSeqs(sequence, moltype=DNA)
		aln = aln1.filtered(just_nucs, motif_length=3)
		names = aln.getSeqNames()
		index_Cap = int([names.index(i) for i in names if 'Capybara' in i][0])
		index_Cpor = int([names.index(i) for i in names if 'Cavia' in i][0])
		index_Rat = int([names.index(i) for i in names if 'ENSRNOP' in i][0])
		CapCDS = aln.takeSeqs([names[index_Cap]])
		Cap_seqs.append(CapCDS.getSeqNames()[0].split('|')[1])
		CporCDS = aln.takeSeqs([names[index_Cpor]])
		Cpor_seqs.append(CporCDS.getSeqNames()[0].split('|')[1])
		RatCDS = aln.takeSeqs([names[index_Rat]])
		Rat_seqs.append(RatCDS.getSeqNames()[0])
		
		if ((len(aln) / len(aln1)) > 0.5): 
		#Solo analizar alineaientos que despues de eliminados los gaps, tengan una longitud mayor al 50% de la longitud original
			CapCDSstr = str(CapCDS.todict().values()[0])
			CapCDSstr = CodonSeq(CapCDSstr)
			CporCDSstr = str(CporCDS.todict().values()[0])
			CporCDSstr = CodonSeq(CporCDSstr)
			RatCDSstr = str(RatCDS.todict().values()[0])
			RatCDSstr = CodonSeq(RatCDSstr)
			try:
				if cal_dn_ds(CapCDSstr,CporCDSstr)[1] == 0.0:
					val1 = 0.0
				elif cal_dn_ds(CporCDSstr,RatCDSstr)[1] == 0.0:
					val2 = 0.0
				else:
					val1 = cal_dn_ds(CapCDSstr,CporCDSstr)[0]/cal_dn_ds(CapCDSstr,CporCDSstr)[1]
					val2 = cal_dn_ds(CporCDSstr,RatCDSstr)[0]/cal_dn_ds(CporCDSstr,RatCDSstr)[1]
				CGP.append(val1)
				GPR.append(val2)
			except:
				print 'Error with: ' + os.path.split(sequences[s])[-1].split('_')[0]
		else:
			CGP.append('NA')
			GPR.append('NA')
		#s = s + 1

Cpor_seqs = np.array(Cpor_seqs)
Cap_seqs = np.array(Cap_seqs)
Rat_seqs = np.array(Rat_seqs)
CGP = np.array(CGP)
GPR = np.array(GPR)
ort_dists=np.column_stack((Cap_seqs,Cpor_seqs,Rat_seqs,CGP,GPR))
ort_dists=ort_dists.tolist()

f.write('\n'.join(['\t'.join(['{:4}'.format(item) for item in row]) for row in ort_dists]))
