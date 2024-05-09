############################################################################

#This script reads a fasta alignment (in this case protein alignments) and identifies sites fixed
# in all rodents except for the Capybara (allowing one unknown 'X' base at each position). 

#The input file is a folder with fasta alignments.

#Use as follows: python extract_unique_residues_allow_oneX.py $PWD/FOLDER

#Read the README file to understand the structure of the output files

# (c) Santiago Herrera-Alvarez
#Herrera-Alvarez et al. (2018). How to make a giant rodent: Genomic basis and tradeoffs of gigantism in the capybara, the world's largest rodent. BioRxiv. doi: https://doi.org/10.1101/424606

############################################################################

# Start of script
from collections import Counter
from glob import glob
import numpy as np
from Bio import AlignIO
import os, sys

try:
        dirseq=sys.argv[1]
        print(os.path.split(dirseq)[-1] + "...OK!")
except:
        print("Please pass a directory with fasta sequences")

os.chdir(dirseq)
path=os.getcwd() + '/*'
sequences = glob(path)

file_1 = 'Unique_residues_index.out'
file_2 = 'Unique_residues_positions.out'

filepath1=os.path.join(dirseq.rsplit('/',1)[0],file_1) #File with unique residue concentration index values
f1=open(filepath1,'w')

filepath2=os.path.join(dirseq.rsplit('/',1)[0],file_2) #File with unique Capybara residues
f2=open(filepath2,'w')

#Function to trimm out gaps and stops (*) from the alignment
def prune_aln(aln):
    align_array = np.array([list(rec) for rec in aln], np.character)
    mask = []
    for column in range(alignment.get_alignment_length()):
        c = str(align_array[:,column])
        if ('-' in c) or ('*' in c):
            mask.append(False)
        else:
            mask.append(True)
    return np.compress(mask,align_array,axis=1)


clusters = []
index = [] # Number of unique Capybara residues normalized by pruned alignment length (after Keane et al. 2015 --> Bowhead whale genome)

for sequence in sequences:
    print('Alignment:    ' + sequence.split('/')[-1].split('_')[0])
    f2.write('\n' + 'Alignment:    ' + sequence.split('/')[-1].split('_')[0] + '\n')
    clusters.append(sequence.split('/')[-1].split('_')[0])
    alignment = AlignIO.read(sequence, 'fasta')
    align_array_trim = prune_aln(alignment)
    single = {}
    #For each fasta file, the script prints the name of the file and creates an alignment (using numpy array).
    
    if ((align_array_trim.shape[1] / alignment.get_alignment_length()) > 0.5): 
    #Only analyze trimmed alignments that have a length >50% of original length.
        seq_ids=[]
        for record in alignment:
            seq_ids.append(record.id)
    
        for column in range(align_array_trim.shape[1]):
            counter = Counter(align_array_trim[:,column])
            main_letter = counter.most_common(1)[0][0]
        #For each column (position in the alignment), 'main_letter' is the most frequent letter, if several letters have the same frequency it will be one of them.
    
            if (len(counter) == 2) or (counter[b'X'] == 1 and (len(counter) - counter[b'X'] == 2)):
                # single polymorphism
                #"(...) a maximum of one unknown residue per position was allowed in species other than the bowhead". (Keane et al. 2015)
                l = list(counter)
                l = [x.decode('UTF8') for x in l]
                if ('X' in l) and (len(counter) > 2):
                    del counter[b'X']
                letter = counter.most_common(2)[1][0]
                for pos in np.where(align_array_trim[:,column] == letter)[0]:
                    change = main_letter.decode('UTF-8') + str(column + 1) + letter.decode('UTF-8')
                    single.setdefault(change, [])
                    single[change].append(seq_ids[int(pos)].split('|')[0])
            #If there is only two different letters in that position it means that there is a single polyorphism in the position.
        
        print("   " + "Single polymorphisms:" + "  " + str(len(single)))
        f2.write("   " + "Alignment length:" + "  " + str(align_array_trim.shape[1]) + '\n')
        f2.write("   " + "Single polymorphisms:" + "  " + str(len(single)) + '\n')
        unique_changes_count = 0
        
        #Extract the aminoacid changes relative to 'main_letter' that occur ONLY in a particular lineage (in this case the capybara)
        #based on the sequence ID.
        for change in single:
            number_spp = len(single.get(change))
            if (number_spp == 1) and (single.get(change)[0] == 'Capybara'):
                unique_changes_count += 1
                f2.write("     " + change + ": " + str(single[change])[1:-1] + '\n')
                f2.write('')
        residue_concentration = unique_changes_count / align_array_trim.shape[1]
        index.append(residue_concentration)
    else:
        index.append('NA')
        f2.write('   ' + 'IGNORED: too many gaps!' + '\n')

clusters = np.array(clusters)
index = np.array(index)
data = np.column_stack((clusters,index))
data=data.tolist()

f1.write('\n'.join(['\t'.join(['{:4}'.format(item) for item in row]) for row in data]))

f1.close()
f2.close()

# End of script


