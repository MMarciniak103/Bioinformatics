import numpy as np
from itertools import combinations

from alignment_algorithms.global_alignment_impl import GlobalAlignment
from data_structures.fasta import  FastaSequence
import copy
from collections import defaultdict
from functools import reduce
from alignment_algorithms.star_alignment_impl import StarAlignment

seq1 = 'MATKA'
seq2 = "MOTINA"
seq3 = "MOTHER"
seq4 = "MUTTER"
seq5 = "MOR"

sequences = [seq1,seq2,seq3,seq4,seq5]
seq_hash = defaultdict()
for i in range(len(sequences)):
    seq_hash[sequences[i]] = i


indices = list(range(len(sequences)))

indent_cost = 2
substitution_cost = 1
match_cost = 0


fasta_sequences = []

for seq in sequences:
    fastaSeq = FastaSequence('custom')
    fastaSeq.set_sequence(seq)
    fasta_sequences.append(fastaSeq)


seq_combinations = combinations(fasta_sequences,2)

print('FASTA SEQUENCES')
print(fasta_sequences)
for fasta in fasta_sequences:
    print(fasta.get_sequence())



# scores_matrix = defaultdict()
# for seq in sequences:
#     scores_matrix[seq] = {}
#     for seq2 in sequences:
#         scores_matrix[seq][seq2] = 0

scores_matrix = np.zeros((len(sequences),len(sequences)))


all_alignments = []
all_scores = []

for combination in seq_combinations:
    ga = GlobalAlignment(combination)
    alignment,score = ga.predict_alignment(indent_cost,indent_cost,substitution_cost,match_cost)
    # scores_matrix[combination[0].get_sequence()][combination[1].get_sequence()] = score
    # scores_matrix[combination[1].get_sequence()][combination[0].get_sequence()] = score
    scores_matrix[seq_hash[combination[0].get_sequence()]][seq_hash[combination[1].get_sequence()]] = score
    scores_matrix[seq_hash[combination[1].get_sequence()]][seq_hash[combination[0].get_sequence()]] = score

    all_alignments.append(alignment)
    all_scores.append(score)

print('-----------------------------------')
print('Alignments')
for i,alignment in enumerate(all_alignments):
    print(alignment,' ',all_scores[i])

print('---------------------------------')
print('SCORES MATRIX')
print(seq_hash)
for row in scores_matrix:
    print(row)

sa = StarAlignment(fasta_sequences)
alignments, msa_score, opt_score = sa.predict_alignment(float(indent_cost), float(substitution_cost), float(match_cost))

print('--------------------------')
for alignment in alignments:
    print(alignment)
print(msa_score)
print(opt_score)


print()
C = [f'x{i}' for i in range(len(sequences))]
print('C ', C)
connections = []


E = scores_matrix.copy()
np.fill_diagonal(E,np.Infinity)
print(E)
print('argmin ',np.argmin(E))
min_id = np.argmin(E)
min_value = np.min(E)
indices = divmod(min_id,E.shape[1])
print('take ',E.take(min_id))


connections.append((min_value/2,[C[i] for i in indices]))
print(connections)