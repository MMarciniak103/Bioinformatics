import numpy as np
from itertools import combinations
from alignment_algorithms.global_alignment_impl import GlobalAlignment
from data_structures.fasta import  FastaSequence
import copy
from collections import defaultdict
from functools import reduce
from alignment_algorithms.star_alignment_impl import StarAlignment
import matplotlib.pyplot as plt

seq1 = 'MATKA'
seq2 = "MOTINA"
seq3 = "MOTHER"
seq4 = "MUTTER"
seq5 = "MOR"

sequences = [seq1,seq2,seq3,seq4,seq5]


indices = list(range(len(sequences)))

indent_cost = 2
substitution_cost = 1
match_cost = 0


fasta_sequences = []

for i,seq in enumerate(sequences):
    fastaSeq = FastaSequence(f'x{i}')
    fastaSeq.set_sequence(seq)
    fasta_sequences.append(fastaSeq)


seq_hash = defaultdict()
for i in range(len(fasta_sequences)):
    seq_hash[fasta_sequences[i].get_sequence_name()] = i



seq_combinations = combinations(fasta_sequences,2)

print('FASTA SEQUENCES')
print(fasta_sequences)
for fasta in fasta_sequences:
    print(fasta.get_sequence() ,' ',fasta.get_sequence_name())



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
    scores_matrix[seq_hash[combination[0].get_sequence_name()]][seq_hash[combination[1].get_sequence_name()]] = score
    scores_matrix[seq_hash[combination[1].get_sequence_name()]][seq_hash[combination[0].get_sequence_name()]] = score

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
# C = [f'x{i}' for i in range(len(sequences))]
C = [seq.get_sequence_name() for seq in fasta_sequences]
print('C ', C)
connections = []

# scores_matrix=[
#     [0,3,4,4,5],
#     [3,0,3,4,4],
#     [4,3,0,2,3],
#     [4,4,2,0,4],
#     [5,4,3,4,0]
# ]
#
# scores_matrix = np.array(scores_matrix,dtype=float)

E = scores_matrix.copy()
np.fill_diagonal(E,np.Infinity)


def get_best_in_E(E,C,connections):
    min_id = np.argmin(E)
    min_value = np.min(E)
    indices = divmod(min_id,E.shape[1])
    #
    # print('MIN ID : ',min_id)
    # print('MIN VALUE: ',min_value)
    # print('INDICES = ',indices)

    connections.append([min_value/2,[C[i] for i in indices]])

    # print('CONNECTIONS ',connections)



    return connections,indices

def get_new_C_matrix(C,indices,connections):
    next_C = []
    for i,elem in enumerate(C):
        if i not in indices:
            next_C.append(elem)

    # for connection in connections:
    #     next_C.append(connection[1])
    next_C.append(connections[-1][1])

    return next_C


def flatten(S):
    if S == []:
        return S
    if isinstance(S[0], list):
        return flatten(S[0]) + flatten(S[1:])
    return S[:1] + flatten(S[1:])


def get_E_matrix(next_C):
    new_E = np.zeros((len(next_C), len(next_C)))
    for k,item in enumerate(next_C):
        for l,item2 in enumerate(next_C):

            if(item == item2):
                new_E[k][l] = np.Infinity
                continue

            row_elements = []
            columns_elements = []
            if isinstance(item,list):
                row_elements = flatten(item)
            else:
                row_elements = [item]
            if isinstance(item2,list):
                columns_elements = flatten(item2)
            else:
                columns_elements = [item2]

            score = 0
            for row_item in row_elements:
                for column_item in columns_elements:
                    i = seq_hash[row_item]
                    j = seq_hash[column_item]
                    score+= scores_matrix[i][j]/(len(row_elements)*len(columns_elements))

            new_E[k][l] = score
    return new_E




while (len(connections) != len(fasta_sequences)-1):
    connections, indices = get_best_in_E(E, C, connections)
    C = get_new_C_matrix(C,indices,connections)
    E = get_E_matrix(C)


print('END')
print('Connections ',connections)
print('C ',C)
print('E ',E)


### ----------------------------------------------- PLOT DENDROGRAM --------------------------------------------------

def mk_fork(x0,x1,y0,y1,new_level):
    points=[[x0,x0,x1,x1],[y0,new_level,new_level,y1]]
    connector=[(x0+x1)/2.,new_level]
    return (points),connector

levels = [connection[0] for connection in connections]
print('levels ',levels)

seen = []
locations = []
for connection in connections:
    # for elem in flatten(connection[1])[:2]:
    #     print(elem)
    # print(flatten(connection[1])[:2])
    # print(tuple(seq_hash[elem] for elem in flatten(connection[1])[:2]))
    first_elem = 0
    second_elem = 0
    if isinstance(connection[1][0],list):
        first_elem = connection[1][0][0]
    else:
        first_elem = connection[1][0]
    if isinstance(connection[1][1],list):
        second_elem = connection[1][1][0]
    else:
        second_elem = connection[1][1]
    # locations.append(tuple(seq_hash[elem] for elem in flatten(connection[1])[:2]))
    locations.append((seq_hash[first_elem],seq_hash[second_elem]))
    seen += flatten([elem for elem in flatten(connection[1])])

label_map = {i:{} for i in range(len(fasta_sequences))}
for i in range(len(fasta_sequences)):
    label_map[i] = {'label':fasta_sequences[i].get_sequence_name(),'xpos':i,'ypos':0}


print('locations ',locations)
print('label_map ',label_map)



fig,ax=plt.subplots()

for i,(new_level,(loc0,loc1)) in enumerate(zip(levels,locations)):

    print('step {0}:\t connecting ({1},{2}) at level {3}'.format(i, loc0, loc1, new_level ))

    x0,y0=label_map[loc0]['xpos'],label_map[loc0]['ypos']
    x1,y1=label_map[loc1]['xpos'],label_map[loc1]['ypos']

    print('\t points are: {0}:({2},{3}) and {1}:({4},{5})'.format(loc0,loc1,x0,y0,x1,y1))

    p,c=mk_fork(x0,x1,y0,y1,new_level)

    ax.plot(*p)
    ax.scatter(*c)
    print('\t connector is at:{0}'.format(c))


    label_map[loc0]['xpos']=c[0]
    label_map[loc0]['ypos']=c[1]
    label_map[loc0]['label']='{0}/{1}'.format(label_map[loc0]['label'],label_map[loc1]['label'])
    print('\t updating label_map[{0}]:{1}'.format(loc0,label_map[loc0]))

    ax.text(*c,label_map[loc0]['label'])

_xticks=np.arange(0,6,1)
_xticklabels=[*[seq.get_sequence_name() for seq in fasta_sequences]]

ax.set_xticks(_xticks)
ax.set_xticklabels(_xticklabels)

ax.set_ylim(0,1.05*np.max(levels))

plt.show()




# seen = []
# print('CONNECTION LOOP: ')
# text_print = ''
# for connection in (connections):
#     print(connection)
#     text_print+='\t'.join([elem for elem in flatten(connection[1]) if elem not in seen])
#     seen += [elem for elem in flatten(connection[1])]
#
# print(text_print)

#
#
# connections,indices = get_best_in_E(E,C,connections)
# print('Connections ',connections)
#
# next_C = get_new_C_matrix(C,indices,connections)
# print('next C: ',next_C)
#
# new_E = np.zeros((len(next_C),len(next_C)))




# print('TEST OF FUNC')
# testing_list = [[5,3,[8,7]],6]
# print(flatten(testing_list))
#
#
# new_E = get_E_matrix(next_C)
# print('new  E')
# print(new_E)
#
# connections,indices=get_best_in_E(new_E,next_C,connections)
# next_C = get_new_C_matrix(next_C,indices,connections)
# print('next C: ',next_C)
# new_E = get_E_matrix(next_C)
#
# print(new_E)
#
# print('-------------------')
# print('-------------------')
# connections,indices=get_best_in_E(new_E,next_C,connections)
# next_C = get_new_C_matrix(next_C,indices,connections)
# print('next C: ',next_C)
# new_E = get_E_matrix(next_C)
#
# print(new_E)
#
# print('-------------------')
# print('-------------------')
# connections,indices=get_best_in_E(new_E,next_C,connections)
# next_C = get_new_C_matrix(next_C,indices,connections)
# print('next C: ',next_C)
# new_E = get_E_matrix(next_C)
#
# print(new_E)
#
#
# print(connections)
