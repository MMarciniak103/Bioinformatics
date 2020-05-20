import numpy as np
from itertools import combinations
from alignment_algorithms.global_alignment_impl import GlobalAlignment
from data_structures.fasta import  FastaSequence
from collections import defaultdict
from alignment_algorithms.star_alignment_impl import StarAlignment
import matplotlib.pyplot as plt

class UPGMA:

    def __init__(self,sequences):
        for seq in sequences:
            assert isinstance(seq, FastaSequence), 'Sequences must be object of FastaSequence'

        self.sequences = sequences

    def format_groups_input(self,text):
        return f'{text}'.replace('[','(').replace(']',')')

    def _get_newick(self,HISTORY_C):
        """
        Returns history of connections formatted in newick's format.
        :param HISTORY_C: history of connections
        :return: history formatted in newick's format
        """
        newick_txt = ''
        for hist in reversed(HISTORY_C[-1]):
            # hist is saved in cost,groups format -> example: 1.0,[x1,x2]
            cost = hist[0]
            groups = hist[1]
            if(newick_txt == ''):
                newick_txt+= f'{self.format_groups_input(groups)}:{cost}'
            else:
                newick_txt = newick_txt.replace(f'{self.format_groups_input(groups)}',f'{self.format_groups_input(groups)}:{cost}')

        return newick_txt
    def _get_combinations_scores(self, combinations, indent, substitution, match, seq_hash):
        """
        Finds global alignments for given combinations of sequences. It also saves its scores into score matrix.
        :param combinations: generator containing all combinations of sequences
        :param indent: indent cost
        :param substitution: substitution cos
        :param match: match cost
        :param seq_hash: hash table containing sequences index in score matrix
        :return: scores matrix filled with alignments scores
        """

        scores_matrix = np.zeros((len(self.sequences),len(self.sequences)))
        all_scores = []

        for combination in combinations:
            ga = GlobalAlignment(combination)
            _,score = ga.predict_alignment(indent,indent,substitution,match)
            name1 = combination[0].get_sequence_name().split('.')[0]
            name2 = combination[1].get_sequence_name().split('.')[0]
            #Score matrix is symetrical !
            scores_matrix[seq_hash[name1]][seq_hash[name2]] = score
            scores_matrix[seq_hash[name2]][seq_hash[name1]] = score

            all_scores.append(score)

        return scores_matrix


    def get_best_in_E(self,E,C,connections):
        """
        Looks for the best alignment in current iteration. It connects sequences that has the lowest alignment cost and
        add them to connections list.
        :param E: Matrix containing alignment's scores for current sets
        :param C: List containing names of sequences
        :param connections: list containing connections in dendrogram
        :return: new connections list and indices of connected sequences
        """

        min_id = np.argmin(E)
        min_value = np.min(E)
        indices = divmod(min_id, E.shape[1])

        connections.append([min_value / 2, [C[i] for i in indices]])

        return connections, indices

    def _get_new_C_matrix(self, C, indices, connections,HISTORY_C):
        """
        Creates new name's list.
        :param C: old name's list
        :param indices: indices of connected sequences for current iteration
        :param connections: connections list
        :param HISTORY_C: history of all connections
        :return: new name's list.
        """
        next_C = []
        for i, elem in enumerate(C):
            if i not in indices:
                next_C.append(elem)


        HISTORY_C.append(connections)
        next_C.append(connections[-1][1])
        return next_C,HISTORY_C



    def _get_E_matrix(self, next_C, seq_hash, scores_matrix):
        """
        Creates new E matrix for current iteration
        :param next_C: list containing sequences names (grouped sequences represented as nested list)
        :param seq_hash: hash table containing information about sequence index in score matrix
        :param scores_matrix: original score matrix
        :return: new E matrix
        """
        new_E = np.zeros((len(next_C), len(next_C)))
        for k, item in enumerate(next_C):
            for l, item2 in enumerate(next_C):

                if (item == item2):
                    new_E[k][l] = np.Infinity
                    continue

                if isinstance(item, list):
                    row_elements = flatten(item)
                else:
                    row_elements = [item]
                if isinstance(item2, list):
                    columns_elements = flatten(item2)
                else:
                    columns_elements = [item2]

                score = 0
                for row_item in row_elements:
                    for column_item in columns_elements:
                        i = seq_hash[row_item]
                        j = seq_hash[column_item]
                        score += scores_matrix[i][j] / (len(row_elements) * len(columns_elements))

                new_E[k][l] = score
        return new_E

    def create_upgma(self,indent_cost,substitution_cost,match_cost):
        """
        Creates upgma dendrogram for given sequences. It utilizes recursion for that purpose.
        :param indent_cost: indent cost
        :param substitution_cost: substitution cost
        :param match_cost: match cost
        :return:
        """
        #Create hash table for sequences (key - seq name, value - its position in sequences list)
        seq_hash = defaultdict()
        for i in range(len(self.sequences)):
            seq_hash[self.sequences[i].get_sequence_name().split('.')[0]] = i

        seq_combinations = combinations(self.sequences,2)

        #Original score matrix
        scores_matrix = self._get_combinations_scores(seq_combinations, indent_cost, substitution_cost, match_cost, seq_hash)

        # C list containing names of sequences
        C = [seq.get_sequence_name().split('.')[0] for seq in self.sequences]
        # list containing connections in our dendrogram
        connections = []

        E = scores_matrix.copy() #Initialize E matrix as a copy of score matrix
        np.fill_diagonal(E,np.Infinity) # Diagnoal is filled with infinity, so that algorithm will ignore it when looking for minimum cost

        HISTORY_C = []

        while (len(connections) != len(self.sequences)-1):
            connections,indices = self.get_best_in_E(E,C,connections)
            C,HISTORY_C = self._get_new_C_matrix(C, indices, connections,HISTORY_C)
            E = self._get_E_matrix(C, seq_hash, scores_matrix)


        newick = self._get_newick(HISTORY_C)

        print('END')
        print('Connections ', connections)
        print('C ', C)
        print('E ', E)
        print('newick ',newick)

        return connections,seq_hash,newick

def flatten(S):
    """
    Helper functions that uses recursion to flatten nested list.
    :param S: nested list
    :return: flattened list
    """
    if S == []:
        return S
    if isinstance(S[0], list):
        return flatten(S[0]) + flatten(S[1:])
    return S[:1] + flatten(S[1:])