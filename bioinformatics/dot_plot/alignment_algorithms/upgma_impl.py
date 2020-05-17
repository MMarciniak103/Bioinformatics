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

    def get_combinations_scores(self,combinations,indent,substitution,match,seq_hash):
        """
        Finds global alignments for given combinations of sequences. It also saves its scores into score matrix.
        :param combinations: generator containing all combinations of sequences
        :param indent: indent cost
        :param substitution: substitution cos
        :param match: match cost
        :param seq_hash: hash table containing sequences index in score matrix
        :return: scores matrix filled with alignments scores
        """

        scores_matrix = np.zeros(len(self.sequences),len(self.sequences))
        all_scores = []

        for combination in combinations:
            ga = GlobalAlignment(combination)
            _,score = ga.predict_alignment(indent,substitution,match)
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

    def get_new_C_matrix(self,C, indices, connections):
        """
        Creates new name's list.
        :param C: old name's list
        :param indices: indices of connected sequences for current iteration
        :param connections: connections list
        :return: new name's list.
        """
        next_C = []
        for i, elem in enumerate(C):
            if i not in indices:
                next_C.append(elem)


        next_C.append(connections[-1][1])
        return next_C

    def flatten(self,S):
        """
        Helper functions that uses recursion to flatten nested list.
        :param S: nested list
        :return: flattened list
        """
        if S == []:
            return S
        if isinstance(S[0], list):
            return self.flatten(S[0]) + self.flatten(S[1:])
        return S[:1] + self.flatten(S[1:])

    def get_E_matrix(self,next_C,seq_hash,scores_matrix):
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

                row_elements = []
                columns_elements = []
                if isinstance(item, list):
                    row_elements = self.flatten(item)
                else:
                    row_elements = [item]
                if isinstance(item2, list):
                    columns_elements = self.flatten(item2)
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
        for i in range(self.sequences):
            seq_hash[self.sequences[i].get_sequence_name().split('.')[0]] = i

        seq_combinations = combinations(self.sequences,2)

        #Original score matrix
        scores_matrix = self.get_combinations_scores(seq_combinations,indent_cost,substitution_cost,match_cost,seq_hash)

        # C list containing names of sequences
        C = [seq.get_sequence_name().split('.')[0] for seq in self.sequences]
        # list containing connections in our dendrogram
        connections = []

        E = scores_matrix.copy() #Initialize E matrix as a copy of score matrix
        np.fill_diagonal(E,np.Infinity) # Diagnoal is filled with infinity, so that algorithm will ignore it when looking for minimum cost


        while (len(connections) != len(self.sequences)-1):
            connections,indices = self.get_best_in_E(E,C,connections)
            C = self.get_new_C_matrix(C,indices,connections)
            E = self.get_E_matrix(C,seq_hash,scores_matrix)



        print('END')
        print('Connections ', connections)
        print('C ', C)
        print('E ', E)

