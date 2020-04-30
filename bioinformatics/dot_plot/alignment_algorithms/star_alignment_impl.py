from data_structures.fasta import FastaSequence
import numpy as np
from itertools import combinations
from alignment_algorithms.global_alignment_impl import GlobalAlignment
import copy
from collections import defaultdict
from functools import reduce

class StarAlignment:

    def __init__(self,sequences):


        for seq in sequences:
            assert isinstance(seq,FastaSequence), 'Sequences must be object of FastaSequence'

        self.sequences = sequences


    def get_element_combinations(self,elem,combinations):
        """
        Check which combinations contain desired element.
        :param elem: Element that should be part of found combinations.
        :param combinations:  List containing all possible combinations
        :return: list of combinations containing desired element
        """
        found = []
        for i,combinations in enumerate(combinations):
            if combinations[0] == elem or combinations[1] == elem:
                found.append(i)
        return found



    def predict_alignment(self,indent_cost,substitution_cost,match_cost):

        indices = list(range(len(self.sequences)))

        seen = []

        #Get combinations of given sequences
        seq_combinations = combinations(self.sequences,2)
        indices_combinations = combinations(indices,2)

        scores = []
        all_alignments = []
        all_scores = []

        #Find global alignment and its score for every combination
        for combination in seq_combinations:

            ga = GlobalAlignment(combination)
            alignments,score = ga.predict_alignment(indent_cost,indent_cost,substitution_cost,match_cost)

            all_alignments.append(alignments)
            all_scores.append(score)


        print(all_alignments)

        all_scores = np.array(all_scores)
        print(all_scores)