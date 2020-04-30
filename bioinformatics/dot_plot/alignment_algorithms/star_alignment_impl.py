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
        self.seen = []


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

    def transform_central(self,central,central_v2):
        """
        Method that is used to get new central sequence alignment by combining version from previous iteration with
        the new one. It is used to propagate new gaps
        :param central: central sequence from previous iteration
        :param central_v2: central sequence acquired in current iteration
        :return: new central sequence containing gaps form both versions.
        """
        gaps = []
        for i,mark in enumerate(central):
            if mark == '-':
                gaps.append(i)

        self.seen = copy.copy(gaps)

        for i,mark in enumerate(central_v2):
            if mark == '-':
                k = i
                for gap in self.seen:
                    if gap < i :
                        k+=1
                        self.seen.remove(gap)
                #don't repeat same gaps
                if k not in gaps:
                    gaps.append(k)

        inc = 0
        clean_seq = central.replace("-","")
        for gap in gaps:
            clean_seq = clean_seq[:(gap+inc)]+'-'+clean_seq[(gap+inc):]
            if gap != 0:
                inc +=1

        return clean_seq

    def propagate_gaps(self,before,after,seq):
        """
        Method used to propagate new gaps in previously aligned sequences
        :param before: old central sequence
        :param after:  new central sequence with propagated new gaps
        :param seq: sequence that is being filled with new gaps
        :return: sequence with propagated new gaps
        """
        i = 0
        for mark in after:
            if i >= len(before):
                seq = seq + mark
            elif mark != before[i]:
                if mark == '-':
                    seq = seq[:i] + mark + seq[i:]

            else:
                i +=1

        return seq

    def append_next_seq(self,central_x,comb_idx,all_scores,indices_combinations,all_alignments,multiple_alignments):
        """
        Method used to make next step in every iteration. It propagates new gaps in previously appended sequences, and
        appends new sequence with propagated gaps to alignments list.
        :param central_x: central sequence
        :param comb_idx: list of combinations idx
        :param all_scores:  list containing scores for every combination
        :param indices_combinations: list containing indices of sequences that make combinations
        :param all_alignments: list containing all global alignments
        :param multiple_alignments: list containing sequences that makes multiple alignment
        :return: multiple_alignments list extended with new sequence
        """
        min_score_id = np.argmin(all_scores) # Get alignment with best score for current iteration
        best_align_id = comb_idx[min_score_id] # id of this alignment in list containing all combinations
        best_align_indices = indices_combinations[best_align_id] # Indices of sequences that make best alignment

        IDX = np.argwhere(best_align_indices == central_x).item() # Get position of central sequence from best indices
        best_alignment = all_alignments[best_align_id]
        align1 = best_alignment[IDX] # Central Sequence
        align3 = best_alignment[2] # Markers
        align2 = best_alignment[0 if IDX == 1 else 1] # Second sequence

        # If list containing alignments is empty we need to populate it (It means it is first iteration)
        if len(multiple_alignments) == 0:
            multiple_alignments = [align1,align2]
        else:
            transformed = self.transform_central(multiple_alignments[0],align1) # Propagate gaps in central sequence

            #Propagate gaps through all previous alignments
            for i,align in enumerate(multiple_alignments[1:]):
                multiple_alignments[i+1] = self.propagate_gaps(multiple_alignments[0],transformed,align)

            newseq = self.propagate_gaps(align1,transformed,align2) #propgate gaps through new sequence (So it contains gaps from previous sequences in alignment)
            multiple_alignments.append(newseq) # Add new sequence to our alignments
            multiple_alignments[0] = transformed # Replace old central seq with the transformed one

        deleted = np.delete(all_scores,[min_score_id]) # Remove items that were used in current iteration
        deleted_comb_idx = np.delete(comb_idx,[min_score_id])

        return deleted,deleted_comb_idx,multiple_alignments


    def extend_sequences(self,multiple_alignments):
        """
        Helper function that is used to extend found sequences so they are all equal in length
        :param multiple_alignments: found multiple sequences alignment
        :return: transformed multiple_alignments
        """
        longest_seq = len(max(multiple_alignments,key=len))
        for i in range(len(multiple_alignments)):
            if len(multiple_alignments[i]) < longest_seq:
                multiple_alignments[i] = multiple_alignments[i] + '-' * (longest_seq - len(multiple_alignments[i]))

        return multiple_alignments


    def get_matches(self,symbol_set,match):
        """
        Get score for matches
        :param symbol_set: Dictionary containing symbols frequency in given column
        :param match: match cost
        :return: score
        """
        score = 0
        for key in symbol_set.keys():
            comb_number = len(list(combinations(range(1, symbol_set[key] + 1), 2)))
            if comb_number > 0:
                score += match * comb_number

        return score

    def get_substitutions(self,symbol_set, no_sequences,substitution):
        """
        Get score for substitutions
        :param symbol_set: Dictionary containing symbols frequency in given column
        :param no_sequences: number of sequences
        :param substitution:  substitution cost
        :return: score
        """
        score = 0
        # print(symbol_set)
        symbol_set = copy.copy(symbol_set)

        if symbol_set.__contains__('-'):
            symbol_set.pop('-')
        # print(symbol_set)
        # number_different_symbols = len(symbol_set.keys())

        no_pos_without_gaps = reduce(lambda x, y: x + y, symbol_set.values())

        for key in symbol_set.keys():
            if symbol_set[key] != no_sequences:
                score += symbol_set[key] * substitution * (no_pos_without_gaps - symbol_set[key])

        score = score / 2

        return score

    def get_gaps(self,symbol_set, no_sequences,deletion):
        """
        Calculate score for gaps
        :param symbol_set: Dictionary containing symbols frequency in given column
        :param no_sequences: number of sequences
        :param deletion: Indent cost
        :return: score
        """
        score = 0
        # print(symbol_set)
        no_gaps = symbol_set['-']
        score += deletion * (no_sequences - no_gaps) * no_gaps
        return score

    def calculate_alignment_score(self,alignment,match,substitution,gap):
        """
        Calculates score of found multiple sequences alignment
        :param alignment: list containing sequences
        :param match: match cost
        :param substitution: substitution cost
        :param gap: gap cost
        :return: score for this alignment
        """
        score = 0
        for i in range(len(alignment[0])):
            symbols = list(map(lambda x: x[i], alignment))

            symbol_set = defaultdict(int)
            for symbol in symbols:
                symbol_set[symbol] += 1

            old_score = score
            score +=    self.get_matches(symbol_set,match)
            # print('matches ',score)
            score += self.get_substitutions(symbol_set, len(alignment),substitution)
            # print('substitution ',score)
            score += self.get_gaps(symbol_set, len(alignment),gap)
            # print('gaps ',score)

            # print('iteration ',i,' ',score-old_score)
            # print('--------------------------')

        return score

    def predict_alignment(self,indent_cost,substitution_cost,match_cost):

        indices = list(range(len(self.sequences)))

        self.seen = []

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


        # print(all_alignments)

        all_scores = np.array(all_scores)
        # print(all_scores)
        saved_all_scores = copy.copy(all_scores)
        indices_combinations = np.array(list(indices_combinations)) # list containing combinations of all indices


        # Get overall score for each sequence
        elements_scores = []
        for elem in indices:
            # Get all combinations which contain given sequence
            comb_idx = self.get_element_combinations(elem,indices_combinations)
            # Get overall sequence score as a sum of its global alignments scores
            elem_score = np.sum(all_scores[comb_idx])
            elements_scores.append(elem_score)

        # print('ELEMENT SCORES: ', elements_scores)


        # Find candidate for central sequence
        central_x = np.argmin(elements_scores)
        # Find all combinations that contains central sequence
        comb_idx = self.get_element_combinations(central_x,indices_combinations)
        multiple_alignments = [] #Array that would be populated with new sequences in every iteration

        all_scores = all_scores[comb_idx]
        while len(multiple_alignments) < len(self.sequences):
            all_scores,comb_idx,multiple_alignments = self.append_next_seq(central_x,comb_idx,all_scores,indices_combinations,all_alignments,multiple_alignments)

        multiple_alignments = self.extend_sequences(multiple_alignments)

        symbols = list(map(lambda x:x[1],multiple_alignments))

        alignment_score = self.calculate_alignment_score(multiple_alignments,match_cost,substitution_cost,indent_cost)

        print('ALIGNMENT SCORE: ',alignment_score)
        print(np.sum(saved_all_scores))

        for align in multiple_alignments:
            print(align)