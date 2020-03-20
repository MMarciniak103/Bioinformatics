from data_structures.fasta import FastaSequence
import numpy as np

class GlobalAlignment:
    def __init__(self,sequences):
        for seq in sequences:
            assert isinstance(seq, FastaSequence), 'Sequenced must be object of FastaSequence'

        self.sequences = sequences

    def _get_score(self,matrix, i, j, insertion, deletion, substitution):
        '''
        Find minimum cost of getting to position i,j from possible 3 states
        :param matrix: matrix containing alignment cost
        :param i: row position
        :param j: column position
        :param insertion: insertion cost value
        :param deletion: deletion cost value
        :param substitution: substitution cost value
        :return:
        '''
        return np.min([matrix[i - 1][j - 1] + substitution * (self.sequences[1].get_sequence()[i - 1] != self.sequences[0].get_sequence()[j - 1]),
                       matrix[i - 1][j] + deletion,
                       matrix[i][j - 1] + insertion])


    def predict_alignment(self,insertion_cost=1,deletion_cost=1,substitution_cost=1):
        '''
        Find Global Alignment of 2 sequences
        :param insertion_cost: cost of insertion
        :param deletion_cost: cost of deletion
        :param substitution_cost: cost of substitution
        :return: alignments, score of alignment, no. gaps and no. identity positions
        '''
        n = len(self.sequences[1].get_sequence())
        m = len(self.sequences[0].get_sequence())

        #create matrix containing alignments scores and populate it with zeros
        matrix = np.zeros((n + 1, m + 1))

        #populate first column with cost values
        for i in range(1, matrix.shape[0]):
            matrix[i][0] = matrix[i - 1][0] + insertion_cost

        #populate first row with cost values
        for j in range(1, matrix.shape[1]):
            matrix[0][j] = matrix[0][j - 1] + deletion_cost

        #Find cost of every possible position -> dynamic programming
        for i in range(1, matrix.shape[0]):
            for j in range(1, matrix.shape[1]):
                matrix[i][j] = self._get_score(matrix, i, j, insertion_cost, deletion_cost, substitution_cost)

        #-------------------------------- TRACE BACK FOR OPTIMAL ALIGNMENT ------------------------------------------
        # This implementation is favoring substitution over insertion or deletion
        i = n
        j = m
        aln1 = ""
        aln2 = ""
        aln3 = ""
        score = 0
        mismatch =0
        gaps = 0
        identity = 0
        while True:
            if i == 0 and j == 0:
                break
            else:
                #Check diagonal way
                if ((i > 0 and j > 0 and (matrix[i][j] == matrix[i - 1][j - 1] + substitution_cost * (self.sequences[1].get_sequence()[i - 1] != self.sequences[0].get_sequence()[j - 1])))):
                    aln1 = self.sequences[0].get_sequence()[j - 1] + aln1
                    aln2 = self.sequences[1].get_sequence()[i - 1] + aln2
                    if (self.sequences[1].get_sequence()[i - 1] == self.sequences[0].get_sequence()[j - 1]):
                        identity += 1
                        aln3 = "|" + aln3
                    else:
                        score += substitution_cost
                        aln3 = "*" + aln3
                    i -= 1
                    j -= 1

                else:
                    #Check horizontal way
                    if ((i > 0) and (matrix[i][j] == matrix[i - 1][j] + insertion_cost)):
                        aln1 = "-" + aln1
                        aln2 = self.sequences[1].get_sequence()[i - 1] + aln2
                        aln3 = " " + aln3
                        score += insertion_cost
                        i -= 1
                        gaps += 1
                    #Check vertical way
                    else:
                        aln1 = self.sequences[0].get_sequence()[j - 1] + aln1
                        aln3 = " " + aln3
                        aln2 = "-" + aln2
                        score += deletion_cost
                        gaps += 1
                        j -= 1

        return ([aln1,aln2,aln3],score,gaps,identity)