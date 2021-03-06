from data_structures.fasta import FastaSequence
import numpy as np

class GlobalAlignment:
    """
    Class that implements methods allowing to find global alignment of 2 sequences.
    It uses Needleman-Wunsch algorithm with distance as an metric.
    """
    def __init__(self,sequences):
        for seq in sequences:
            assert isinstance(seq, FastaSequence), 'Sequenced must be object of FastaSequence'

        self.sequences = sequences

    def _get_score(self,matrix, i, j, insertion, deletion, substitution, match):
        '''
        Find minimum cost of getting to position i,j from possible 3 states
        :param matrix: matrix containing alignment cost
        :param i: row position
        :param j: column position
        :param insertion: insertion cost value
        :param deletion: deletion cost value
        :param substitution: substitution cost value
        :return:cost of cell with position i,j
        '''
        diagonal_cost = substitution if self.sequences[1].get_sequence()[i - 1] != self.sequences[0].get_sequence()[j - 1] else match
        return np.min([matrix[i - 1][j - 1] +  diagonal_cost,
                       matrix[i - 1][j] + deletion,
                       matrix[i][j - 1] + insertion])


    def predict_alignment(self,insertion_cost=1,deletion_cost=1,substitution_cost=1,match_cost=0,return_matricies=False):
        '''
        Find Global Alignment of 2 sequences
        :param insertion_cost: cost of insertion
        :param deletion_cost: cost of deletion
        :param substitution_cost: cost of substitution
        :param return_matricies:  flag indicating if score matrix and path matrix should be returned
        :return: alignments, score of alignment, no. gaps, no. identity positions and visual traced path
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
                matrix[i][j] = self._get_score(matrix, i, j, insertion_cost, deletion_cost, substitution_cost,match_cost)

        #-------------------------------- TRACE BACK FOR OPTIMAL ALIGNMENT ------------------------------------------
        # This implementation is favoring substitution over insertion or deletion
        path_mark = 'x'
        if return_matricies:
            path_matrix = np.array([['' for i in range(m+1)] for j in range(n+1)]) #this matrix contains optimal alignment path
        i = n
        j = m
        if return_matricies:
            path_matrix[i][j]= path_mark
        aln1 = ""
        aln2 = ""
        aln3 = ""
        score = 0
        gaps = 0
        identity = 0
        while True:
            if i == 0 and j == 0:
                break
            else:
                #Check diagonal way
                diagonal_cost = substitution_cost if self.sequences[1].get_sequence()[i - 1] !=  self.sequences[0].get_sequence()[j - 1] else match_cost
                if ((i > 0 and j > 0 and (matrix[i][j] == matrix[i - 1][j - 1] + diagonal_cost))):
                    aln1 = self.sequences[0].get_sequence()[j - 1] + aln1
                    aln2 = self.sequences[1].get_sequence()[i - 1] + aln2
                    if (self.sequences[1].get_sequence()[i - 1] == self.sequences[0].get_sequence()[j - 1]):
                        identity += 1
                        aln3 = "|" + aln3
                        score += match_cost
                    else:
                        score += substitution_cost
                        aln3 = "*" + aln3
                    if return_matricies:
                        path_matrix[i-1][j-1]=path_mark
                    i -= 1
                    j -= 1

                else:
                    #Check horizontal way
                    if ((i > 0) and (matrix[i][j] == matrix[i - 1][j] + insertion_cost)):
                        aln1 = "-" + aln1
                        aln2 = self.sequences[1].get_sequence()[i - 1] + aln2
                        aln3 = " " + aln3
                        score += insertion_cost
                        if return_matricies:
                            path_matrix[i - 1][j] = path_mark
                        i -= 1
                        gaps += 1
                    #Check vertical way
                    else:
                        aln1 = self.sequences[0].get_sequence()[j - 1] + aln1
                        aln3 = " " + aln3
                        aln2 = "-" + aln2
                        score += deletion_cost
                        if return_matricies:
                            path_matrix[i][j-1] = path_mark
                        gaps += 1
                        j -= 1
        if return_matricies:
            return ([aln1,aln2,aln3],score,gaps,identity,matrix,path_matrix)


        return ([aln1,aln2,aln3],score)