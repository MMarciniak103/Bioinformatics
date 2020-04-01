from data_structures.fasta import FastaSequence
import numpy as np

class LocalAlignment:
    """
     Class that implements methods allowing to find local alignment of 2 sequences.
     It uses Smith-Waterman algorithm with similarity as an metric.
     """
    def __init__(self,sequences,substitution_matrix):
        """
        Initialize LocalAlignment Object
        :param sequences: sequences that would be used to alignment
        :param substitution_matrix: substitution matrix with shape (6x6) in a form of:
            A  C  T  G  U  -
         A  x  x  x  x  x  x
         C  x  x  x  x  x  x
         T  x  x  x  x  x  x
         G  x  x  x  x  x  x
         U  x  x  x  x  x  x
         -  x  x  x  x  x  x
        """
        for seq in sequences:
            assert isinstance(seq, FastaSequence), 'Sequenced must be object of FastaSequence'

        self.sequences = sequences
        self.y = self.sequences[1].get_sequence()
        self.x = self.sequences[0].get_sequence()
        self.substitution_matrix = substitution_matrix

    def get_score(self,matrix, i, j):
        '''
        Find maximum score of getting to position i,j from possible 3 states
        :param matrix: matrix containing alignment cost
        :param i: row position
        :param j: column position
        :return:score of cell with position i,j
        '''
        return np.max([matrix[i - 1][j - 1] + self.substitution_matrix[self.y[i - 1]][self.x[j - 1]],
                       matrix[i - 1][j] + self.substitution_matrix['-'][self.y[i - 1]],
                       matrix[i][j - 1] + self.substitution_matrix['-'][self.x[j - 1]],
                       0])

    def find_alignment(self,matrix, starting_pos):
        """
        Find local alignment of 2 sequences starting from given position.
        It is done with Smith-Waterman algorithm
        :matrix: matrix containing score of each possible state.
        :param starting_pos: cell that is a starting point for algorithm
        :return: local alignment of 2 sequences for given starting point
        """
        aln1 = ''
        aln2 = ''
        i, j = starting_pos
        #Follow path while cell values are bigger than 0.
        while matrix[i][j] > 0:
            if ((i > 0 and j > 0 and (
                    matrix[i][j] == np.max([matrix[i - 1][j - 1] + self.substitution_matrix[self.y[i - 1]][self.x[j - 1]], 0])))):
                aln1 = self.x[j - 1] + aln1
                aln2 = self.y[i - 1] + aln2
                i -= 1
                j -= 1
            else:
                # Check horizontal way
                if ((i > 0) and (matrix[i][j] == matrix[i - 1][j] + np.max([0, self.substitution_matrix['-'][self.y[i - 1]]]))):
                    aln1 = "-" + aln1
                    aln2 = self.y[i - 1] + aln2
                    i -= 1
                # Check vertical way
                else:
                    aln1 = self.x[j - 1] + aln1
                    aln2 = "-" + aln2
                    j -= 1
        return [aln1, aln2]



    def predict_alignment(self):
        """
        Find possible local alignments of 2 sequences
        :return: list containing all alignments found by algorithm.
        """
        n = len(self.y)
        m = len(self.x)

        matrix = np.zeros((n + 1, m + 1))

        biggest_value = 0 #placeholder for biggest value seen
        biggest_value_pos = [] #list containing positions of cells with values equal to biggest value

        for i in range(1, matrix.shape[0]):
            for j in range(1, matrix.shape[1]):
                matrix[i][j] = self.get_score(matrix, i, j)

                #If current value is bigger than biggest_value, change reference
                if matrix[i][j] > biggest_value:
                    biggest_value = matrix[i][j]
                    #clear positions of previous biggest values
                    biggest_value_pos.clear()
                    biggest_value_pos.append((i, j))
                #If it is equal to biggest value append it list, so it can be referenced in backtracking
                elif matrix[i][j] == biggest_value:
                    biggest_value_pos.append((i, j))

        alignments = []
        #traceback local alignment for every cell that has the biggest value
        for pos in biggest_value_pos:
            alignments.append(self.find_alignment(matrix, pos))

        return alignments