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
        self.path_mark = 'X'
        self.path_matrix = None

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
        aln3 = ''
        i, j = starting_pos
        findings_position = [(j,i)] #positions of alignment fragments in original sequences. ( i - 2nd seq, j - 1st one)
        ending_pos = None
        self.path_matrix[i][j] = self.path_mark
        #Follow path while cell values are bigger than 0.
        while matrix[i][j] > 0:
            if ((i > 0 and j > 0 and (
                    matrix[i][j] == np.max([matrix[i - 1][j - 1] + self.substitution_matrix[self.y[i - 1]][self.x[j - 1]], 0])))):
                aln1 = self.x[j - 1] + aln1
                aln2 = self.y[i - 1] + aln2
                if self.x[j-1] == self.y[i-1]:
                    aln3 = "|" + aln3
                else:
                    aln3 = "*"+aln3
                self.path_matrix[i-1][j-1] = self.path_mark
                i -= 1
                j -= 1
                ending_pos = (j,i)
            else:
                # Check horizontal way
                if ((i > 0) and (matrix[i][j] == np.max([matrix[i - 1][j] +  self.substitution_matrix['-'][self.y[i - 1]],0]))):
                    aln1 = "-" + aln1
                    aln2 = self.y[i - 1] + aln2
                    aln3 = " "+ aln3
                    self.path_matrix[i - 1][j] = self.path_mark
                    i -= 1
                    ending_pos = (j,i)
                # Check vertical way
                else:
                    aln1 = self.x[j - 1] + aln1
                    aln2 = "-" + aln2
                    aln3 = " " + aln3
                    self.path_matrix[i][j-1] = self.path_mark
                    j -= 1
                    ending_pos =  (j,i)
        findings_position.append(ending_pos)
        return [aln1, aln3, aln2],findings_position

    def _generate_score_matrix(self, n, m):
        """
        Generates score matrix for linear gap penalty - Smith-Waterman Algorithm
        :param n: rows number - length of 1 sequence
        :param m: columns number - length of 2 sequence
        :return: calculated score matrix ,list containing positions of biggest values, score of alignment
        """
        matrix = np.zeros((n + 1, m + 1))

        biggest_value = 0  # placeholder for biggest value seen
        biggest_value_pos = []  # list containing positions of cells with values equal to biggest value

        for i in range(1, matrix.shape[0]):
            for j in range(1, matrix.shape[1]):
                matrix[i][j] = self.get_score(matrix, i, j)

                # If current value is bigger than biggest_value, change reference
                if matrix[i][j] > biggest_value:
                    biggest_value = matrix[i][j]
                    # clear positions of previous biggest values
                    biggest_value_pos.clear()
                    biggest_value_pos.append((i, j))
                # If it is equal to biggest value append it list, so it can be referenced in backtracking
                elif matrix[i][j] == biggest_value:
                    biggest_value_pos.append((i, j))

        return (matrix, biggest_value_pos,biggest_value)

    def _affine_penalty(self,k,gap_open,gap_enlargement):
        return gap_open + gap_enlargement * k

    def _generate_score_matrix_affine(self,n,m,gap_enlargement):
        """
        Generates score matrix for affine gap penalty  - Gotoh Algorithm (Local)
        :param n: rows number - length of 1 sequence
        :param m: columns number  - length of 2 sequence
        :param gap_enlargement: gap enlargement cost value
        :return: 3 matrices used by algorithm to calculate final score matrix,list containing positions of biggest
        score, alignment score
        """
        matrix = np.zeros((n + 1, m + 1))
        P = np.zeros((n + 1, m + 1))
        Q = np.zeros((n + 1, m + 1))


        for i in range(1, P.shape[1]):
            P[0][i] = -np.inf

        for i in range(1, Q.shape[0]):
            Q[i][0] = -np.inf

        biggest_value = 0
        biggest_value_pos = []
        for i in range(1, matrix.shape[0]):
            for j in range(1, matrix.shape[1]):

                P_gap_open = self.substitution_matrix[self.y[i-1]]['-']

                P[i][j] = np.max([matrix[i - 1][j] + self._affine_penalty(1, P_gap_open,gap_enlargement),
                                  P[i - 1][j] + gap_enlargement])

                Q_gap_open =  self.substitution_matrix['-'][self.x[j-1]]
                Q[i][j] = np.max([matrix[i][j - 1] + self._affine_penalty(1,Q_gap_open,gap_enlargement),
                                  Q[i][j - 1] + gap_enlargement])

                matrix[i][j] = np.max([
                    matrix[i - 1][j - 1] + self.substitution_matrix[self.y[i - 1]][self.x[j - 1]],
                    P[i][j],
                    Q[i][j],
                    0
                ])

                if matrix[i][j] > biggest_value:
                    biggest_value = matrix[i][j]
                    biggest_value_pos.clear()
                    biggest_value_pos.append((i, j))
                elif matrix[i][j] == biggest_value:
                    biggest_value_pos.append((i, j))

        return (P,Q,matrix,biggest_value_pos,biggest_value)

    def predict_alignment(self,affine_gaping = 0,gap_enlargement=0):
        """
        Find possible local alignments of 2 sequences
        :return: list containing all alignments found by algorithm and score value of this alignment. It also returns
        score matrix and path_matrix containing alignments path.
        """
        n = len(self.y)
        m = len(self.x)
        self.path_matrix =  np.array([['' for i in range(m+1)] for j in range(n+1)]) ##this matrix contains optimal alignments paths

        if affine_gaping ==0:
            matrix,biggest_value_pos,score = self._generate_score_matrix(n,m)
        else:
            P,Q,matrix,biggest_value_pos,score = self._generate_score_matrix_affine(n,m,gap_enlargement)

        alignments = []
        findings_positions = []
        #traceback local alignment for every cell that has the biggest value
        for pos in biggest_value_pos:
            alignment, findings_pos =self.find_alignment(matrix, pos)
            alignments.append(alignment)
            findings_positions.append(findings_pos)

        return alignments,score,matrix,self.path_matrix,findings_positions