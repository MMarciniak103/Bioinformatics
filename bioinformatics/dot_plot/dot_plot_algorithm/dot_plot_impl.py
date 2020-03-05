from data_structures.fasta import FastaSequence
import numpy as np

class DotPlot():
	def __init__(self,sequences):
		super().__init__()
		for seq in sequences:
			assert isinstance(seq,FastaSequence) , 'Sequenced must be object of FastaSequence'

		self.sequences = sequences
		self._xlen= len(sequences[0].get_sequence())
		self._ylen = len(sequences[1].get_sequence())


		assert (self._xlen > 1 and self._ylen > 1) , 'Sequences must be longer than 1'

	def change_sequences(self,sequences):
		self.sequences = sequences

	def make_dot_plot(self):
		dot_plot = np.zeros((self._xlen,self._ylen))
		dot_plot2 =  [[self.compare_pair(self.sequences[0].get_sequence()[i],self.sequences[1].get_sequence()[j]) for i in range(self._xlen)] for j in range(self._ylen)]
		return dot_plot2

	def compare_pair(self,element1,element2):
		if element1 == element2:
			return 1
		return 0

	def run_window(self,matrix,K=4,S=1):
		matrix = np.array(matrix)
		matrix2 = np.zeros(np.shape(matrix))
		#Find out how many steps would filter make in a given direction
		row_iters = (np.shape(matrix)[0] - K) + 1
		col_iters = (np.shape(matrix)[1] - K) + 1
		#Check diagonal of your window and decide if its allowed to stay or be replaced with zeros.
		for i in range(int(row_iters)):
			for j in range(int(col_iters)):
				window = matrix[i:i+K,j:j+K]
				#Get elements inside diagonal
				diagonal = np.diagonal(window)
				#Sum them to check if there is more than a given threshold
				sum_of_elemnents = np.sum(diagonal)
				#If diagonal contains more elements than given threshold, we clear all elements that lay outside it
				if sum_of_elemnents >= K-S:
					matrix2[i:i+K,j:j+K] = np.zeros(np.shape(window)) + np.diagflat(diagonal)
		return matrix2

