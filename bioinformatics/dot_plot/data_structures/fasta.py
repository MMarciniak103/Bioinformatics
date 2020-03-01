class FastaSequence:
	def __init__(self,name):
		self.name = name
		self.sequence = []

	def get_sequence(self):
		return self.sequence[0]
		
	def set_sequence(self,sequence):
		self.sequence.append(sequence)

	def get_sequence_length(self):
		return len(self.sequence)

	def __str__(self):
		return f'{self.name}-----{self.sequence}'