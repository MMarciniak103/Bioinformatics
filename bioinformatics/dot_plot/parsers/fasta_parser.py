from data_structures.fasta import FastaSequence
import re 

class FastaParser:
	def __init__(self):
		self.sequences = []

	def parse_data(self,content):
		"""
		Method that process raw content of string and maps it into FastaSequence Object
		:param content: Sequence string in fasta format
		:return: list of FastaSequence Objects.
		"""
		new_seq = True
		start_reading = False
		# seq_name = ''
		seq_string = ''
		content_len = len(content)
		for i,line in enumerate(content):

			#Check if its already another one sequention and save the previous one
			if len(line) > 1:
				if line[0] == '>':
					new_seq = True
					if len(seq_string) != 0:
						sequence.set_sequence(re.sub(r"\s+", "", seq_string, flags=re.UNICODE))
						self.sequences.append(sequence)

			if new_seq:
				sequence = FastaSequence(line)
				seq_string = ''
				new_seq = False
				continue
				
			#Check if it is end of an content list		
			if i == content_len-1:
				sequence.set_sequence(re.sub(r"\s+", "", seq_string, flags=re.UNICODE))
				self.sequences.append(sequence)
			seq_string +=line

		if len(self.sequences) != 2:
			raise InvalidSequenceException("REQUIERED 2 SEQUENCES!")

		return self.sequences



class InvalidSequenceException(Exception):
	def __init__(self,message):
		super().__init__(message)