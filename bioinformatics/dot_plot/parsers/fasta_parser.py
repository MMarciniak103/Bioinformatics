from data_structures.fasta import FastaSequence
from collections import OrderedDict
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

		#Iterate over string content and parse sequences
		#Save them temporaly in OrderedDict
		seq_id = None
		seq_string = ''
		seqs = OrderedDict()
		for line in content:
			if line.startswith('>'):
				if seq_id is not None:
					seqs[seq_id] = seq_string
				seq_id = line[1:].strip()
			else:
				seq_string +=line.strip()

		seqs[seq_id] = seq_string

		#Iterate over parsed sequences and map them into FastaSequence Objects
		for k,v in seqs.items():
			sequence = FastaSequence(k)
			sequence.set_sequence(re.sub(r"\s+", "", v, flags=re.UNICODE))
			self.sequences.append(sequence)

		#Check if there was exactly 2 sequences, if no then throw an exception
		if len(self.sequences) != 2:
			raise InvalidSequenceException("REQUIERED 2 SEQUENCES!")

		return self.sequences



class InvalidSequenceException(Exception):
	def __init__(self,message):
		super().__init__(message)