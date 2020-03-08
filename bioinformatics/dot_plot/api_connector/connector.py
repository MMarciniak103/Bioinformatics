import requests
from parsers.fasta_parser import FastaParser,InvalidSequenceException
import tkinter as tk

class APIConnector:
	def __init__(self):
		self.content = []
		self.parser = FastaParser()

	def request(self,url):
		"""
		Method that makes request to API of ncbi database
		:param url: url of api request
		:return: processed content of result
		"""
		try:
			result = requests.get(url)
			result_list = result.content.decode().split('\n')
			for i,v in enumerate(result_list):
				if v == '':
					del result_list[i]
			self.content = self.parser.parse_data(result_list)
		except InvalidSequenceException:
			tk.messagebox.showerror("Error","Requiered 2 sequences!")
		except Exception:
			tk.messagebox.showerror("Error","Invalid URL")

		return self.content
