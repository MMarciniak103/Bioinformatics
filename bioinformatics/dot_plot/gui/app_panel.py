import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
from parsers.fasta_parser import FastaParser,InvalidSequenceException
from data_structures.fasta import FastaSequence
from api_connector.connector import APIConnector
import re 

class AppPanel(tk.Tk):


	def __init__(self,height=800,width=1000,*args,**kwargs):
		tk.Tk.__init__(self,*args,**kwargs)
		self.height = height
		self.width = width

		self.sequences = []

		font_type = ('Tahoma',8,'bold')

		canvas = tk.Canvas(self,height=self.height,width=self.width)
		canvas.pack()

		#ComboBox that is associated with choosing data aquiring method
		self.combo_box_frame = tk.Frame(self)
		self.combo_box_frame.place(relx=0.5, rely=0.05, relwidth=0.1, relheight=0.1, anchor='n')
		self.comboBox = ttk.Combobox(self.combo_box_frame, values=['File', 'Api Request', 'Custom'])
		self.comboBox.place(relwidth=1,relheight=0.5)
		self.comboBox.current(1)
		self.comboBox.bind("<<ComboboxSelected>>",self.cb_selection)

		self.frame = tk.Frame(self,bd=5)
		self.frame.place(relx=0.5,rely=0.15,relwidth=0.75,relheight=0.1,anchor='n')

		#Widget associated with api request url 
		self.url_widgets = []
		self.url_label_1 = tk.Label(self.frame, font=font_type)
		self.url_label_1.place(relx=0.2, rely=0.1, relwidth=0.3, relheight=0.3, anchor='n')
		self.url_label_1['text'] = '1st sequence ID:'
		self.url_entry_1 = tk.Entry(self.frame, font=font_type, bd=1)
		self.url_entry_1.place(relx=0.2, rely=0.4, relwidth=0.3, relheight=0.4, anchor='n')
		
		self.url_label_2 = tk.Label(self.frame, font=font_type)
		self.url_label_2.place(relx=0.45, rely=0.1, relwidth=0.3, relheight=0.3, anchor='n')
		self.url_label_2['text'] = '2nd sequence ID:'
		self.url_entry_2 = tk.Entry(self.frame, font=font_type, bd=1)
		self.url_entry_2.place(relx=0.45, rely=0.4, relwidth=0.3, relheight=0.4, anchor='n')

		self.url_get_btn = tk.Button(self.frame,text='Request',font=font_type,command=lambda:self.make_request())
		self.url_get_btn.place(relx=0.8,rely=0.4,relwidth=0.15,relheight=0.4,anchor='n')

		self.url_widgets.extend([self.url_label_1,self.url_label_2,self.url_entry_1,self.url_entry_2,self.url_get_btn])


		#Button that opens filechooser dialog
		self.filechooser_btn = tk.Button(self.frame,text='Read File',font=font_type,command=lambda:self.read_file())
		# self.filechooser_btn.place(relx=0.5,rely=0.3,relwidth=0.15,relheight=0.4,anchor='n')

		#Custom sequences entries
		self.first_seq_label = tk.Label(self.frame,font=font_type)
		self.first_seq_label['text'] = 'First Sequence'
		self.first_seq_entry = tk.Entry(self.frame,font=font_type)
		self.second_seq_label = tk.Label(self.frame,font=font_type)
		self.second_seq_label['text'] = 'Second Sequence'
		self.second_seq_entry = tk.Entry(self.frame,font=font_type)
		self.custom_seq_btn = tk.Button(self.frame,text='Pass Sequences',font=font_type,command=lambda:self.load_custom_seq())
		self.custom_seq_widgets = []
		self.custom_seq_widgets.extend([self.first_seq_label,self.first_seq_entry,self.second_seq_label,self.second_seq_entry,self.custom_seq_btn])

		#Chart Frame 
		self.chart_frame = tk.Frame(self,bd=5)
		self.chart_frame.place(relx=0.5,rely=0.3,relwidth=0.5,relheight=0.5,anchor='n')
		figure= plt.Figure(figsize=(20,20),dpi=100)
		ax = figure.add_subplot(111)
		self.chart_type = FigureCanvasTkAgg(figure,self.chart_frame)
		self.chart_type.get_tk_widget().pack()
		
		#Plot Button 
		self.plot_btn_frame = tk.Frame(self)
		self.plot_btn_frame.place(relx=0.5,rely=0.8,relwidth=0.75,relheight=0.1,anchor='n')
		self.plot_btn = tk.Button(self.plot_btn_frame,font=font_type,text='Plot',command=lambda:self.plot())
		self.plot_btn.place(relx=0.5,rely=0.7,relwidth=0.2,relheight=0.3,anchor='n')


	def plot(self):
		print(len(self.sequences))

		if len(self.sequences) !=2 :
			tk.messagebox.showerror("ERROR","Requiered 2 sequences!")

		for sequence in self.sequences:
			if sequence.get_sequence()  == '':
				tk.messagebox.showerror("ERROR","Requiered 2 sequences!")


	def read_file(self):
		"""
		Method that reads selected file content and prepocess it before passing it to parser object.Selected item must be in fasta format.
		:return:
		"""
		filename = askopenfilename()
		if filename:
			if filename[-5:] != 'fasta':
				tk.messagebox.showerror("ERROR","Requiered FASTA file!")
				return 
			with open(filename,'r') as f:
				content = f.readlines()
				parser = FastaParser()
				for i,line in enumerate(content):
					content[i] = line.strip()
					if content[i] == '':
						del content[i]
				try:
					self.sequences = parser.parse_data(content)
				except InvalidSequenceException:
					tk.messagebox.showerror("ERROR","Requiered 2 sequences in file!")

		self._hide(self.url_label_1)


	def make_request(self):
		"""
		Method that enables connection to API of  ncbi database and requests data for a given pair of id.
		"""
		api_conn = APIConnector()

		url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={self.url_entry_1.get()},{self.url_entry_2.get()}&rettype=fasta&retmode=text"

		self.sequences = api_conn.request(url)
		

	def _hide(self,widget):
		widget.place_forget()

	def load_custom_seq(self):
		'''
		Method that converts input data into fasta sequence and appends it to sequences list.
		'''
		first_seq = self.first_seq_entry.get().upper()
		second_seq = self.second_seq_entry.get().upper()

		first_seq_fasta = FastaSequence('custom first')
		second_seq_fasta = FastaSequence('custom second')

		first_seq_fasta.set_sequence(re.sub(r"\s+","",first_seq,flags=re.UNICODE))
		second_seq_fasta.set_sequence(re.sub(r"\s+","",second_seq,flags=re.UNICODE))

		self.sequences = [first_seq_fasta,second_seq_fasta]



	def cb_selection(self,event):
		"""
		Method that is handling user interaction with combo box. It checks which method to aquire data is chosen by user.
		:param event:
		"""
		if(self.comboBox.current()==0):
	
			for item in self.custom_seq_widgets:
				self._hide(item)

			for item in self.url_widgets:
				self._hide(item)

			self.filechooser_btn.place(relx=0.5,rely=0.3,relwidth=0.15,relheight=0.4,anchor='n')

		elif(self.comboBox.current()==1):
			self._hide(self.filechooser_btn)
			
			for item in self.custom_seq_widgets:
				self._hide(item)

			self.url_label_1.place(relx=0.2, rely=0.1, relwidth=0.3, relheight=0.3, anchor='n')
			self.url_entry_1.place(relx=0.2, rely=0.4, relwidth=0.3, relheight=0.4, anchor='n')
			self.url_label_2.place(relx=0.45, rely=0.1, relwidth=0.3, relheight=0.3, anchor='n')
			self.url_entry_2.place(relx=0.45, rely=0.4, relwidth=0.3, relheight=0.4, anchor='n')
			self.url_get_btn.place(relx=0.8,rely=0.4,relwidth=0.15,relheight=0.4,anchor='n')

			for item in self.url_widgets:
				print(item)

		elif(self.comboBox.current()==2):

			for item in self.url_widgets:
				self._hide(item)

			self._hide(self.filechooser_btn)

			self.first_seq_entry.place(relx=0.2,rely=0.4,relwidth=0.4,relheight=0.4,anchor='n')
			self.first_seq_label.place(relx=0.2,rely=0.1,relwidth=0.3,relheight=0.3,anchor='n')

			self.second_seq_entry.place(relx=0.8,rely=0.4,relwidth=0.4,relheight=0.4,anchor='n')
			self.second_seq_label.place(relx=0.8,rely=0.1,relwidth=0.3,relheight=0.3,anchor='n')

			self.custom_seq_btn.place(relx=0.5,rely=0.5,relwidth=0.2,relheight=0.4,anchor='n')