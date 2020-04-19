import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import ttk
import matplotlib.pyplot as plt
from parsers.fasta_parser import FastaParser,InvalidSequenceException,InvalidCharsInSequenceException
from data_structures.fasta import FastaSequence
from api_connector.connector import APIConnector
from dot_plot_algorithm.dot_plot_impl import  DotPlot
from gui.global_alignment_window import GlobalAlignmentWindow
from gui.local_alignment_window import LocalAlignmentWindow
import re
from textwrap import wrap


class AppPanel(tk.Tk):


	def __init__(self,height=800,width=600,*args,**kwargs):
		tk.Tk.__init__(self,*args,**kwargs)
		self.height = height
		self.width = width

		self.resizable(False,False)

		self.matrix = None

		self.sequences = []
		self.top_windows = []

		self.font_type = ('Tahoma',8,'bold')
		self.BG_COLOR = "#212121"
		self.LIME = "#76FF03"

		canvas = tk.Canvas(self,height=self.height,width=self.width)
		canvas.pack()
		canvas.configure(background=self.BG_COLOR)

		#ComboBox that is associated with choosing data aquiring method
		self.combo_box_frame = tk.Frame(self,bg=self.BG_COLOR)
		self.combo_box_frame.place(relx=0.5, rely=0.05, relwidth=0.3, relheight=0.1, anchor='n')
		self.comboBox = ttk.Combobox(self.combo_box_frame, values=['File', 'Api Request', 'Custom'])
		self.comboBox.place(relwidth=1,relheight=0.5)
		self.comboBox.current(1)
		self.comboBox.bind("<<ComboboxSelected>>",self.cb_selection)

		self.frame = tk.Frame(self,bd=5,bg=self.BG_COLOR)
		self.frame.place(relx=0.5,rely=0.10,relwidth=0.75,relheight=0.85,anchor='n')

		#Widgets associated with api request url 
		self.url_widgets = []
		self.url_label_3 = tk.Label(self.frame, font=self.font_type, bg=self.LIME)
		self.url_label_3.place(relx=0.5,rely=0.05,relwidth=0.5,relheight=0.05,anchor='n')
		self.url_label_3['text'] = 'DB type:'
		self.url_entry_3 = tk.Entry(self.frame,font=self.font_type,bd=1)
		self.url_entry_3.place(relx=0.5,rely=0.1,relwidth=0.5,relheight=0.075,anchor='n')
		self.url_label_1 = tk.Label(self.frame, font=self.font_type,bg=self.LIME)
		self.url_label_1.place(relx=0.5, rely=0.2, relwidth=0.5, relheight=0.05, anchor='n')
		self.url_label_1['text'] = '1st sequence ID:'
		self.url_entry_1 = tk.Entry(self.frame, font=self.font_type, bd=1)
		self.url_entry_1.place(relx=0.5, rely=0.25, relwidth=0.5, relheight=0.1, anchor='n')
		
		self.url_label_2 = tk.Label(self.frame, font=self.font_type,bg=self.LIME)
		self.url_label_2.place(relx=0.5, rely=0.4, relwidth=0.5, relheight=0.05, anchor='n')
		self.url_label_2['text'] = '2nd sequence ID:'
		self.url_entry_2 = tk.Entry(self.frame, font=self.font_type, bd=1)
		self.url_entry_2.place(relx=0.5, rely=0.45, relwidth=0.5, relheight=0.1, anchor='n')

		self.url_get_btn = tk.Button(self.frame,text='Request',font=self.font_type,bg=self.LIME,command=lambda:self.make_request())
		self.url_get_btn.place(relx=0.5,rely=0.6,relwidth=0.4,relheight=0.1,anchor='n')

		self.url_widgets.extend([self.url_label_1,self.url_label_2,self.url_entry_1,self.url_entry_2,self.url_get_btn,self.url_entry_3,self.url_label_3])


		#Button that opens filechooser dialog
		self.filechooser_btn = tk.Button(self.frame,text='Read File',bg=self.LIME,font=self.font_type,command=lambda:self.read_file())
		# self.filechooser_btn.place(relx=0.5,rely=0.3,relwidth=0.15,relheight=0.4,anchor='n')

		#Custom sequences entries
		self.first_seq_label = tk.Label(self.frame,font=self.font_type,bg=self.LIME)
		self.first_seq_label['text'] = 'First Sequence'
		self.first_seq_entry = tk.Entry(self.frame,font=self.font_type)
		self.second_seq_label = tk.Label(self.frame,font=self.font_type,bg=self.LIME)
		self.second_seq_label['text'] = 'Second Sequence'
		self.second_seq_entry = tk.Entry(self.frame,font=self.font_type)
		self.custom_seq_btn = tk.Button(self.frame,text='Pass Sequences',bg=self.LIME,font=self.font_type,command=lambda:self.load_custom_seq())
		self.custom_seq_widgets = []
		self.custom_seq_widgets.extend([self.first_seq_label,self.first_seq_entry,self.second_seq_label,self.second_seq_entry,self.custom_seq_btn])


		#Plot Button 
		self.plot_btn_frame = tk.Frame(self,bg=self.BG_COLOR)
		self.plot_btn_frame.place(relx=0.5,rely=0.75,relwidth=0.90,relheight=0.2,anchor='n')
		self.plot_btn = tk.Button(self.plot_btn_frame,bg=self.LIME,font=self.font_type,text='Plot',command=lambda:self.plot())
		self.plot_btn.place(relx=0.5,rely=0.15,relwidth=0.2,relheight=0.2,anchor='n')
		self.var1 = tk.IntVar()
		self.windowed_check_box =tk.Checkbutton(self.plot_btn_frame, text="use window", variable=self.var1,onvalue=1,offvalue=0,command=lambda:self.show_window_widgets())
		self.windowed_check_box.place(relx=0.85,rely=0.4,relwidth=0.2,relheight=0.3,anchor='n')

		#Window size and threshold values widgets
		self.window_size_label = tk.Label(self.plot_btn_frame,font=self.font_type,bg=self.LIME)
		self.window_size_label['text'] = 'Window size'
		self.window_size_entry = tk.Entry(self.plot_btn_frame,font=self.font_type)

		self.threshold_label = tk.Label(self.plot_btn_frame,font=self.font_type,bg=self.LIME)
		self.threshold_label['text'] = 'Threshold'
		self.threshold_entry = tk.Entry(self.plot_btn_frame,font=self.font_type)

		self.window_widgets = []
		self.window_widgets.extend([self.window_size_label,self.window_size_entry,self.threshold_label,self.threshold_entry])

		#Alignment Buttons
		self.global_alignment_btn = tk.Button(self.plot_btn_frame,bg=self.LIME,font=self.font_type,text="Global Alignment",command = lambda:self.show_alignment_window('global'))
		self.global_alignment_btn.place(relx=0.5,rely=0.4,relwidth=0.3,relheight=0.2,anchor='n')


		self.local_alignment_btn = tk.Button(self.plot_btn_frame,bg=self.LIME,font=self.font_type,text="Local Alignment",command = lambda:self.show_alignment_window('local'))
		self.local_alignment_btn.place(relx=0.5,rely=0.65,relwidth=0.3,relheight=0.2,anchor='n')


	def show_window_widgets(self):
		if(self.var1.get()==1):
			self.window_size_label.place(relx=0.05, rely=0.3, relheight=0.1, relwidth=0.2)
			self.window_size_entry.place(relx=0.05, rely=0.4, relheight=0.2, relwidth=0.2)
			self.threshold_label.place(relx=0.05,rely=0.65,relheight=0.1,relwidth=0.2)
			self.threshold_entry.place(relx=0.05,rely=0.75,relheight=0.2,relwidth=0.2)
		else:
			for item in self.window_widgets:
				self._hide(item)



	def show_alignment_window(self,type):

		if not self.check_sequences():
			return

		for top_window in self.top_windows:
			top_window.destroy()
		if type == 'global':
			alignment_window = GlobalAlignmentWindow(self)
		elif type == 'local':
			alignment_window = LocalAlignmentWindow(self)
		self.top_windows.clear()
		self.top_windows.append(alignment_window)


	def check_sequences(self):
		'''
		Method that is used to check if provided sequences are valid.
		:return: boolean value of sequences valid status.
		'''
		if len(self.sequences) != 2:
			tk.messagebox.showerror("ERROR", "Requiered 2 sequences!")
			return False

		for sequence in self.sequences:
			if sequence.get_sequence() == '':
				tk.messagebox.showerror("ERROR", "Requiered 2 sequences!")
				return False
		return  True

	def plot(self):
		"""
		Method that makes dot plot for given sequences. If given sequences are invalid it would show popup window with exception.
		:return:
		"""
		plt.close('all')

		if not self.check_sequences():
			return

		dot_plot_manager = DotPlot(self.sequences)
		if self.matrix is None:
			self.matrix = dot_plot_manager.make_dot_plot()
		if self.var1.get() == 1:
			window_size = self.window_size_entry.get()
			window_threshold = self.threshold_entry.get()
			if window_size=='' or window_threshold=='':
				tk.messagebox.showerror('ERROR','Provide window size and threshold level!')
				return
			shorter_seq = self.sequences[0].get_sequence_length() if self.sequences[0].get_sequence_length() < self.sequences[1].get_sequence_length() else self.sequences[1].get_sequence_length()
			if int(window_size) > shorter_seq:
				tk.messagebox.showerror('ERROR','Choose smaller window size!')
				return
			if int(window_threshold) > int(window_size):
				tk.messagebox.showerror('ERROR', 'Threshold value cannot be bigger than window size!')
				return
			windowed = dot_plot_manager.run_window(self.matrix,K = int(window_size),S = int(window_threshold))
			fig2 = plt.figure()
			plt.imshow(windowed,cmap='binary')
			plt.title(f'Window size = {window_size} and threshold = {window_threshold}')
			plt.xlabel('\n'.join(wrap(self.sequences[0].get_sequence_name(),50)))
			plt.ylabel('\n'.join(wrap(self.sequences[1].get_sequence_name(), 50)))
		fig  = plt.figure()
		plt.title('Dot Plot')
		plt.imshow(self.matrix,cmap='binary')
		plt.xlabel('\n'.join(wrap(self.sequences[0].get_sequence_name(),50)))
		plt.ylabel('\n'.join(wrap(self.sequences[1].get_sequence_name(),50)))
		fig.tight_layout()
		plt.show()

	def read_file(self):
		"""
		Method that reads selected file content and prepocess it before passing it to parser object.Selected item must be in fasta format.
		:return:
		"""
		filename = askopenfilename(filetypes=[("FASTA files","*.fasta")])
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

		self.matrix = None

		self._hide(self.url_label_1)


	def make_request(self):
		"""
		Method that enables connection to API of  ncbi database and requests data for a given pair of id and db type.
		Default db type is nuccore
		"""
		api_conn = APIConnector()

		db_type = self.url_entry_3.get()
		if db_type == '':
			db_type = 'nuccore'
		url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={db_type}&id={self.url_entry_1.get()},{self.url_entry_2.get()}&rettype=fasta&retmode=text"

		self.matrix = None
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

		first_seq_set = set(list(first_seq_fasta.get_sequence()))
		second_seq_set = set(list(second_seq_fasta.get_sequence()))

		#Validate User Input - Check if it contains only A,C,T,G chars
		if len(first_seq_set) > 4 or len(second_seq_set) > 4:
			tk.messagebox.showerror("ERROR", "Allowed only chars: A,T,C,G !")
			return
		for char in first_seq_set:
			try:
				self.validate_char(char)
			except InvalidCharsInSequenceException:
				tk.messagebox.showerror("ERROR", "Allowed only chars: A,T,C,G !")
				return
		for char in second_seq_set:
			try:
				self.validate_char(char)
			except InvalidCharsInSequenceException:
				tk.messagebox.showerror("ERROR", "Allowed only chars: A,T,C,G !")
				return
		self.matrix = None
		self.sequences = [first_seq_fasta,second_seq_fasta]

	def validate_char(self,char):
		if char not in ['A','T','C','G']:
			raise InvalidCharsInSequenceException('Only A,T,C,G are allowed !')

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

			self.filechooser_btn.place(relx=0.5,rely=0.3,relwidth=0.5,relheight=0.1,anchor='n')

		elif(self.comboBox.current()==1):
			self._hide(self.filechooser_btn)
			
			for item in self.custom_seq_widgets:
				self._hide(item)

			self.url_label_1.place(relx=0.5, rely=0.2, relwidth=0.5, relheight=0.05, anchor='n')
			self.url_entry_1.place(relx=0.5, rely=0.25, relwidth=0.5, relheight=0.1, anchor='n')
			self.url_label_2.place(relx=0.5, rely=0.4, relwidth=0.5, relheight=0.05, anchor='n')
			self.url_entry_2.place(relx=0.5, rely=0.45, relwidth=0.5, relheight=0.1, anchor='n')
			self.url_get_btn.place(relx=0.5, rely=0.6, relwidth=0.4, relheight=0.1, anchor='n')
			self.url_label_3.place(relx=0.5,rely=0.05,relwidth=0.5,relheight=0.05,anchor='n')
			self.url_entry_3.place(relx=0.5,rely=0.1,relwidth=0.5,relheight=0.075,anchor='n')


		elif(self.comboBox.current()==2):

			for item in self.url_widgets:
				self._hide(item)

			self._hide(self.filechooser_btn)

			self.first_seq_entry.place(relx=0.5,rely=0.25,relwidth=0.9,relheight=0.1,anchor='n')
			self.first_seq_label.place(relx=0.5,rely=0.20,relwidth=0.9,relheight=0.05,anchor='n')

			self.second_seq_entry.place(relx=0.5,rely=0.45,relwidth=0.9,relheight=0.1,anchor='n')
			self.second_seq_label.place(relx=0.5,rely=0.4,relwidth=0.9,relheight=0.05,anchor='n')

			self.custom_seq_btn.place(relx=0.5,rely=0.6,relwidth=0.4,relheight=0.1,anchor='n')
