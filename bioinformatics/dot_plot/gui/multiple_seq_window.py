import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import ttk

from api_connector.connector import APIConnector
from gui.star_alignment_window import StarAlignmentWindow
from gui.app_panel import read_file
from gui.upgma_window import UpgmaWindow


class MultipleSeqWindow(tk.Toplevel):
    def __init__(self, master, *args, **kwargs):
        tk.Toplevel.__init__(self, master, *args, **kwargs)
        self.height = 700
        self.width = 700

        self.resizable(False,False)

        self.windows = []

        self.matrix = None
        self.sequences = []


        self.font_type = master.font_type
        self.BG_COLOR = master.BG_COLOR
        self.LIME = master.LIME


        canvas = tk.Canvas(self,height=self.height,width=self.width)
        canvas.pack()
        canvas.configure(background=master.BG_COLOR)

        self.combo_box_frame = tk.Frame(self, bg=self.BG_COLOR)
        self.combo_box_frame.place(relx=0.5, rely=0.05, relwidth=0.3, relheight=0.1, anchor='n')
        self.comboBox = ttk.Combobox(self.combo_box_frame, values=['File', 'Api Request'])
        self.comboBox.place(relwidth=1, relheight=0.5)
        self.comboBox.current(0)
        self.comboBox.bind("<<ComboboxSelected>>", self.cb_selection)

        self.frame = tk.Frame(self, bd=5, bg=self.BG_COLOR)
        self.frame.place(relx=0.5, rely=0.15, relwidth=0.75, relheight=0.85, anchor='n')

        #--------------------------------- READING FASTA FILES - BUTTON WIDGET -----------------------------------------
        self.filechooser_btn = tk.Button(self.frame, text='Read File', bg=self.LIME, font=self.font_type,command=lambda: read_file(self))
        self.filechooser_btn.place(relx=0.5, rely=0.1, relwidth=0.5, relheight=0.1, anchor='n')

        self.url_label = tk.Label(self.frame, text='URL KEYS', bg=self.LIME, font=self.font_type)
        self.url_entry = tk.Entry(self.frame,font = self.font_type)
        self.db_type_label =  tk.Label(self.frame, text='DB TYPE', bg=self.LIME, font=self.font_type)
        self.db_entry = tk.Entry(self.frame,font = self.font_type)
        self.url_widgets = [self.url_entry,self.url_label,self.db_entry,self.db_type_label]

        self.star_alignment = tk.Button(self.frame, text='STAR', bg=self.LIME, font=self.font_type,command = lambda :self.open_window('STAR'))
        self.star_alignment.place(relx=0.5,rely=0.4,relwidth=0.5,relheight=0.1,anchor='n')

        self.upgma = tk.Button(self.frame,text='UPGMA',bg=self.LIME,font=self.font_type,command=lambda:self.open_window('UPGMA'))
        self.upgma.place(relx=0.5,rely=0.6,relwidth=0.5,relheight=0.1,anchor='n')

    def _hide(self, widget):
        widget.place_forget()


    def check_sequences(self):
        '''
        Method that is used to check if provided sequences are valid.
        :return: boolean value of sequences valid status.
        '''
        if (self.comboBox.current() == 1):
            api_conn = APIConnector()

            db_type = self.db_entry.get()
            if db_type == '':
                db_type = 'nuccore'

            seq_idx = ','.join([id for id in self.url_entry.get().split(';')])

            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={db_type}&id={seq_idx}&rettype=fasta&retmode=text"
            self.sequences = api_conn.request(url)


        if len(self.sequences) <= 1:
            tk.messagebox.showerror("ERROR", "You must provide sequences")
            return False

        for sequence in self.sequences:
            if sequence.get_sequence() == '':
                tk.messagebox.showerror("ERROR", "Provided sequences are invalid!")
                return False
        return True

    def open_window(self,mode):
        if not self.check_sequences():
            return

        for window in self.windows:
            window.destroy()

        if mode == 'STAR':
            alignment_window = StarAlignmentWindow(self)
        if mode == 'UPGMA':
            alignment_window = UpgmaWindow(self)

        self.windows.clear()
        self.windows.append(alignment_window)

    def cb_selection(self, event):
        """
        Method that is handling user interaction with combo box. It checks which method to aquire data is chosen by user.
        :param event:
        """
        if (self.comboBox.current() == 0):
            for item in self.url_widgets:
                self._hide(item)

            self.filechooser_btn.place(relx=0.5, rely=0.1, relwidth=0.5, relheight=0.1, anchor='n')

        elif (self.comboBox.current() == 1):
            self._hide(self.filechooser_btn)

            self.db_type_label.place(relx=0.5,rely=0.05,relwidth=0.5,relheight=0.05,anchor='n')
            self.db_entry.place(relx=0.5,rely=0.1,relwidth=0.5,relheight=0.1,anchor='n')
            self.url_label.place(relx=0.5, rely=0.2, relwidth=0.8, relheight=0.05, anchor='n')
            self.url_entry.place(relx=0.5, rely=0.25, relwidth=0.8, relheight=0.1, anchor='n')

