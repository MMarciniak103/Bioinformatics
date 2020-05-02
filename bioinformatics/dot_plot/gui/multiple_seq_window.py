import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import ttk
from gui.star_alignment_window import StarAlignmentWindow
from gui.app_panel import read_file


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

        self.frame = tk.Frame(self, bd=5, bg=self.BG_COLOR)
        self.frame.place(relx=0.5, rely=0.10, relwidth=0.75, relheight=0.85, anchor='n')

        #--------------------------------- READING FASTA FILES - BUTTON WIDGET -----------------------------------------
        self.filechooser_btn = tk.Button(self.frame, text='Read File', bg=self.LIME, font=self.font_type,command=lambda: read_file(self))
        self.filechooser_btn.place(relx=0.5, rely=0.1, relwidth=0.5, relheight=0.1, anchor='n')


        self.star_alignment = tk.Button(self.frame, text='STAR', bg=self.LIME, font=self.font_type,command = lambda :self.open_window('STAR'))
        self.star_alignment.place(relx=0.5,rely=0.25,relwidth=0.5,relheight=0.1,anchor='n')


    def check_sequences(self):
        '''
        Method that is used to check if provided sequences are valid.
        :return: boolean value of sequences valid status.
        '''
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

        self.windows.clear()
        self.windows.append(alignment_window)