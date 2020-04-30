import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import ttk
from gui.app_panel import AppPanel
from gui.multiple_seq_window import MultipleSeqWindow


class MainWindow(tk.Tk):

    def __init__(self,height=700,width=700,*args,**kwargs):
        tk.Tk.__init__(self,*args,**kwargs)
        self.height = height
        self.width = width

        self.font_type = ('Tahoma', 8, 'bold')
        self.BG_COLOR = "#212121"
        self.LIME = "#76FF03"

        self.windows = []

        self.resizable(False,False)

        canvas = tk.Canvas(self,height=self.height,width=self.width)
        canvas.pack()
        canvas.configure(background=self.BG_COLOR)

        #-------------------------- BUTTONS SECTION ------------------------------------------------
        self.buttons_frame = tk.Frame(self,bg=self.BG_COLOR)
        self.buttons_frame.place(relwidth=1,relheight=1)

        self.pair_seq_btn = tk.Button(self.buttons_frame,bg = self.LIME,font=self.font_type,text='Pair Of Sequences',command=lambda:self.open_window('PAIR'))
        self.pair_seq_btn.place(relx=0.5,rely=0.1,relwidth=0.3,relheight=0.1,anchor = 'n')

        self.multiple_seq_btn = tk.Button(self.buttons_frame,bg = self.LIME,font=self.font_type,text='Multiple Sequences',command=lambda:self.open_window('MULTI'))
        self.multiple_seq_btn.place(relx=0.5,rely=0.3,relwidth=0.3,relheight=0.1,anchor = 'n')



    def open_window(self,type):

        for window in self.windows:
            window.destroy()


        if type == 'PAIR':
            window = AppPanel(self)

        elif type == 'MULTI':
            window = MultipleSeqWindow(self)

        self.windows.clear()
        self.windows.append(window)
