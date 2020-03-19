import tkinter as tk
from tkinter import ttk

class AlignmentWindow(tk.Toplevel):
    def __init__(self,master,*args,**kwargs):
        tk.Toplevel.__init__(self,master,*args,*kwargs)
        self.height = 500
        self.width = 600

        self.resizable(False, False)



        canvas = tk.Canvas(self, height=self.height, width=self.width)
        canvas.pack()
        canvas.configure(background=master.BG_COLOR)

        self.table_frame = tk.Frame(self,bg = '#76AA03')
        self.table_frame.place(relx=0.5,rely=0.2,relwidth=0.7,relheight=0.7,anchor='n')

        self.table_label = tk.Label(self.table_frame, bg=master.LIME)
        self.table_label.place(relwidth=1.0, relheight=1.0)
