import tkinter as tk
from tkinter import ttk

class AlignmentWindow(tk.Toplevel):
    def __init__(self,master,*args,**kwargs):
        tk.Toplevel.__init__(self,master,*args,*kwargs)
        self.height = 500
        self.width = 400

        self.resizable(False, False)


        canvas = tk.Canvas(self, height=self.height, width=self.width)
        canvas.pack()

