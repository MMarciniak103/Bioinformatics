import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import ttk


class StarAlignmentWindow(tk.Toplevel):
    def __init__(self,master,*args,**kwargs):
        tk.Toplevel.__init__(self,master,*args,**kwargs)
        self.height = 700
        self.width = 700

        self.resizable(False,False)

        self.alignments = []

        canvas = tk.Canvas(self,height=self.height,width=self.width)
        canvas.pack()
        canvas.configure(background = master.BG_COLOR)

        self.cost_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.cost_frame.place(relx=0.5, rely=0.05, relwidth=0.8, relheight=0.1, anchor='n')
        # INDENT
        self.indent_cost_label = tk.Label(self.cost_frame, bg=master.LIME, font=master.font_type)
        self.indent_cost_label.place(relx=0.25, rely=0.2, relwidth=0.2, relheight=0.3, anchor='n')
        self.indent_cost_label['text'] = "INDENT COST"
        self.indent_cost_entry = tk.Entry(self.cost_frame, font=master.font_type)
        self.indent_cost_entry.place(relx=0.25, rely=0.5, relwidth=0.2, relheight=0.4, anchor='n')
       # SUBSTITUTION
        self.substitution_cost_label = tk.Label(self.cost_frame, bg=master.LIME, font=master.font_type)
        self.substitution_cost_label.place(relx=0.5, rely=0.2, relwidth=0.25, relheight=0.3, anchor='n')
        self.substitution_cost_label['text'] = "SUBSTITUTION COST"
        self.substitution_cost_entry = tk.Entry(self.cost_frame, font=master.font_type)
        self.substitution_cost_entry.place(relx=0.5, rely=0.5, relwidth=0.25, relheight=0.4, anchor='n')
        # MATCH
        self.match_cost_label = tk.Label(self.cost_frame, bg=master.LIME, font=master.font_type)
        self.match_cost_label.place(relx=0.75, rely=0.2, relwidth=0.2, relheight=0.3, anchor='n')
        self.match_cost_label['text'] = "MATCH COST"
        self.match_cost_entry = tk.Entry(self.cost_frame, font=master.font_type)
        self.match_cost_entry.place(relx=0.75, rely=0.5, relwidth=0.2, relheight=0.4, anchor='n')

