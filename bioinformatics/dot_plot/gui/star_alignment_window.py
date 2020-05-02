import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import ttk
from alignment_algorithms.star_alignment_impl import StarAlignment

class StarAlignmentWindow(tk.Toplevel):
    def __init__(self,master,*args,**kwargs):
        tk.Toplevel.__init__(self,master,*args,**kwargs)
        self.height = 700
        self.width = 700

        self.resizable(False,False)

        self.alignments = []
        self.msa_score = None
        self.opt_score = None
        self.print_text = ""

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

        # ------------------------------------Table containing alignment information-----------------------------------
        self.table_frame = tk.Frame(self, bg=master.BG_COLOR)
        self.table_frame.place(relx=0.5, rely=0.3, relwidth=0.7, relheight=0.35, anchor='n')
        self.table_label = tk.Label(self.table_frame, bg=master.LIME, font=master.font_type)
        self.table_label.place(relwidth=1.0, relheight=1.0)


        self.alignment_btn_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.alignment_btn_frame.place(relx=0.5, rely=0.8, relwidth=0.9, relheight=0.15, anchor='n')

        self.star_alignment_btn = tk.Button(self.alignment_btn_frame, text='ALIGNMENT', font=master.font_type, bg=master.LIME, command=lambda:self.star_alignment(master))
        self.star_alignment_btn.place(relx=0.5, rely=0.15, relwidth=0.3, relheight=0.4, anchor='n')



    def star_alignment(self,master):

        indent_cost = self.indent_cost_entry.get()
        substitution_cost = self.substitution_cost_entry.get()
        match_cost = self.match_cost_entry.get()
        if (indent_cost == '' or substitution_cost == '' or match_cost == ''):
            tk.messagebox.showerror("ERROR", "Provide cost value for all possible states!")
            return

        self.alignments = []
        self.msa_score = None
        self.opt_score = None
        self.print_text = ""

        sa = StarAlignment(master.sequences)
        self.alignments,self.msa_score,self.opt_score = sa.predict_alignment(float(indent_cost),float(substitution_cost),float(match_cost))
        self.print_text = f"STAR algorithm found optimal alignment\n" \
                          f"Estimated cost is in range:\n" \
                          f"[{self.opt_score} - {self.msa_score}]"

        self.table_label['text'] = self.print_text