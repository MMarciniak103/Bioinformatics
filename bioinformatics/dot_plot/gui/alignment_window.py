import tkinter as tk
from tkinter import ttk

class AlignmentWindow(tk.Toplevel):
    def __init__(self,master,*args,**kwargs):
        tk.Toplevel.__init__(self,master,*args,*kwargs)
        self.height = 600
        self.width = 700

        self.resizable(False, False)



        canvas = tk.Canvas(self, height=self.height, width=self.width)
        canvas.pack()
        canvas.configure(background=master.BG_COLOR)

        #---------------------------------Widgets for cost values selection---------------------------------------------
        self.cost_frame = tk.Frame(self,bg=master.BG_COLOR,bd=5)
        self.cost_frame.place(relx=0.5,rely=0.05,relwidth=0.8,relheight=0.1,anchor='n')
        #INSERTION
        self.insertion_cost_label = tk.Label(self.cost_frame,bg=master.LIME,font=master.font_type)
        self.insertion_cost_label.place(relx=0.2,rely=0.2,relwidth=0.25,relheight=0.3,anchor='n')
        self.insertion_cost_label['text'] = "INSERTION COST"
        self.insertion_cost_entry = tk.Entry(self.cost_frame,font=master.font_type)
        self.insertion_cost_entry.place(relx=0.2,rely=0.5,relwidth=0.25,relheight=0.4,anchor='n')
        #DELETION
        self.deletion_cost_label = tk.Label(self.cost_frame, bg=master.LIME, font=master.font_type)
        self.deletion_cost_label.place(relx=0.5, rely=0.2, relwidth=0.25, relheight=0.3, anchor='n')
        self.deletion_cost_label['text'] = "DELETION COST"
        self.deletion_cost_entry = tk.Entry(self.cost_frame, font=master.font_type)
        self.deletion_cost_entry.place(relx=0.5, rely=0.5, relwidth=0.25, relheight=0.4, anchor='n')
        #SUBSTITUTION
        self.substitution_cost_label = tk.Label(self.cost_frame, bg=master.LIME, font=master.font_type)
        self.substitution_cost_label.place(relx=0.8, rely=0.2, relwidth=0.25, relheight=0.3, anchor='n')
        self.substitution_cost_label['text'] = "SUBSTITUTION COST"
        self.substitution_cost_entry = tk.Entry(self.cost_frame, font=master.font_type)
        self.substitution_cost_entry.place(relx=0.8, rely=0.5, relwidth=0.25, relheight=0.4, anchor='n')

        #------------------------------------Table containing alignment informations------------------------------------
        self.table_frame = tk.Frame(self,bg = master.BG_COLOR)
        self.table_frame.place(relx=0.5,rely=0.2,relwidth=0.7,relheight=0.6,anchor='n')
        self.table_label = tk.Label(self.table_frame, bg=master.LIME)
        self.table_label.place(relwidth=1.0, relheight=1.0)

        #------------------------------------Widgets for alignment type selection --------------------------------------
        self.alignment_btn_frame = tk.Frame(self,bg=master.BG_COLOR,bd=5)
        self.alignment_btn_frame.place(relx=0.5,rely=0.85,relwidth=0.9,relheight=0.1,anchor='n')

        self.global_alignment_btn = tk.Button(self.alignment_btn_frame,text='GLOBAL',font=master.font_type,bg=master.LIME)
        self.global_alignment_btn.place(relx=0.3,rely=0.4,relwidth=0.3,relheight=0.5,anchor='n')

        self.local_alignment_btn = tk.Button(self.alignment_btn_frame,text="LOCAL",font=master.font_type,bg=master.LIME)
        self.local_alignment_btn.place(relx=0.7,rely=0.4,relwidth=0.3,relheight=0.5,anchor='n')



    def global_alignment(self):
        pass
