import tkinter as tk
from tkinter.constants import *
from tkinter import ttk, filedialog
from alignment_algorithms.star_alignment_impl import StarAlignment

class StarAlignmentWindow(tk.Toplevel):
    def __init__(self,master,*args,**kwargs):
        tk.Toplevel.__init__(self,master,*args,**kwargs)
        self.height = 700
        self.width = 700

        self.resizable(False,False)

        self.clustal_windows = []

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
        self.table_frame.place(relx=0.5, rely=0.3, relwidth=0.7, relheight=0.25, anchor='n')
        self.table_label = tk.Label(self.table_frame, bg=master.LIME, font=master.font_type)
        self.table_label.place(relwidth=1.0, relheight=1.0)


        self.alignment_btn_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.alignment_btn_frame.place(relx=0.5, rely=0.6, relwidth=0.9, relheight=0.35, anchor='n')

        self.star_alignment_btn = tk.Button(self.alignment_btn_frame, text='ALIGNMENT', font=master.font_type, bg=master.LIME, command=lambda:self.star_alignment(master))
        self.star_alignment_btn.place(relx=0.5, rely=0.15, relwidth=0.3, relheight=0.2, anchor='n')

        self.show_clustal_btn = tk.Button(self.alignment_btn_frame,text = 'SHOW CLUSTAL', font=master.font_type,bg=master.LIME,command = lambda :self.show_clustal(master))
        self.show_clustal_btn.place(relx=0.5,rely=0.4,relwidth=0.3,relheight=0.2,anchor='n')


        self.save_btn = tk.Button(self.alignment_btn_frame, text="SAVE", font=master.font_type, bg=master.LIME,command=lambda :self.save_alignment(master))
        self.save_btn.place(relx=0.5, rely=0.65, relwidth=0.3, relheight=0.2, anchor='n')


    def star_alignment(self,master):
        """
        Looks for optimal multiple sequences alignment of given sequences. It uses STAR algorithm
        :param master: master window
        """
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

    def _make_clustal_markers(self,i):
        """
        Provides marker's path that includes '*' symbol if given column is fully conserved.
        :param i: starting position of fragment that is used to construct marker's path
        :return: created marker's path
        """
        markers = ''

        for k in range(i,(i+50)):
            if k >= len(self.alignments[0]):
                break
            bad = False
            symbol = self.alignments[0][k]
            for j in range(len(self.alignments)):
                if self.alignments[j][k] != symbol:
                    bad = True
                    break

            if bad:
                markers = markers + '|'
            else:
                markers = markers + '*'

        return markers

    def _prepare_clustal_content(self, master):
        """
        Prepares text content in a clustal format.
        :param master: mater window
        :return: transformed text containing alignments
        """
        sequences_len = len(self.alignments[0])
        seq_names = [seq.get_sequence_name().split('.')[0] for seq in master.sequences]
        # seq_names = []
        beggining_space = len(seq_names[0])+2

        content = ""

        for i in range(0,sequences_len,50):
            for j in range(len(seq_names)):
                content += f"{seq_names[j]}: {self.alignments[j][i:i+50]}\n"

            for _ in range(beggining_space):
                content += "|"
            content += self._make_clustal_markers(i)
            content += "\n"

        content = content.replace("|"," ")

        return content

    def _prepare_save_content(self,master):
        seq_names = [seq.get_sequence_name() for seq in master.sequences]
        content = ">Multiple Sequences Alignment "
        for i,name in enumerate(seq_names):
            content += f",seq{i+1}: {name} "

        content += "\n"

        for align in self.alignments:
            content += f"{align}\n"

        return content


    def show_clustal(self,master):
        if len(self.alignments) == 0:
            tk.messagebox.showerror("ERROR", "There is no alignment!")
            return

        for window in self.clustal_windows:
            window.destroy()

        save_text = self._prepare_clustal_content(master)

        win = tk.Toplevel()
        win.wm_title("MSA CLUSTAL")

        scroll = ttk.Scrollbar(win)
        scroll.pack(side=RIGHT,fill=Y)

        listbox = tk.Listbox(win,yscrollcommand=scroll.set)

        for line in save_text.split('\n'):
            listbox.insert(END,line)
        listbox.pack(fill=BOTH,expand=1,padx=100)
        listbox.configure(justify=RIGHT,font='Courier')
        scroll.config(command=listbox.yview)



        self.clustal_windows.append(win)

    def save_alignment(self,master):
        """
        Saves multiple sequences alignment in a fasta format
        :param master: master window
        """
        if len(self.alignments) == 0:
            tk.messagebox.showerror("ERROR", "There is no alignment!")
            return

        f = filedialog.asksaveasfile(mode='w', defaultextension=".fasta")
        if f is None:  # If user didn't choose any file
            return

        save_text = self._prepare_save_content(master)

        for line in save_text:
            f.write(line)
        f.close()

