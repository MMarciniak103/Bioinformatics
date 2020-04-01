import tkinter as tk
from alignment_algorithms.alignment_impl import GlobalAlignment
from tkinter import filedialog
from textwrap import wrap
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D


class GlobalAlignmentWindow(tk.Toplevel):
    def __init__(self, master, *args, **kwargs):
        tk.Toplevel.__init__(self, master, *args, *kwargs)
        self.height = 700
        self.width = 700

        self.resizable(False, False)

        self.print_text = ""
        self.alignments = []
        self.alignment_matrix = None
        self.path_matrix = None

        canvas = tk.Canvas(self, height=self.height, width=self.width)
        canvas.pack()
        canvas.configure(background=master.BG_COLOR)

        # ---------------------------------Widgets for cost values selection---------------------------------------------
        self.cost_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.cost_frame.place(relx=0.5, rely=0.05, relwidth=0.8, relheight=0.1, anchor='n')
        # INSERTION
        self.insertion_cost_label = tk.Label(self.cost_frame, bg=master.LIME, font=master.font_type)
        self.insertion_cost_label.place(relx=0.1, rely=0.2, relwidth=0.2, relheight=0.3, anchor='n')
        self.insertion_cost_label['text'] = "INSERTION COST"
        self.insertion_cost_entry = tk.Entry(self.cost_frame, font=master.font_type)
        self.insertion_cost_entry.place(relx=0.1, rely=0.5, relwidth=0.2, relheight=0.4, anchor='n')
        # DELETION
        self.deletion_cost_label = tk.Label(self.cost_frame, bg=master.LIME, font=master.font_type)
        self.deletion_cost_label.place(relx=0.35, rely=0.2, relwidth=0.2, relheight=0.3, anchor='n')
        self.deletion_cost_label['text'] = "DELETION COST"
        self.deletion_cost_entry = tk.Entry(self.cost_frame, font=master.font_type)
        self.deletion_cost_entry.place(relx=0.35, rely=0.5, relwidth=0.2, relheight=0.4, anchor='n')
        # SUBSTITUTION
        self.substitution_cost_label = tk.Label(self.cost_frame, bg=master.LIME, font=master.font_type)
        self.substitution_cost_label.place(relx=0.6, rely=0.2, relwidth=0.25, relheight=0.3, anchor='n')
        self.substitution_cost_label['text'] = "SUBSTITUTION COST"
        self.substitution_cost_entry = tk.Entry(self.cost_frame, font=master.font_type)
        self.substitution_cost_entry.place(relx=0.6, rely=0.5, relwidth=0.25, relheight=0.4, anchor='n')
        # MATCH
        self.match_cost_label = tk.Label(self.cost_frame, bg=master.LIME, font=master.font_type)
        self.match_cost_label.place(relx=0.85, rely=0.2, relwidth=0.2, relheight=0.3, anchor='n')
        self.match_cost_label['text'] = "MATCH COST"
        self.match_cost_entry = tk.Entry(self.cost_frame, font=master.font_type)
        self.match_cost_entry.place(relx=0.85, rely=0.5, relwidth=0.2, relheight=0.4, anchor='n')

        # ------------------------------------Table containing alignment informations------------------------------------
        self.table_frame = tk.Frame(self, bg=master.BG_COLOR)
        self.table_frame.place(relx=0.5, rely=0.2, relwidth=0.7, relheight=0.6, anchor='n')
        self.table_label = tk.Label(self.table_frame, bg=master.LIME, font=master.font_type)
        self.table_label.place(relwidth=1.0, relheight=1.0)

        # ------------------------  Widgets for alignment plot and save button selection -------------------------------
        self.alignment_btn_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.alignment_btn_frame.place(relx=0.5, rely=0.8, relwidth=0.9, relheight=0.15, anchor='n')

        self.global_alignment_btn = tk.Button(self.alignment_btn_frame, text='ALIGNMENT', font=master.font_type,
                                              bg=master.LIME, command=lambda: self.global_alignment(master))
        self.global_alignment_btn.place(relx=0.5, rely=0.15, relwidth=0.3, relheight=0.4, anchor='n')


        self.save_btn = tk.Button(self.alignment_btn_frame, text="SAVE", font=master.font_type, bg=master.LIME,
                                  command=lambda: self.save_alignment())
        self.save_btn.place(relx=0.5, rely=0.65, relwidth=0.3, relheight=0.4, anchor='n')

    def wrap_text(self,text):
        return '\n'.join(wrap(text,50))

    def format_table_data(self, alignments, score, gaps, identity, cost_array):
        '''
        Formats text that would be printed in statistics window. This formatted text is also part of report that is able
        to be saved in txt file
        :param alignments: found alignments
        :param score: score of alignment (for distance metric we want to minimize it)
        :param gaps: number of gaps in alignment
        :param identity: number of matching positions in alignment
        :param cost_array: array containing cost values for different scenarios.
        :return:
        '''
        self.print_text = "STATISTICS:\n " \
                          f"seq1: {self.master.sequences[0].get_sequence_name()}\n" \
                          f"seq2: {self.master.sequences[1].get_sequence_name()}\n" \
                          f"Mode: distance\n " \
                          f"Insertion cost : {cost_array[0]}\n" \
                          f"Deletion cost : {cost_array[1]}\n" \
                          f"Substitution cost : {cost_array[2]}\n" \
                          f"Match cost : {cost_array[3]}\n" \
                          f"Score {score}\n" \
                          f"Length: {len(alignments[0])}\n" \
                          f"Gaps: {gaps}/{len(alignments[0])} ({(gaps/len(alignments[0])*100):.2f}%)\n" \
                          f"Identity: {identity}/{len(alignments[0])} ({(identity/len(alignments[0])*100):.2f}%)\n"


    def prepare_txt_to_show(self):
        print_text = self.print_text.split('\n')
        print_text[1] = f"seq1: {self.wrap_text(self.master.sequences[0].get_sequence_name())}\n"
        print_text[2] = f"seq2: {self.wrap_text(self.master.sequences[1].get_sequence_name())}\n"
        show = '\n'.join([e for e in print_text])
        return show

    def global_alignment(self, master):
        """
        Find global alignment of 2 given sequences with specified cost values for different scenarios.
        After the alignment is found, it would print statistics of calculations and show score heatmap with
        optimal alignment path.
        :param master: master window
        :return:
        """
        self.print_text = ""
        self.alignments = []
        self.alignment_matrix = None
        self.path_matrix = None

        insertion_cost = self.insertion_cost_entry.get()
        deletion_cost = self.deletion_cost_entry.get()
        substitution_cost = self.substitution_cost_entry.get()
        match_cost = self.match_cost_entry.get()
        if (insertion_cost == '' or deletion_cost == '' or substitution_cost == '' or match_cost == ''):
            tk.messagebox.showerror("ERROR", "Provide cost value for all possible states!")
            return

        ga = GlobalAlignment(master.sequences)
        self.alignments, score, gaps, identity, self.alignment_matrix, self.path_matrix = ga.predict_alignment(
            insertion_cost=float(insertion_cost), deletion_cost=float(deletion_cost),
            substitution_cost=float(substitution_cost), match_cost=float(match_cost))
        if len(self.alignments[0]) != 0:
            cost_array = [insertion_cost, deletion_cost, substitution_cost, match_cost]
            self.format_table_data(self.alignments, score, gaps, identity, cost_array)
            self.table_label['text'] = self.prepare_txt_to_show()

        plt.close('all')
        fig = plt.figure()
        sns.heatmap(self.alignment_matrix, annot=self.path_matrix, fmt='')
        path_mark = Line2D([0], [0], color='w', marker='x', linewidth=0)
        plt.legend([path_mark], ['alignment path'])
        plt.title('Global Alignment')
        plt.xlabel(self.wrap_text(self.master.sequences[0].get_sequence_name()))
        plt.ylabel(self.wrap_text(self.master.sequences[1].get_sequence_name()))
        plt.tight_layout()
        plt.show()

    def save_alignment(self):
        '''
        Method that is used to save calculated alignment with its statistics.
        :return:
        '''
        if len(self.alignments) != 3:
            tk.messagebox.showerror("ERROR", "There is no alignment!")
            return

        f = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
        if f is None:  # If user didn't choose any file
            return
        save_text = self.print_text + "\n" \
                                      f"{self.alignments[0]}\n" \
                                      f"{self.alignments[2]}\n" \
                                      f"{self.alignments[1]}"
        for line in save_text:
            f.write(line)
        f.close()
