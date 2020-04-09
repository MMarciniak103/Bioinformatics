import tkinter as tk
import numpy as np
from alignment_algorithms.local_alignment_impl import LocalAlignment
from collections import defaultdict

class LocalAlignmentWindow(tk.Toplevel):
    def __init__(self, master, *args, **kwargs):
        tk.Toplevel.__init__(self, master, *args, *kwargs)
        self.height = 700
        self.width = 700

        self.resizable(False, False)

        # This list is used as an container for substitution matrix entries
        self.smatrix_entries = []

        canvas = tk.Canvas(self, height=self.height, width=self.width)
        canvas.pack()
        canvas.configure(background=master.BG_COLOR)

        # -----------------------------------SUBSTITUTION MATRIX COMPONENTS ---------------------------------------------
        self.smatrix_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.smatrix_frame.place(relx=0.5, rely=0.1, relwidth=0.8, relheight=0.4, anchor='n')

        self._create_smatrix_labels(frame=self.smatrix_frame, bg_color=master.LIME, font=master.font_type)
        self._create_smatrix_entries(frame=self.smatrix_frame, font=master.font_type)

        self.affine_flag = tk.IntVar()
        self.affine_checkbox = tk.Checkbutton(self.smatrix_frame, text="affine penalty",
                                              variable=self.affine_flag,
                                              onvalue=1,
                                              offvalue=0,
                                              command=lambda: self._show_affine_widgets())
        self.affine_checkbox.place(relx=0.8, rely=0.85, relwidth=0.2, relheight=0.1, anchor='n')

        self.affine_label = tk.Label(self.smatrix_frame, text='Enlargement cost', bg=master.LIME, font=master.font_type)
        self.affine_entry = tk.Entry(self.smatrix_frame, font=master.font_type)

        # ------------------------------------Table containing alignment informations-----------------------------------
        self.table_frame = tk.Frame(self, bg=master.BG_COLOR)
        self.table_frame.place(relx=0.5, rely=0.5, relwidth=0.7, relheight=0.35, anchor='n')
        self.table_label = tk.Label(self.table_frame, bg=master.LIME, font=master.font_type)
        self.table_label.place(relwidth=1.0, relheight=1.0)

        # ------------------------  Widgets for alignment plot and save button selection -------------------------------
        self.alignment_btn_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.alignment_btn_frame.place(relx=0.5, rely=0.8, relwidth=0.9, relheight=0.15, anchor='n')

        self.global_alignment_btn = tk.Button(self.alignment_btn_frame, text='ALIGNMENT', font=master.font_type,
                                              bg=master.LIME,command=lambda:self.predict_local_alignemnt(master))
        self.global_alignment_btn.place(relx=0.5, rely=0.15, relwidth=0.3, relheight=0.4, anchor='n')

        self.save_btn = tk.Button(self.alignment_btn_frame, text="SAVE", font=master.font_type, bg=master.LIME)
        self.save_btn.place(relx=0.5, rely=0.65, relwidth=0.3, relheight=0.4, anchor='n')


    def predict_local_alignemnt(self,master):
        mappings = {0:'A',1:'C',2:'T',3:'G',4:'U',5:'-'}
        substitution_matrix = {
            'A': {'A': 1, 'C': 1, 'T': 1, 'G': 1, 'U': 1, '-': 1},

            'C': {'A': 1, 'C':1, 'T': 1, 'G':  1, 'U': 1, '-':1},

            'T': {'A': 1, 'C': 1, 'T': 1, 'G': 1, 'U': 1, '-': 1},

            'G': {'A': 1, 'C': 1, 'T': 1, 'G': 1, 'U': 1, '-': 1},

            'U': {'A': 1, 'C': 1, 'T': 1, 'G': 1, 'U': 1, '-': 1},

            '-': {'A':1, 'C': 1, 'T': 1, 'G': 1, 'U': 1, '-':1}
        }

        try:
            for i in range(len(self.smatrix_entries)):
                for j in range(len(self.smatrix_entries[i])):
                    if self.smatrix_entries[i][j].get() == '':
                        tk.messagebox.showerror("ERROR", "You must provide sub cost for every possible state!")
                        return
                    else:
                        substitution_matrix[mappings[i]][mappings[j]] = float(self.smatrix_entries[i][j].get())
                        # substitution_matrix[mappings[i]][mappings[j]]= self.smatrix_entries[i][j].get()
        except Exception:
            tk.messagebox.showerror("ERROR", "There was an error while parsing matrix values!")
            return

        enlargement_cost = 0

        affine_flag = self.affine_flag.get()

        if affine_flag == 1:
            try:
                enlargement_cost = float(self.affine_entry.get())
            except Exception:
                tk.messagebox.showerror("ERROR","Ivalid value for enlargement cost!")
                return

        LA = LocalAlignment(master.sequences,substitution_matrix)

        alignments = LA.predict_alignment(affine_flag,enlargement_cost)

        print(alignments)


    def _create_smatrix_labels(self, frame, bg_color, font):
        """
        Method that creates labels for substitution matrix
        :param frame: frame upon which labels would be placed
        :param bg_color: background color
        :param font: font type
        :return:
        """
        symbols = ['A', 'C', 'T', 'G', 'U', '-']
        for i in range(6):
            # Create label for given symbol row-wise
            label = tk.Label(frame, bg=bg_color, font=font)
            label.place(relx=0.2, rely=0.2 + i * 0.1, relwidth=0.1, relheight=0.1, anchor='n')
            label['text'] = symbols[i]
            # Create label for given symbol column wise
            label = tk.Label(frame, bg=bg_color, font=font)
            label.place(relx=0.3 + i * 0.1, rely=0.1, relwidth=0.1, relheight=0.1, anchor='n')
            label['text'] = symbols[i]

    def _create_smatrix_entries(self, frame, font):
        """
        Creates entry for every cell in substitution matrix. It also put their references in list with a shape (6,6)
        :param frame: fram upon which entries would be placed
        :param font: font type
        :return:
        """
        for i in range(6):
            row = []
            for j in range(6):
                entry = tk.Entry(frame, font=font)
                entry.place(relx=0.3 + j * 0.1, rely=0.2 + i * 0.1, relwidth=0.1, relheight=0.1, anchor='n')
                row.append(entry)
            self.smatrix_entries.append(row)

    def _show_affine_widgets(self):
        if (self.affine_flag.get() == 1):
            self.affine_label.place(relx=0.3, rely=0.85, relwidth=0.2, relheight=0.1, anchor='n')
            self.affine_entry.place(relx=0.5, rely=0.85, relwidth=0.2, relheight=0.1, anchor='n')
        else:
            self.affine_label.place_forget()
            self.affine_entry.place_forget()
