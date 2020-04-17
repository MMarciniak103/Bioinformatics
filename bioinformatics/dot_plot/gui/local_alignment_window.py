import tkinter as tk
from alignment_algorithms.local_alignment_impl import LocalAlignment
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap
from matplotlib.lines import Line2D

class LocalAlignmentWindow(tk.Toplevel):
    def __init__(self, master, *args, **kwargs):
        tk.Toplevel.__init__(self, master, *args, *kwargs)
        self.height = 700
        self.width = 700

        self.alignments = []
        self.print_text = ""
        self.score_matrix = None
        self.path_matrix = None

        self.mappings = {0:'A',1:'C',2:'T',3:'G',4:'U',5:'-'}
        # Default substitution matrix values
        self.substitution_matrix = {
            'A': {'A': 1, 'C': -1, 'T': -1, 'G': -1, 'U': -1, '-': -3},

            'C': {'A': -1, 'C':1, 'T': -1, 'G':  -1, 'U': -1, '-':-3},

            'T': {'A': -1, 'C': -1, 'T': 1, 'G': -1, 'U': -1, '-': -3},

            'G': {'A': -1, 'C': -1, 'T': -1, 'G': 1, 'U': -1, '-': -3},

            'U': {'A': -1, 'C': -1, 'T': -1, 'G': -1, 'U': 1, '-': -3},

            '-': {'A':-3, 'C': -3, 'T': -3, 'G': -3, 'U': -3, '-':-3}
        }



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

        # ------------------------------------Table containing alignment information-----------------------------------
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

        try:
            for i in range(len(self.smatrix_entries)):
                for j in range(len(self.smatrix_entries[i])):
                    if self.smatrix_entries[i][j].get() == '':
                        tk.messagebox.showerror("ERROR", "You must provide sub cost for every possible state!")
                        return
                    else:
                        self.substitution_matrix[self.mappings[i]][self.mappings[j]] = float(self.smatrix_entries[i][j].get())
                        # substitution_matrix[mappings[i]][mappings[j]]= self.smatrix_entries[i][j].get()
        except Exception:
            tk.messagebox.showerror("ERROR", "There was an error while parsing matrix values!")
            return

        self.alignments = []

        enlargement_cost = 0

        affine_flag = self.affine_flag.get()

        # Check if user wants to use affine cost for gaps or linear one.
        if affine_flag == 1:
            try:
                enlargement_cost = float(self.affine_entry.get())
            except Exception:
                tk.messagebox.showerror("ERROR","Ivalid value for enlargement cost!")
                return

        LA = LocalAlignment(master.sequences,self.substitution_matrix)

        #Find all optimal alignments and score associated with them
        self.alignments,score,self.score_matrix,self.path_matrix = LA.predict_alignment(affine_flag,enlargement_cost)
        align_text_placeholder = "alignment" if len(self.alignments) == 1 else "alignments"
        self.print_text = f"Found {len(self.alignments)} optimal local {align_text_placeholder}.\n" \
                          f"Score value is: {score}"
        self.table_label['text'] = self.print_text
        plt.close('all')
        fig = plt.figure()
        sns.heatmap(self.score_matrix, annot=self.path_matrix,fmt='',annot_kws={"color": 'aqua'})
        path_mark = Line2D([0], [0], color='w', marker='x', linewidth=0)
        plt.legend([path_mark], ['alignment path'])
        plt.title('Local Alignment')
        plt.xlabel(self.wrap_text(master.sequences[0].get_sequence_name()))
        plt.ylabel(self.wrap_text(master.sequences[1].get_sequence_name()))
        plt.tight_layout()
        plt.show()



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
        Creates entry for every cell in substitution matrix. It also put their references in list with a shape (6,6).
        It sets default values to be equal with substitution matrix contained in program memory.
        :param frame: fram upon which entries would be placed
        :param font: font type
        :return:
        """
        for i in range(6):
            row = []
            for j in range(6):
                entry = tk.Entry(frame, font=font)
                entry.place(relx=0.3 + j * 0.1, rely=0.2 + i * 0.1, relwidth=0.1, relheight=0.1, anchor='n')
                value = self.substitution_matrix[self.mappings[i]][self.mappings[j]]
                entry.insert(value,str(value))
                row.append(entry)
            self.smatrix_entries.append(row)

    def _show_affine_widgets(self):
        if (self.affine_flag.get() == 1):
            self.affine_label.place(relx=0.3, rely=0.85, relwidth=0.2, relheight=0.1, anchor='n')
            self.affine_entry.place(relx=0.5, rely=0.85, relwidth=0.2, relheight=0.1, anchor='n')
        else:
            self.affine_label.place_forget()
            self.affine_entry.place_forget()


    def wrap_text(self,text):
        return '\n'.join(wrap(text,50))