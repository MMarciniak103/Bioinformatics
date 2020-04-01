import tkinter as tk
import numpy as np

class LocalAlignmentWindow(tk.Toplevel):
    def __init__(self, master, *args, **kwargs):
        tk.Toplevel.__init__(self, master, *args, *kwargs)
        self.height = 700
        self.width = 700

        self.resizable(False, False)

        #This list is used as an container for substitution matrix entries
        self.smatrix_entries = []

        canvas = tk.Canvas(self, height=self.height, width=self.width)
        canvas.pack()
        canvas.configure(background=master.BG_COLOR)

        #-----------------------------------SUBSTITUTION MATRIX COMPONENTS ---------------------------------------------
        self.smatrix_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.smatrix_frame.place(relx=0.5, rely=0.1, relwidth=0.8, relheight=0.4, anchor='n')

        self._create_smatrix_labels(frame=self.smatrix_frame,bg_color=master.LIME,font=master.font_type)
        self._create_smatrix_entries(frame=self.smatrix_frame,font=master.font_type)

        # ------------------------------------Table containing alignment informations-----------------------------------
        self.table_frame = tk.Frame(self, bg=master.BG_COLOR)
        self.table_frame.place(relx=0.5, rely=0.5, relwidth=0.7, relheight=0.35, anchor='n')
        self.table_label = tk.Label(self.table_frame, bg=master.LIME, font=master.font_type)
        self.table_label.place(relwidth=1.0, relheight=1.0)

        # ------------------------  Widgets for alignment plot and save button selection -------------------------------
        self.alignment_btn_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.alignment_btn_frame.place(relx=0.5, rely=0.8, relwidth=0.9, relheight=0.15, anchor='n')

        self.global_alignment_btn = tk.Button(self.alignment_btn_frame, text='ALIGNMENT', font=master.font_type,
                                              bg=master.LIME)
        self.global_alignment_btn.place(relx=0.5, rely=0.15, relwidth=0.3, relheight=0.4, anchor='n')

        self.save_btn = tk.Button(self.alignment_btn_frame, text="SAVE", font=master.font_type, bg=master.LIME)
        self.save_btn.place(relx=0.5, rely=0.65, relwidth=0.3, relheight=0.4, anchor='n')



    def _create_smatrix_labels(self,frame,bg_color,font):
        """
        Method that creates labels for substitution matrix
        :param frame: frame upon which labels would be placed
        :param bg_color: background color
        :param font: font type
        :return:
        """
        symbols = ['A','C','T','G','U','-']
        for i in range(6):
            #Create label for given symbol row-wise
            label = tk.Label(frame,bg=bg_color,font=font)
            label.place(relx=0.2,rely=0.2+i*0.1,relwidth=0.1,relheight=0.1,anchor='n')
            label['text'] = symbols[i]
            #Create label for given symbol column wise
            label = tk.Label(frame,bg=bg_color,font=font)
            label.place(relx=0.3+i*0.1,rely=0.1,relwidth=0.1,relheight=0.1,anchor='n')
            label['text'] = symbols[i]


    def _create_smatrix_entries(self,frame,font):
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
                entry.place(relx=0.3+j*0.1, rely=0.2+i*0.1, relwidth=0.1, relheight=0.1, anchor='n')
                row.append(entry)
            self.smatrix_entries.append(row)
