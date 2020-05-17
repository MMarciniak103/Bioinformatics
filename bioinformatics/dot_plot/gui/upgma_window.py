import tkinter as tk
from tkinter.constants import *
from tkinter import ttk, filedialog
import matplotlib.pyplot as plt
from alignment_algorithms.upgma_impl import UPGMA,flatten
import numpy as np


class UpgmaWindow(tk.Toplevel):
    def __init__(self,master,*args,**kwargs):
        tk.Toplevel.__init__(self,master,*args,**kwargs)
        self.height = 700
        self.width = 700

        self.resizable(False,False)
        self.connections = None
        self.seq_hash = None

        canvas = tk.Canvas(self, height=self.height, width=self.width)
        canvas.pack()
        canvas.configure(background=master.BG_COLOR)

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

        self.upgma_frame = tk.Frame(self, bg=master.BG_COLOR, bd=5)
        self.upgma_frame.place(relx=0.5, rely=0.6, relwidth=0.9, relheight=0.35, anchor='n')

        self.upgma_btn = tk.Button(self.upgma_frame, text='UPGMA', font=master.font_type,
                                   bg=master.LIME, command=lambda: self.upgma(master))
        self.upgma_btn.place(relx=0.5, rely=0.15, relwidth=0.3, relheight=0.2, anchor='n')

        self.show_clustal_btn = tk.Button(self.upgma_frame, text='DENDROGRAM', font=master.font_type,
                                          bg=master.LIME,command=lambda : self.show_dendrogram(master))
        self.show_clustal_btn.place(relx=0.5, rely=0.4, relwidth=0.3, relheight=0.2, anchor='n')
        #
        # self.save_btn = tk.Button(self.upgma_frame, text="SAVE", font=master.font_type, bg=master.LIME,
        #                           command=lambda: self.save_alignment(master))
        # self.save_btn.place(relx=0.5, rely=0.65, relwidth=0.3, relheight=0.2, anchor='n')


    def upgma(self,master):
        """
        Creates UPGMA dendrogram for given sequences
        :param master: master window
        :return:
        """
        indent_cost = self.indent_cost_entry.get()
        substitution_cost = self.substitution_cost_entry.get()
        match_cost = self.match_cost_entry.get()
        if (indent_cost == '' or substitution_cost == '' or match_cost == ''):
            tk.messagebox.showerror("ERROR", "Provide cost value for all possible states!")
            return

        upgma = UPGMA(master.sequences)
        self.connections,self.seq_hash = upgma.create_upgma(float(indent_cost),float(substitution_cost),float(match_cost))

    def _mk_fork(self,x0, x1, y0, y1, new_level):
        points = [[x0, x0, x1, x1], [y0, new_level, new_level, y1]]
        connector = [(x0 + x1) / 2., new_level]
        return (points), connector


    def show_dendrogram(self,master):
        if self.connections is None or self.seq_hash is None:
            tk.messagebox.showerror("ERROR!","You must first create upgma dendrogram")
            return

        levels = [connection[0] for connection in self.connections]
        seen = []
        locations = []
        for connection in self.connections:
            if isinstance(connection[1][0], list):
                first_elem = connection[1][0][0]
            else:
                first_elem = connection[1][0]
            if isinstance(connection[1][1], list):
                second_elem = connection[1][1][0]
            else:
                second_elem = connection[1][1]
            # locations.append(tuple(seq_hash[elem] for elem in flatten(connection[1])[:2]))
            locations.append((self.seq_hash[first_elem], self.seq_hash[second_elem]))
            seen += flatten([elem for elem in flatten(connection[1])])

        label_map = {i: {} for i in range(len(master.sequences))}
        for i in range(len(master.sequences)):
            label_map[i] = {'label': master.sequences[i].get_sequence_name().split('.')[0], 'xpos': i, 'ypos': 0}

        print('locations ', locations)
        print('label_map ', label_map)

        plt.close('all')
        fig, ax = plt.subplots()

        for i, (new_level, (loc0, loc1)) in enumerate(zip(levels, locations)):
            print('step {0}:\t connecting ({1},{2}) at level {3}'.format(i, loc0, loc1, new_level))

            x0, y0 = label_map[loc0]['xpos'], label_map[loc0]['ypos']
            x1, y1 = label_map[loc1]['xpos'], label_map[loc1]['ypos']

            print('\t points are: {0}:({2},{3}) and {1}:({4},{5})'.format(loc0, loc1, x0, y0, x1, y1))

            p, c = self._mk_fork(x0, x1, y0, y1, new_level)

            ax.plot(*p)
            ax.scatter(*c)
            print('\t connector is at:{0}'.format(c))

            label_map[loc0]['xpos'] = c[0]
            label_map[loc0]['ypos'] = c[1]
            label_map[loc0]['label'] = '{0}/{1}'.format(label_map[loc0]['label'], label_map[loc1]['label'])
            print('\t updating label_map[{0}]:{1}'.format(loc0, label_map[loc0]))

            ax.text(*c, label_map[loc0]['label'])

        _xticks = np.arange(0, 6, 1)
        _xticklabels = [*[seq.get_sequence_name().split('.')[0] for seq in master.sequences]]

        ax.set_xticks(_xticks)
        ax.set_xticklabels(_xticklabels)

        ax.set_ylim(0, 1.05 * np.max(levels))

        plt.show()
