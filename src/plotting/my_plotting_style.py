#!/usr/bin/env python3
''' This should improve the graphics in the document by using TeX for text,
by using the same font as the report, and using reasonable text sizes.
Import into any script that plots Figures for the report.
Incurs a large time penalty, which is worth it for final plots.
Blank entries should cause plots to inherit fonts from the document.'''

import matplotlib

pgf_latex = {"pgf.texsystem": "pdflatex",
             "text.usetex": True,
             "font.family": "serif",
             "font.serif": [],
             "font.sans-serif": [],
             "font.monospace": [],
             "lines.linewidth": 0.75,
             "axes.labelsize": 10,
             "font.size": 10,
             "legend.fontsize": 8,
             "xtick.labelsize": 10,
             "ytick.labelsize": 10,
             "pgf.preamble": [r"\usepackage[utf8x]{inputenc}",
                              r"\usepackage[T1]{fontenc}",
                              r"\usepackage{siunitx}",
                              ]
             }

matplotlib.rcParams['text.latex.preamble'] =\
    ['\\usepackage[separate -uncertainty=true]{siunitx}']

matplotlib.rcParams.update(pgf_latex)

#---text sizes are meaningless if LaTeX squashes figure---


'''
usage:
    fig = plt.gcf()
    fig = single_column_width(fig)
    plt.tight_layout()
or:
    fig = plt.gcf()
    fig = defined_size(fig, size='double_column_square')
    plt.tight_layout()
or 
'''

def single_column_width(fig):
    fig = defined_size(fig, size='single_column')
    return fig


def defined_size(fig, size='single_column', **kwargs):
    size_dict = {'single_column':[3.125, 2.5],
                 'double_column_short':[6.5, 2.5],
                 'double_column_square':[6.5, 6.5],
                 'gal_an_size':[7, 6.5]}
    fig.set_size_inches(size_dict[size])
    return fig

