from matplotlib.font_manager import FontProperties
import os
import subprocess
import pwd

import matplotlib.pyplot as plt
import numpy as np
try:
    from matplotlib2tikz import save as tikz_save
except:
    pass


def get_username():
    """
    returns the user name
    """
    return pwd.getpwuid(os.getuid())[0]


def load_custom_pgf_preamble(dtype="standard", macros="thesis"):
    """
    Load preamble for latex style plots. It includes amsmath, scientific units,
    font style for springer publications, etc. If the macros parameter
    is defined, then all the provided latex commands in there will be loaded
    in the preamble as well.

    @param dtype: string that specifies the math style (standard, springer)
    @param macros: path to latex file that contains latex commands
    """
    pysgpp_uq_font = load_font()

    pgf_preamble = {"font.family": pysgpp_uq_font["family"],  # use serif/main font for text elements
                    "text.usetex": True,  # use inline math for ticks
                    "text.latex.preamble": [r"\usepackage[utf8x]{inputenc}",
                                            r'\usepackage{amsmath}',
                                            r"\usepackage{amssymb}",
                                            r"\usepackage{tikz}",
                                            r"\usepackage{pgfplots}",
                                            r'\usepackage[scientific-notation=true]{siunitx}'
                                            ],
                    'axes.labelsize': pysgpp_uq_font["size"],
                    'font.size': pysgpp_uq_font["size"],
                    'legend.fontsize': pysgpp_uq_font["size"],
                    'xtick.labelsize': pysgpp_uq_font["size"],
                    'ytick.labelsize': pysgpp_uq_font["size"],
                    'axes.unicode_minus': True,
                    #                     'figure.figsize': (5, 4.5),
                    'image.cmap': load_default_color_map(dtype="string")
                    #                     'axes.titlepad': 25
                    }

    if dtype == "springer":
        pgf_preamble["text.latex.preamble"] += [r"\usepackage{mathptmx}",
                                                r"\usepackage{rsfso}"]
    else:
        pgf_preamble["text.latex.preamble"] += [r"\usepackage[T1]{fontenc}"]

    if get_username() == "franzefn":
        if macros == "thesis":
            cmd_filename = r"/home/franzefn/Promotion/UQ/repos/dissertation/thesis/commands.tex"
        elif macros == "l2leja":
            cmd_filename = r"/home/franzefn/Promotion/Paper/Awesome-CT-Leja-Papers-of-Dabian-Holzelin/l2-leja/paper/commands.tex"
        elif macros == "l2leja_david":
            cmd_filename = r".../Awesome-CT-Leja-Papers-of-Dabian-Holzelin/l2-leja/paper/commands.tex"
        else:
            cmd_filename = r"/home/franzefn/Promotion/Paper/repos/SGA16/paper/commands.tex"
    elif get_username() == "holzmudd":
        if macros == "l2leja":
            cmd_filename = r".../Awesome-CT-Leja-Papers-of-Dabian-Holzelin/l2-leja/paper/commands.tex"
        else:
            cmd_filename = "plots/commands.tex"
    elif get_username() == "rehmemk":
        cmd_filename = "plots/commands.tex"
    else:
        cmd_filename = r"/home/franzefn/Promotion/Paper/repos/SGA16/paper/commands.tex"

    if os.path.exists(cmd_filename):
        fd = open(cmd_filename, "r")
        for paramstring in fd.readlines():
            decoded_paramstring = paramstring.decode('utf8')
            if decoded_paramstring[0] not in ["%", "\n"]:
                pgf_preamble["text.latex.preamble"].append(decoded_paramstring)

    return pgf_preamble


def initialize_plotting_style(dtype="standard",
                              macros="thesis",
                              style="seaborn-paper"):
    """
    Initialize the plotting style for the current python session.

    @param dtype: font style type (standard, springer)
    @param macros: path to commands latex file -> containing (newcommand)
    @param style: plot style
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    plt.style.use(style)

    # Include packages `amssymb` and `amsmath` in LaTeX preamble
    # as they include extended math support (symbols, envisonments etc.)
    pgf_with_custom_preamble = load_custom_pgf_preamble(dtype=dtype,
                                                        macros=macros)
    mpl.rcParams.update(pgf_with_custom_preamble)


def intToRGB(i):
    blue = i & 255
    green = (i >> 8) & 255
    red = (i >> 16) & 255
    return red, green, blue


def rgbTpInt(rgb):
    red, green, blue = rgb
    i = (red << 16) + (green << 8) + blue
    return i


def load_default_color_map(dtype="cmap"):
    """
    Load default color map -> viridis
    @param dtype: color map identifier
    """
    if dtype == "cmap":
        return plt.get_cmap('viridis')
    else:
        return 'viridis'


def load_color(i):
    """
    Provide color via index
    @param i: index
    """
    colors = list(plt.rcParams['axes.prop_cycle'])
    return colors[i % len(colors)]["color"]


def load_marker(i):
    """
    Provide marker via index
    @param i: index
    """
    markers = ["o", "v", "D", "s", "*", "d", "^", "x", "+"]
    return markers[i % len(markers)]


def load_linestyle(i):
    """
    Provide linestyle via index.
    @param i: index
    """
    linestyles = ["-", ":", "-.", "--"]
    return linestyles[i % len(linestyles)]


def load_bw_color(i, nmax=11):
    """
    Load black and white color via index

    @param i: index
    @param nmax: maximum number of cycle
    """
    colors = np.linspace(0, 0.75, nmax, endpoint=True, dtype=str)
    return colors[i % len(colors)]


def load_font():
    """
    Load standard font
    """
    return {'family': 'serif',
            'size': 20}


def load_font_properties(size=None,
                         family=None):
    fontP = FontProperties()
    fontP.set_size(load_font()['size'] if size is None else size)
    fontP.set_family(load_font()["family"] if family is None else family)
    return fontP

#     colors = [None] * n
#     np.random.seed(123)
#
#     for i in range(n):
#         colors[i] = '#%06X' % np.random.randint(0, 0xFFFFFF)
#
#     return colors


def savefig(fig,
            filename,
            lgd=None,
            tikz=False,
            mpl3d=False,
            crop=False,
            sizes=(5.25, 5)):
    """
    Save a figure in various formats.

    @param fig: figure object
    @param filename: filename (without extension) where the file should be stored
    @param lgd: legend object
    @param tikz: save result as tikz file using matplotlib2tikz
    @param mpl3d: is it a 3d plot? Then one should not use a tight layout
    @param crop: crop the resulting pdf -> pdfcrop needs to be installed
    @param sizes: tuple containing the figures size in x and y direction in inches
    """
    if sizes is not None:
        sizex, sizey = sizes
        fig.set_size_inches(sizex, sizey, forward=True)
    if mpl3d:
        fig.savefig("%s.png" % filename)
        fig.savefig("%s.pdf" % filename)
    else:
        # fig.tight_layout()
        if lgd is None:
            fig.savefig("%s.png" % filename, bbox_inches='tight')
            fig.savefig("%s.pdf" % filename, bbox_inches='tight')
        else:
            fig.savefig("%s.png" % filename,
                        bbox_extra_artists=(lgd,),
                        bbox_inches='tight')
            fig.savefig("%s.pdf" % filename,
                        bbox_extra_artists=(lgd,),
                        bbox_inches='tight')
        if tikz:
            try:
                tikz_save("%s.tex" % filename, fig)
            except:
                pass

    if crop:
        subprocess.call(["pdfcrop",
                         "%s.pdf" % filename,
                         "%s.pdf" % filename])

    plt.close(fig)


def insert_legend(fig, loc="right", ncol=3, has_axis=True, shift=0.0):
    """
    Inserts a legend outside of the plotting area. It is possible that the
    exact position (bbox_to_anchor) needs to be adjusted but once the
    size of the figures is fixed, this works fine.

    @param fig: figure object
    @param loc: location of the legend (right, bottom, top, left)
    @param ncol: number of columns for the labels
    @param has_axis: defines if xticks labels are shown
    @param shift: additional offset of the legend from the inner plotting area
    """
    if loc == "right":
        lgd = plt.legend(loc='upper left',
                         bbox_to_anchor=(1.02 + shift, 1),
                         borderaxespad=0,
                         prop=load_font_properties())
    elif loc == "bottom":
        lgd = plt.legend(loc='upper center',
                         ncol=ncol,
                         bbox_to_anchor=(0.5, -0.3 + shift) if has_axis else (0.5, -0.08 + shift),
                         borderaxespad=0,
                         prop=load_font_properties())
    elif loc == "top":
        lgd = plt.legend(loc='upper center',
                         bbox_to_anchor=(0.5, 1.25 + shift),
                         ncol=ncol,
                         borderaxespad=0,
                         prop=load_font_properties())
    elif loc == "left":
        lgd = plt.legend(loc='upper right',
                         bbox_to_anchor=(-0.33 + shift, 1),
                         borderaxespad=0,
                         prop=load_font_properties())
    else:
        lgd = plt.legend(loc=loc, prop=load_font_properties())

    try:
        plt.setp(lgd.get_title(),
                 multialignment='left')
        for txt in lgd.get_texts():
            txt.set_ha('left')  # ha is alias for "horizontal alignment"
    except:
        pass

    plt.draw()  # to know size of legend
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.67, box.height])

    return lgd
