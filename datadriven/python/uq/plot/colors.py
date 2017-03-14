import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
try:
    from matplotlib2tikz import save as tikz_save
except:
    pass

def intToRGB(i):
    blue = i & 255
    green = (i >> 8) & 255
    red = (i >> 16) & 255
    return red, green, blue

def rgbTpInt(rgb):
    red, green, blue = rgb
    i = (red << 16) + (green << 8) + blue
    return i


def load_color(i):
    colors = list(plt.rcParams['axes.prop_cycle'])
    return colors[i % len(colors)]["color"]

def load_marker(i):
    markers = ["o", "v", "D", "s", "*", "d", "^"]
    return markers[i % len(markers)]


def load_font():
    return {'family':'serif',
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

def savefig(fig, filename, lgd=None, tikz=True, mpl3d=False):
    fig.tight_layout()
    if mpl3d:
        fig.savefig("%s.png" % filename)
        fig.savefig("%s.pdf" % filename)
    else:
        if lgd is None:
            fig.savefig("%s.png" % filename, bbox_inches='tight')
            fig.savefig("%s.pdf" % filename, bbox_inches='tight')
            if tikz:
                try:
                    tikz_save("%s.tex" % filename)
                except:
                    pass
        else:
            fig.savefig("%s.png" % filename,
                        bbox_extra_artists=(lgd,),
                        bbox_inches='tight')
            fig.savefig("%s.pdf" % filename,
                        bbox_extra_artists=(lgd,),
                        bbox_inches='tight')
            if tikz:
                try:
                    tikz_save("%s.tex" % filename)
                except:
                    pass

    plt.close(fig)


def insert_legend(fig, loc="right", ncol=3, has_axis=True):
    if loc == "right":
        lgd = plt.legend(loc='upper left',
                         bbox_to_anchor=(1.02, 1),
                         borderaxespad=0,
                         prop=load_font_properties())
    elif loc == "bottom":
        lgd = plt.legend(loc='upper center',
                         ncol=ncol,
                         bbox_to_anchor=(0.5, -0.3) if has_axis else (0.5, -0.08),
                         borderaxespad=0,
                         prop=load_font_properties())
    elif loc == "top":
        lgd = ax.legend(loc='upper center',
                        bbox_to_anchor=(0.5, 1.25),
                        ncol=ncol,
                        borderaxespad=0,
                        prop=load_font_properties())
    elif loc == "left":
        lgd = plt.legend(loc='upper right',
                         bbox_to_anchor=(-0.1, 1),
                         borderaxespad=0,
                         prop=load_font_properties())
    else:
        raise AttributeError("loc '%s' not known" % loc)

    try:
        plt.setp(lgd.get_title(),
                 multialignment='left')
        for txt in lgd.get_texts():
            txt.set_ha('left')  # ha is alias for horizontalalignment
    except:
        pass

    plt.draw()  # to know size of legend
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.67, box.height])

    return lgd
