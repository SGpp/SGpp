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

def load_font_properties():
    fontP = FontProperties()
    fontP.set_size(load_font()['size'])
    fontP.set_family(load_font()["family"])
    return fontP

#     colors = [None] * n
#     np.random.seed(123)
#
#     for i in range(n):
#         colors[i] = '#%06X' % np.random.randint(0, 0xFFFFFF)
#
#     return colors

def savefig(fig, filename, lgd=None):
    fig.tight_layout()
    if lgd is None:
        fig.savefig("%s.png" % filename, bbox_inches='tight')
        fig.savefig("%s.pdf" % filename, bbox_inches='tight')
        tikz_save("%s.tex" % filename)
    else:
        fig.savefig("%s.png" % filename,
                    bbox_extra_artists=(lgd,),
                    bbox_inches='tight')
        fig.savefig("%s.pdf" % filename,
                    bbox_extra_artists=(lgd,),
                    bbox_inches='tight')
        try:
            tikz_save("%s.tex" % filename)
        except:
            pass

    plt.close(fig)
