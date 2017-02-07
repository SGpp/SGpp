import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

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
    return list(plt.rcParams['axes.prop_cycle'])[i]["color"]

def load_marker(i):
    return ["o", "v", "D", "s", "*", "d", "^"][i]


def load_font():
    return {'family':'serif',
            'size': 12}

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
