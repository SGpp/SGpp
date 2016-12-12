import numpy as np

def intToRGB(i):
    blue = i & 255
    green = (i >> 8) & 255
    red = (i >> 16) & 255
    return red, green, blue

def rgbTpInt(rgb):
    red, green, blue = rgb
    i = (red << 16) + (green << 8) + blue
    return i


def loadColorSequence(n):
    colors = [None] * n
    np.random.seed(1234567)

    for i in range(n):
        colors[i] = '#%06X' % np.random.randint(0, 0xFFFFFF)

    return colors
