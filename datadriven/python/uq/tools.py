import os
import re
import gzip

from pysgpp import DataMatrix, DataVector, Grid
import numpy as np
from math import log, trunc, exp


def addValidSequenceNumber(names=None, n=5):
    k = 0
    while names:
        s = str(k).zfill(n)
        dumps = dict([(key, name % s) for key, name in names.items()])

        if all([not os.path.exists(dump) for dump in dumps.values()]):
            return dumps

        k += 1


def gzOpen(filename, mode="r"):
    # gzip-file
    if re.match(".*\.gz$", filename):
        # mode set for binary data?
        if not mode[-1] == "b":
            mode += "b"
        fd = gzip.open(filename, mode)
    # non gzip-file
    else:
        fd = open(filename, mode)
    return fd


def createGrid(dim, level, borderType, isFull=False):
    from pysgpp.extensions.datadriven.learner.Types import BorderTypes
    if borderType == BorderTypes.NONE:
        grid = Grid.createLinearGrid(dim)
    elif borderType == BorderTypes.TRAPEZOIDBOUNDARY:
        grid = Grid.createLinearTrapezoidBoundaryGrid(dim)
    elif borderType == BorderTypes.COMPLETEBOUNDARY:
        grid = Grid.createLinearBoundaryGrid(dim, 0)
    else:
        raise Exception('Unknown border type')

    # create regular grid of level accLevel
    gridGen = grid.getGenerator()
    if isFull:
        gridGen.full(level)
    else:
        gridGen.regular(level)

    return grid


def natural_sort(accLevel):
    def convert(text): return \
        int(text) if text.isdigit()else text.lower()

    def alphanum_key(key): return \
        [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(accLevel, key=alphanum_key)


def unique(accLevel):
    seen = set()
    seen_add = seen.add
    return [x for x in accLevel if x not in seen and not seen_add(x)]


# # @brief write ARFF data
#
# @param data data to write
# @param merge flag to merge data (default False)
def writeDataARFF(data, merge=False):
    if len(data) == 0:
        return
    # is data a dataset or is it a list of datasets?
    if isinstance(data, dict):
        data = [data]

    hasHeader = False
    for dataset in data:
        if hasHeader is False or merge is False:
            hasHeader = True
            fout = gzOpen(dataset["filename"], "w")
            fout.write("@RELATION \"%s\"\n\n" % dataset["filename"])
            fstring = ""

            if isinstance(dataset["data"], DataMatrix):
                dim = dataset["data"].getNcols()
            else:
                dim = len(dataset["data"])

            for i in xrange(dim):
                if 'names' not in dataset:
                    fout.write("@ATTRIBUTE x%d NUMERIC\n" % i)
                else:
                    fout.write("@ATTRIBUTE %s NUMERIC\n" % dataset['names'][i])
                fstring = fstring + "%s,"

            hasclass = False
            if 'classes' in dataset:
                hasclass = True
                fout.write("@ATTRIBUTE class NUMERIC\n")
                fstring = fstring + "%s"
            else:
                fstring = fstring.strip(',')

            fstring = fstring + "\n"
            fout.write("\n@DATA\n")

        if isinstance(dataset["data"], DataMatrix):
            num_rows = dataset["data"].getNrows()
            for row in xrange(num_rows):
                lout = []
                for column in xrange(dim):
                    lout.append(dataset["data"].get(row, column))
                if hasclass:
                    lout.append(dataset["classes"][row])
                fout.write(fstring % tuple(lout))
        else:
            num_rows = len(dataset["data"][0])
            for row in xrange(num_rows):
                lout = []
                for column in xrange(dim):
                    lout.append(dataset["data"][column][row])
                if hasclass:
                    lout.append(dataset["classes"][row])
                fout.write(fstring % tuple(lout))

        if merge == False:
            fout.close()

    if merge == True:
        fout.close()
    return


def check(n, dim, nmax=50000):
    if dim > log(nmax) / log(n):
        n = trunc(exp(log(nmax) / dim))

    return n, n ** dim


def eval_linear(ppd, dim, norm_offset=0.0):
    n, nrows = check(ppd, dim)
    out = np.zeros(nrows * dim).reshape(nrows, dim)

    if nrows > 1:
        # Mesh width
        mw = (1 - 2 * norm_offset) / (1. * n - 1)
        # Generate a linear grid over the function range
        for j in range(dim):
            offset = n ** (dim - j - 1)
            i = 0
            while i < nrows:
                for k in xrange(n):
                    for accLevel in xrange(offset):
                        v = norm_offset + mw * k
                        out[i + k * offset + accLevel, j] = v
                i += n * offset

    return out


def eval_fullGrid(level, dim, border=True):
    if border:
        grid = Grid.createLinearBoundaryGrid(dim, 1)
    else:
        grid = Grid.createLinearGrid(dim)

    grid.getGenerator().full(level)
    gs = grid.getStorage()
    ans = np.ndarray((gs.getSize(), dim))
    p = DataVector(dim)

    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        ans[i, :] = p.array()

    return ans


def writeCSV(filename, samples, delim=' '):
    fd = open(filename, 'w')
    p = DataVector(samples.getNcols())
    for i in xrange(samples.getNrows()):
        samples.getRow(i, p)
        fd.write(delim.join(str(p)[1:-1].split(', ')) + '\n')
    fd.close()
