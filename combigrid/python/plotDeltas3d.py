# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# inspired by https://stackoverflow.com/questions/18859935/painting-a-cube

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib import colors, cm
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.pyplot import figure, show

    def quad(plane='xy', origin=None, width=1, height=1, depth=0):
        """returns a rectangle"""
        u, v = (0, 0) if origin is None else origin

        plane = plane.lower()
        if plane == 'xy':
            vertices = ((u, v, depth),
                        (u + width, v, depth),
                        (u + width, v + height, depth),
                        (u, v + height, depth))
        elif plane == 'xz':
            vertices = ((u, depth, v),
                        (u + width, depth, v),
                        (u + width, depth, v + height),
                        (u, depth, v + height))
        elif plane == 'yz':
            vertices = ((depth, u, v),
                        (depth, u + width, v),
                        (depth, u + width, v + height),
                        (depth, u, v + height))
        else:
            raise ValueError('"{0}" is not a supported plane!'.format(plane))
        return np.array(vertices)

    def subcube(origin=(0, 0, 0), width=1, height=1, depth=1,):
        u, v, w = origin
        quads = [
            quad('xy', (u, v), width, height, w),
            quad('xy', (u, v), width, height, w+depth),
            quad('yz', (v, w), height, depth, u),
            quad('yz', (v, w), height, depth, u+width),
            quad('xz', (u, w), width, depth, v),
            quad('xz', (u, w), width, depth, v+height)
        ]
        return np.array(quads)

    def subcube_at(midpoint=(1, 1, 1)):
        halfwidth = 0.4
        width = 2*halfwidth
        return subcube((midpoint[0]-halfwidth, midpoint[1]-halfwidth, midpoint[2]-halfwidth), width, width, width)

    def get_projected_cubes_and_colors(dataframe, dim, scalarMap):
        dframe = dataframe.copy()
        # select those with the maximum delta for that projection
        dframe = dframe.groupby(dim, group_keys=False).apply(
            lambda x: x.loc[x.delta.idxmax()])
        dframe.drop_duplicates(subset=dim, inplace=True, keep='last')
        subcubes = []
        colors = []
        for i, l in dframe.iterrows():
            # generate a cube for each
            #         print((l[dim[0]], l[dim[1]], l[dim[2]]))
            # print((l[dim[0]], l[dim[1]]))
            if len(dim) == 2:
                p = subcube_at((l[dim[0]], l[dim[1]], 1))
            else:
                p = subcube_at((l[dim[0]], l[dim[1]], l[dim[2]]))
            subcubes.append(p)
            color = scalarMap.to_rgba(l['delta'], alpha=0.3)
            for i in range(6):
                # add color for each side of the cube
                colors.append(color)
        return subcubes, colors

    def adaptiveGeneratorToDataframe(adaptiveGenerator):
        df = pd.DataFrame()
        levels = adaptiveGenerator.getLevels()
        deltas = adaptiveGenerator.getDeltas(levels)
        levels = list(levels)
        for i in range(len(levels)):
            lInfo = {"l_" + str(j): (levels[i][j])
                     for j in range(len(levels[i]))}
            lInfo["delta"] = float(deltas[i])
            df = df.append(lInfo, ignore_index=True)
        # drop potential nan values
        df = df.dropna()
        return df

    def plotDeltas3D(dataframe, dim=['l_0', 'l_1', 'l_2']):
        """print 3d cube representation of adaptive scheme"""
        """expects a pandas dataframe with column 'delta' and level indices 'l_*'"""
        canvas = figure()
        axes = Axes3D(canvas)

        if not isinstance(dataframe, pd.core.frame.DataFrame):
            dataframe = adaptiveGeneratorToDataframe(dataframe)

        max_delta = dataframe['delta'].max()
        cmap = plt.get_cmap('magma_r')
        norm = colors.Normalize(vmin=0., vmax=max_delta)
        scalarMap = cm.ScalarMappable(norm=norm, cmap=cmap)
        scalarMap.set_array(dataframe['delta'].values)

        q = subcube(origin=(0.5, 0.5, 0.5),
                    width=0.5,
                    height=0.5,
                    depth=0.5,)
        p = subcube_at((0., 0., 0.))

        a, c = get_projected_cubes_and_colors(dataframe, dim, scalarMap)

        allcubes = np.concatenate(a, axis=0)
        collection = Poly3DCollection(allcubes)
        collection.set_color(c)
        axes.set_xlim(dataframe[dim[0]].min() - 0.5,
                      dataframe[dim[0]].max() + 0.5)
        axes.set_xlabel(dim[0])
        axes.set_ylim(dataframe[dim[1]].max() + 0.5,
                      dataframe[dim[1]].min() - 0.5)
        axes.set_ylabel(dim[1])
        if len(dim) > 2:
            axes.set_zlim(dataframe[dim[2]].min() - 0.5,
                          dataframe[dim[2]].max() + 0.5)
            axes.set_zlabel(dim[2])
        else:
            range_zero = dataframe[dim[0]].max() - dataframe[dim[0]].min()
            range_one = dataframe[dim[1]].max() - dataframe[dim[1]].min()
            axes.set_zlim(0.6, 1+max(range_zero, range_one))
        axes.add_collection3d(collection)

        cbar = canvas.colorbar(scalarMap)

        show()


except ImportError as e:
    print(e)
    print("Cannot use pysgpp 3D delta visualization")
    pass
