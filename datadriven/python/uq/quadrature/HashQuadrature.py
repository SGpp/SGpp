from strategies.QuadratureFactory import QuadratureFactory


class HashQuadrature(object):
    """
    Generic object for quadrature strategies
    """
    def __init__(self, gridType=None, U=None, T=None, opQuad=None):
        """
        Constructor
        """
        self._map = HashQuadratureMap()
        self._U = U
        self._T = T
        self._gridType = gridType

        if opQuad is None:
            opQuad = QuadratureFactory.findQuadratureStrategyByMeasure(U)
        self.__opQuad = opQuad

    def setDistributionAndTransformation(self, U, T):
        self._U = U
        self._T = T

    def setGridType(self, gridType):
        self._gridType = gridType

    def quad(self, *args, **kws):
        return self.__opQuad.quad(*args, **kws)


# -----------------------------------------------------------------
class HashQuadratureMap(object):

    def __init__(self):
        """
        Constructor
        """
        self._map = {}

    def getKey(self, dist, gps, d=None):
        """
        Generates a unique key for a given list of grid points
        @param gps: list of HashGridPoint
        @param d: int dimension
        """
        if d is None:
            return tuple([str(dist)] + [(gp.getLevel(d), gp.getIndex(d))
                                        for gp in gps
                                        for d in xrange(gp.getDimension())])
        else:
            return tuple([str(dist)] + [(gp.getLevel(d), gp.getIndex(d))
                                        for gp in gps])

    def __getitem__(self, key):
        return self._map[key]

    def __setitem__(self, key, value):
        self._map[key] = value

    def __contains__(self, key):
        return key in self._map
