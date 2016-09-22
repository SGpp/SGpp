from strategies.QuadratureFactory import QuadratureFactory


class HashQuadrature(object):
    """
    Generic object for quadrature strategies
    """
    def __init__(self, U=None, T=None, opQuad=None):
        """
        Constructor
        """
        self._map = HashQuadratureMap()
        self._U = U
        self._T = T

        if opQuad is None:
            opQuad = QuadratureFactory.findQuadratureStrategyByMeasure(U)
        self.__opQuad = opQuad

    def quad(self, *args, **kws):
        return self.__opQuad.quad(*args, **kws)


# -----------------------------------------------------------------
class HashQuadratureMap(object):

    def __init__(self):
        """
        Constructor
        """
        self._map = {}

    def getKey(self, gps, d):
        """
        Generates a unique key for a given list of grid points
        @param gps: list of HashGridPoint
        @param d: int dimension
        """
        return tuple([(gp.getLevel(d), gp.getIndex(d)) for gp in gps])

    def __getitem__(self, key):
        return self._map[key]

    def __setitem__(self, key, value):
        self._map[key] = value

    def __contains__(self, key):
        return key in self._map
