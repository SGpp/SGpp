from BilinearQuadratureStrategy import BilinearQuadratureStrategy
from scipy.integrate import quad, dblquad


class DiscreteBilinearScipyQuadratureStrategy(BilinearQuadratureStrategy):
    """
    Use Scipy for quadrature
    """

    def __init__(self):
        """
        Constructor
        """
        super(self.__class__, self).__init__()

    def computeBilinearFormEntry(self, gpi, basisi, gpj, basisj):
        # if not, compute it
        val = 1
        err = 0.

        # run over all dimensions
        d = 0
        while d < gpi.dim() and abs(val) > 1e-15:
            # get level index
            lid, iid = gpi.getLevel(d), gpi.getIndex(d)
            ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

            # compute left and right boundary of the support of both
            # basis functions
            xlow = max([(iid - 1) * 2 ** -lid, (ijd - 1) * 2 ** -ljd])
            xhigh = min([(iid + 1) * 2 ** -lid, (ijd + 1) * 2 ** -ljd])

            # same level, different index
            if lid == ljd and iid != ijd and lid > 0:
                val = err = 0
            # the support does not overlap
            elif lid != ljd and xlow >= xhigh:
                val = err = 0
            else:
                # ----------------------------------------------------
                # use scipy for integration
                def f(p):
                    return basisi.eval(lid, iid, p) * \
                        basisj.eval(ljd, ijd, p)

                s, err1d = quad(f, xlow, xhigh, epsabs=1e-8)
                # ----------------------------------------------------
                val *= s
                err += err1d

            d += 1

        return val, err
