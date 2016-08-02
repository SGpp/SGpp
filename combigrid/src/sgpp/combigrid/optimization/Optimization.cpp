// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include "Optimization.h"

namespace sgpp {
namespace combigrid {

namespace optimize {

double glomin(double a, double b, double c, double m, double e, double t, func_base& f, double& x) {
  //    GLOMIN seeks a global minimum of a function F(X) in an interval [A,B].
  //
  //  Discussion:
  //
  //    This function assumes that F(X) is twice continuously differentiable
  //    over [A,B] and that F''(X) <= M for all X in [A,B].
  //
  //  Parameters:
  //
  //    Input, double A, B, the endpoints of the interval.
  //    It must be the case that A < B.
  //
  //    Input, double C, an initial guess for the global
  //    minimizer.  If no good guess is known, C = A or B is acceptable.
  //
  //    Input, double M, the bound on the second derivative.
  //
  //    Input, double E, a positive tolerance, a bound for the
  //    absolute error in the evaluation of F(X) for any X in [A,B].
  //
  //    Input, double T, a positive error tolerance.
  //
  //    Input, func_base& F, a user-supplied c++ functor whose
  //    global minimum is being sought.  The input and output
  //    of F() are of type double.
  //
  //    Output, double &X, the estimated value of the abscissa
  //    for which F attains its global minimum value in [A,B].
  //
  //    Output, double GLOMIN, the value F(X).
  //
  double a0;
  double a2;
  double a3;
  double d0;
  double d1;
  double d2;
  double h;
  int k;
  double m2;
  double macheps;
  double p;
  double q;
  double qs;
  double r;
  double s;
  double sc;
  double y;
  double y0;
  double y1;
  double y2;
  double y3;
  double yb;
  double z0;
  double z1;
  double z2;

  a0 = b;
  x = a0;
  a2 = a;
  y0 = f(b);
  yb = y0;
  y2 = f(a);
  y = y2;

  if (y0 < y) {
    y = y0;
  } else {
    x = a;
  }

  if (m <= 0.0 || b <= a) {
    return y;
  }

  macheps = r8_epsilon();

  m2 = 0.5 * (1.0 + 16.0 * macheps) * m;

  if (c <= a || b <= c) {
    sc = 0.5 * (a + b);
  } else {
    sc = c;
  }

  y1 = f(sc);
  k = 3;
  d0 = a2 - sc;
  h = 9.0 / 11.0;

  if (y1 < y) {
    x = sc;
    y = y1;
  }
  for (;;) {
    d1 = a2 - a0;
    d2 = sc - a0;
    z2 = b - a2;
    z0 = y2 - y1;
    z1 = y2 - y0;
    r = d1 * d1 * z0 - d0 * d0 * z1;
    p = r;
    qs = 2.0 * (d0 * z1 - d1 * z0);
    q = qs;

    if (k < 1000000 || y2 <= y) {
      for (;;) {
        if (q * (r * (yb - y2) + z2 * q * ((y2 - y) + t)) < z2 * m2 * r * (z2 * q - r)) {
          a3 = a2 + r / q;
          y3 = f(a3);

          if (y3 < y) {
            x = a3;
            y = y3;
          }
        }
        k = ((1611 * k) % 1048576);
        q = 1.0;
        r = (b - a) * 0.00001 * static_cast<double>(k);

        if (z2 <= r) {
          break;
        }
      }
    } else {
      k = ((1611 * k) % 1048576);
      q = 1.0;
      r = (b - a) * 0.00001 * static_cast<double>(k);

      while (r < z2) {
        if (q * (r * (yb - y2) + z2 * q * ((y2 - y) + t)) < z2 * m2 * r * (z2 * q - r)) {
          a3 = a2 + r / q;
          y3 = f(a3);

          if (y3 < y) {
            x = a3;
            y = y3;
          }
        }
        k = ((1611 * k) % 1048576);
        q = 1.0;
        r = (b - a) * 0.00001 * static_cast<double>(k);
      }
    }

    r = m2 * d0 * d1 * d2;
    s = sqrt(((y2 - y) + t) / m2);
    h = 0.5 * (1.0 + h);
    p = h * (p + 2.0 * r * s);
    q = q + 0.5 * qs;
    r = -0.5 * (d0 + (z0 + 2.01 * e) / (d0 * m2));

    if (r < s || d0 < 0.0) {
      r = a2 + s;
    } else {
      r = a2 + r;
    }

    if (0.0 < p * q) {
      a3 = a2 + p / q;
    } else {
      a3 = r;
    }

    for (;;) {
      a3 = r8_max(a3, r);

      if (b <= a3) {
        a3 = b;
        y3 = yb;
      } else {
        y3 = f(a3);
      }

      if (y3 < y) {
        x = a3;
        y = y3;
      }

      d0 = a3 - a2;

      if (a3 <= r) {
        break;
      }

      p = 2.0 * (y2 - y3) / (m * d0);

      if ((1.0 + 9.0 * macheps) * d0 <= r8_abs(p)) {
        break;
      }

      if (0.5 * m2 * (d0 * d0 + p * p) <= (y2 - y) + (y3 - y) + 2.0 * t) {
        break;
      }
      a3 = 0.5 * (a2 + a3);
      h = 0.9 * h;
    }

    if (b <= a3) {
      break;
    }

    a0 = sc;
    sc = a2;
    a2 = a3;
    y0 = y1;
    y1 = y2;
    y2 = y3;
  }

  return y;
}

double local_min(double a, double b, double t, func_base& f, double& x) {
  //    LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
  //
  //  Discussion:
  //
  //    The method used is a combination of golden section search and
  //    successive parabolic interpolation.  Convergence is never much slower
  //    than that for a Fibonacci search.  If F has a continuous second
  //    derivative which is positive at the minimum (which is not at A or
  //    B), then convergence is superlinear, and usually of the order of
  //    about 1.324....
  //
  //    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
  //    F is never evaluated at two points closer than TOL.
  //
  //    If F is a unimodal function and the computed values of F are always
  //    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
  //    LOCAL_MIN approximates the abscissa of the global minimum of F on the
  //    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
  //
  //    If F is not unimodal, then LOCAL_MIN may approximate a local, but
  //    perhaps non-global, minimum to the same accuracy.
  //
  //  Parameters:
  //
  //    Input, double A, B, the endpoints of the interval.
  //
  //    Input, double T, a positive absolute error tolerance.
  //
  //    Input, func_base& F, a user-supplied c++ functor whose
  //    local minimum is being sought.  The input and output
  //    of F() are of type double.
  //
  //    Output, double &X, the estimated value of an abscissa
  //    for which F attains a local minimum value in [A,B].
  //
  //    Output, double LOCAL_MIN, the value F(X).
  //
  double c;
  double d;
  double e;
  double eps;
  double fu;
  double fv;
  double fw;
  double fx;
  double m;
  double p;
  double q;
  double r;
  double sa;
  double sb;
  double t2;
  double tol;
  double u;
  double v;
  double w;
  //
  //  C is the square of the inverse of the golden ratio.
  //
  c = 0.5 * (3.0 - sqrt(5.0));

  eps = sqrt(r8_epsilon());

  sa = a;
  sb = b;
  x = sa + c * (b - a);
  w = x;
  v = w;
  e = 0.0;
  fx = f(x);
  fw = fx;
  fv = fw;

  for (;;) {
    m = 0.5 * (sa + sb);
    tol = eps * r8_abs(x) + t;
    t2 = 2.0 * tol;
    //
    //  Check the stopping criterion.
    //
    if (r8_abs(x - m) <= t2 - 0.5 * (sb - sa)) {
      break;
    }
    //
    //  Fit a parabola.
    //
    r = 0.0;
    q = r;
    p = q;

    if (tol < r8_abs(e)) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (0.0 < q) {
        p = -p;
      }
      q = r8_abs(q);
      r = e;
      e = d;
    }

    if (r8_abs(p) < r8_abs(0.5 * q * r) && q * (sa - x) < p && p < q * (sb - x)) {
      //
      //  Take the parabolic interpolation step.
      //
      d = p / q;
      u = x + d;
      //
      //  F must not be evaluated too close to A or B.
      //
      if ((u - sa) < t2 || (sb - u) < t2) {
        if (x < m) {
          d = tol;
        } else {
          d = -tol;
        }
      }
    } else {
      //
      //  A golden-section step.
      //
      if (x < m) {
        e = sb - x;
      } else {
        e = sa - x;
      }
      d = c * e;
    }
    //
    //  F must not be evaluated too close to X.
    //
    if (tol <= r8_abs(d)) {
      u = x + d;
    } else if (0.0 < d) {
      u = x + tol;
    } else {
      u = x - tol;
    }

    fu = f(u);
    //
    //  Update A, B, V, W, and X.
    //
    if (fu <= fx) {
      if (u < x) {
        sb = x;
      } else {
        sa = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    } else {
      if (u < x) {
        sa = u;
      } else {
        sb = u;
      }

      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  return fx;
}

double r8_abs(double x) {
  //    R8_ABS returns the absolute value of an R8.
  //
  //  Parameters:
  //
  //    Input, double X, the quantity whose absolute value is desired.
  //
  //    Output, double R8_ABS, the absolute value of X.
  double value;

  if (0.0 <= x) {
    value = x;
  } else {
    value = -x;
  }
  return value;
}

double r8_epsilon() {
  //    R8_EPSILON returns the R8 roundoff unit.
  //
  //  Discussion:
  //
  //    The roundoff unit is a number R which is a power of 2 with the
  //    property that, to the precision of the computer's arithmetic,
  //      1 < 1 + R
  //    but
  //      1 = ( 1 + R / 2 )
  //
  //  Parameters:
  //
  //    Output, double R8_EPSILON, the R8 round-off unit.
  const double value = 2.220446049250313E-016;

  return value;
}

double r8_max(double x, double y) {
  //    R8_MAX returns the maximum of two R8's.
  //
  //  Parameters:
  //
  //    Input, double X, Y, the quantities to compare.
  //
  //    Output, double R8_MAX, the maximum of X and Y.
  double value;

  if (y < x) {
    value = x;
  } else {
    value = y;
  }
  return value;
}

double r8_sign(double x) {
  //    R8_SIGN returns the sign of an R8.
  //
  //  Parameters:
  //
  //    Input, double X, the number whose sign is desired.
  //
  //    Output, double R8_SIGN, the sign of X.
  //
  double value;

  if (x < 0.0) {
    value = -1.0;
  } else {
    value = 1.0;
  }
  return value;
}

// ======================================================================
// === Simple wrapper functions
// === for convenience and/or compatibility.
//
// === The three functions are the same as above,
// === except that they take a plain function F
// === instead of a c++ functor.  In all cases, the
// === input and output of F() are of type double.

typedef double doubleOfdouble(double);

class func_wrapper : public func_base {
  doubleOfdouble* func;

 public:
  explicit func_wrapper(doubleOfdouble* f) { func = f; }
  virtual ~func_wrapper() {}

  virtual double operator()(double x) { return func(x); }
};

// ****************************************************************************80

double glomin(double a, double b, double c, double m, double e, double t, double f(double x),
              double& x) {
  func_wrapper foo(f);
  return glomin(a, b, c, m, e, t, foo, x);
}

// ****************************************************************************80

double local_min(double a, double b, double t, double f(double x), double& x) {
  func_wrapper foo(f);
  return local_min(a, b, t, foo, x);
}

// ======================================================================
// Generally useful functor to evaluate a monic polynomial.
// For details, see class definition in Optimize.h

double monicPoly::operator()(double x) {
  double rslt(1);
  for (int ii = static_cast<int>(coeff.size()) - 1; ii >= 0; ii--) {
    rslt *= x;
    rslt += coeff[ii];
  }
  return rslt;
}

// Similarly, evaluate a general polynomial (not necessarily monic):
double Poly::operator()(double x) {
  double rslt(0);
  for (int ii = static_cast<int>(coeff.size()) - 1; ii >= 0; ii--) {
    rslt *= x;
    rslt += coeff[ii];
  }
  return rslt;
}

}  // namespace optimize
}  // namespace combigrid
} /* namespace sgpp*/
