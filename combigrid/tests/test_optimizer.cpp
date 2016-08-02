// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/combigrid/optimization/Optimization.hpp>

#include <cstdlib>
#include <cmath>
#include <ctime>

const double check_tolerance(1e-15);

void test_local_min_all();
void test_glomin_all();

template <class T>
void test_local_min_one(double a, double b, double t, T f, double expected_xval,
                        double expected_yval);
void test_glomin_one(double a, double b, double c, double m, double e, double t, double f(double x),
                     double expected_xval, double expected_yval);

double g_01(double x);
double g_02(double x);
double g_03(double x);
double g_04(double x);
double g_05(double x);
double h_01(double x);
double h_02(double x);
double h_03(double x);
double h_04(double x);
double h_05(double x);

// in %
double tol = 0.01;

BOOST_AUTO_TEST_CASE(testOptimizer) {
  test_local_min_all();
  test_glomin_all();
}

void test_local_min_all() {
  double a;
  double b;
  double t;

  t = r8_epsilon();

  a = 0.0;
  b = 3.141592653589793;

  test_local_min_one(a, b, t, g_01, 2.0, 1.0);

  a = 0.0;
  b = 1.0;

  test_local_min_one(a, b, t, g_02, 0.351734, 0.827184);

  const int degree(4);
  double coefflist[1 + degree] = {3, 1, 2, 0, 1};
  Poly quadric(coefflist, degree);

  test_local_min_one(-2.0, 2.0, t, quadric, -0.236733, 2.87849);

  a = 0.0001;
  b = 1.0;

  test_local_min_one(a, b, t, g_04, 0.0953446, 1.20492);

  a = 0.0002;
  b = 2.0;

  test_local_min_one(a, b, t, g_05, 0.703205, 0.628026);

  return;
}

void test_glomin_all() {
  double a;
  double b;
  double c;
  double e;
  double m;
  double t;

  e = sqrt(r8_epsilon());
  t = sqrt(r8_epsilon());

  a = 7.0;
  b = 9.0;
  c = (a + b) / 2.0;
  m = 0.0;

  test_glomin_one(a, b, c, m, e, t, h_01, 9.0, -7.0);

  a = 7.0;
  b = 9.0;
  c = (a + b) / 2.0;
  m = 100.0;

  test_glomin_one(a, b, c, m, e, t, h_01, 9.0, -7.0);

  a = -1.0;
  b = +2.0;
  c = (a + b) / 2.0;
  m = 2.0;

  test_glomin_one(a, b, c, m, e, t, h_02, 0.0, 0.0);

  a = -1.0;
  b = +2.0;
  c = (a + b) / 2.0;
  m = 2.1;

  test_glomin_one(a, b, c, m, e, t, h_02, 0.0, 0.0);

  a = -0.5;
  b = +2.0;
  c = (a + b) / 2.0;
  m = 14.0;

  test_glomin_one(a, b, c, m, e, t, h_03, 0.0000000030521847430027881,
                  0.0000000000000000093158317340058548);

  a = -0.5;
  b = +2.0;
  c = (a + b) / 2.0;
  m = 28.0;

  test_glomin_one(a, b, c, m, e, t, h_03, 0.0000000062226175770753174,
                  0.000000000000000038720969751472479);

  a = -10.0;
  b = +10.0;
  c = (a + b) / 2.0;
  m = 72.0;

  test_glomin_one(a, b, c, m, e, t, h_04, -0.679579, -0.824239);

  a = -10.0;
  b = +10.0;
  c = (a + b) / 2.0;
  m = 72.0;

  test_glomin_one(a, b, c, m, e, t, h_05, -1.19514, -0.0634905);

  return;
}

template <class T>
void test_local_min_one(double a, double b, double t, T f, double expected_xval,
                        double expected_yval) {
  //  Parameters:
  //
  //    Input, double A, B, the endpoints of the interval.
  //
  //    Input, double T, a positive absolute error tolerance.
  //
  //    Input, double F ( double x ), the name of a user-supplied
  //    function, whose local minimum is being sought.
  //
  //    Input, expected_xval, y_val: expected values
  double fx;
  double x;

  fx = local_min(a, b, t, f, x);

  BOOST_CHECK_CLOSE(x, expected_xval, tol);
  BOOST_CHECK_CLOSE(fx, expected_yval, tol);

  return;
}

void test_glomin_one(double a, double b, double c, double m, double e, double t, double f(double x),
                     double expected_xval, double expected_yval) {
  //  Parameters:
  //
  //    Input, double A, B, the endpoints of the interval.
  //
  //    Input, double C, an initial guess for the global
  //    minimizer.  If no good guess is known, C = A or B is acceptable.
  //
  //    Input, double M, the bound on the second derivative.
  //
  //    Input, double E, a positive tolerance, a bound for the
  //    absolute error in the evaluation of F(X) for any X in [A,B].
  //
  //    Input, double T, a positive absolute error tolerance.
  //
  //    Input, double F ( double x ), the name of a user-supplied
  //    function whose global minimum is being sought.
  //
  //    Input, expected_xval, y_val: expected values
  double fx;
  double x;

  fx = glomin(a, b, c, m, e, t, f, x);

  BOOST_CHECK_CLOSE(x, expected_xval, tol);
  BOOST_CHECK_CLOSE(fx, expected_yval, tol);

  return;
}

double g_01(double x) {
  // g_01(x) = (x-2)^2 + 1
  return (x - 2.0) * (x - 2.0) + 1.0;
}

double g_02(double x) {
  // g_02(0) = x^2 + exp ( - x )
  return x * x + exp(-x);
}

double g_03(double x) {
  // g_03(x) = x^4+2x^2+x+3
  return ((x * x + 2.0) * x + 1.0) * x + 3.0;
}

double g_04(double x) {
  // g_04(x) = exp(x)+1/(100X)
  return exp(x) + 0.01 / x;
}

double g_05(double x) {
  // g_05(x) = exp(x) - 2x + 1/(100x) - 1/(1000000x^2)
  return exp(x) - 2.0 * x + 0.01 / x - 0.000001 / x / x;
}

double h_01(double x) {
  // h_01(x) = 2 - x.
  return 2.0 - x;
}

double h_02(double x) {
  // h_02(x) = x^2.
  return x * x;
}

double h_03(double x) {
  // h_03(x) = x^3+x^2.
  return x * x * (x + 1.0);
}

double h_04(double x) {
  // h_04(x) = ( x + sin ( x ) ) * exp ( - x * x ).
  return (x + sin(x)) * exp(-x * x);
}

double h_05(double x) {
  // h_05(x) = ( x - sin ( x ) ) * exp ( - x * x ).
  return (x - sin(x)) * exp(-x * x);
}
