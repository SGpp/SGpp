// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/L2LejaPointDistribution.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <cmath>
#include <vector>
#include <iostream>

/**
 * Warning: This tests are made for a epsilon of 0.00001 in the leja class (for the optimizer)!
 */

// in %
const double tolerance = 0.01;

// BOOST_AUTO_TEST_CASE(testLejaPoints) {
//  std::function<double(double)> weight_func = [](double d) { return 1; };
//  double start = 0.5;
//  double lower_bound = 0;
//  double upper_bound = 1;
//  int number = 8;
//
//  std::vector<double> leja_points;
//
//  // start point
//  leja_points.push_back(start);
//  auto points = leja_points;
//
//  // calc leja points
//  sgpp::combigrid::LejaPointDistribution::calc_leja_points(leja_points, points, number,
//  lower_bound,
//                                                           upper_bound, weight_func);
//
//  // correct solution:
//  std::vector<double> correct_leja_points;
//  correct_leja_points.push_back(1.0428174199820221e-05);
//  correct_leja_points.push_back(0.065002567709224635);
//  correct_leja_points.push_back(0.17065421163314368);
//  correct_leja_points.push_back(0.34719663104769344);
//  correct_leja_points.push_back(0.5);
//  correct_leja_points.push_back(0.66085325446028498);
//  correct_leja_points.push_back(0.78867243975866386);
//  correct_leja_points.push_back(0.91962400339636541);
//  correct_leja_points.push_back(0.99998958672635252);
//
//  // test number of points
//  BOOST_CHECK_EQUAL(leja_points.size(), correct_leja_points.size());
//
//  for (size_t i = 0; i < correct_leja_points.size(); ++i) {
//    BOOST_CHECK_CLOSE(leja_points.at(i), correct_leja_points.at(i), tolerance);
//  }
//}

BOOST_AUTO_TEST_CASE(testLejaDistributionStartingPoint) {
  sgpp::combigrid::LejaPointDistribution leja;

  double first_point = leja.compute(0, 0);

  BOOST_CHECK_CLOSE(first_point, 0.5, tolerance);
}

/*
 * the starting point for this function is 1.0
 */
double quadratic_func(double x) { return 42 * x * x; }

BOOST_AUTO_TEST_CASE(testLejaOtherStartingPoint) {
  auto w = sgpp::combigrid::SingleFunction(quadratic_func);
  sgpp::combigrid::LejaPointDistribution leja(w);

  // check if first point is correct, it should be 1.0
  BOOST_CHECK_CLOSE(leja.compute(0, 0), 1.0, tolerance);
}

double sinusweight(double x) { return std::sin(x * 3.14159265358979323846); }

BOOST_AUTO_TEST_CASE(testLejaSinusWithNormalDistribution) {
  auto w = sgpp::combigrid::SingleFunction(sinusweight);
  sgpp::combigrid::LejaPointDistribution leja(w);

  // check if first point is correct, it should be 0.5
  BOOST_CHECK_CLOSE(leja.compute(0, 0), 0.5, tolerance);

  // check the next few points
  BOOST_CHECK_CLOSE(leja.compute(1, 1), 1 - 0.22614641471525609, tolerance);
  BOOST_CHECK_CLOSE(leja.compute(1, 2), 1 - 0.81695850918945612, tolerance);
  BOOST_CHECK_CLOSE(leja.compute(1, 3), 1 - 0.088181946121442964, tolerance);
  BOOST_CHECK_CLOSE(leja.compute(1, 4), 1 - 0.92906611726964328, tolerance);
}

BOOST_AUTO_TEST_CASE(testL2LejaPoints) {
  sgpp::combigrid::L2LejaPointDistribution l2leja;
  sgpp::combigrid::PolynomialQuadratureEvaluator quad;

  std::vector<double> points_all{
      0.5,       0.875,       0.101776,    0.281491,   0.973826,  0.0198389, 0.705177,   0.388289,
      0.994366,  0.174955,    0.793351,    0.00414731, 0.605696,  0.930801,  0.0569999,  0.229032,
      0.835673,  0.445066,    0.998797,    0.0365097,  0.655543,  0.136822,  0.953319,   0.334706,
      0.750183,  0.000857035, 0.552322,    0.904909,   0.0788431, 0.985365,  0.202712,   0.41685,
      0.814875,  0.0111596,   0.579507,    0.963983,   0.306817,  0.118908,  0.681575,   0.890453,
      0.0279353, 0.999744,    0.254479,    0.474651,   0.770641,  0.0471021, 0.941875,   0.156159,
      0.630539,  0.361311,    0.000183025, 0.990173,   0.855379,  0.525842,  0.0684998,  0.727994,
      0.189255,  0.9185,      0.00742122,  0.979864,   0.320771,  0.0907109, 0.459864,   0.996876,
      0.668604,  0.0154277,   0.82542,     0.241806,   0.565906,  0.958724,  0.127814,   0.402164,
      0.865481,  0.00230663,  0.7816,      0.268797,   0.969429,  0.0416694, 0.618096,   0.214889,
      0.91183,   0.512777,    0.0238889,   0.992424,   0.716941,  0.109749,  0.373747,   0.845787,
      0.0627104, 0.999946,    0.431229,    0.165108,   0.740011,  0.936467,  0.00568931, 0.592392,
      0.294422,  0.883438,    0.0846337,   0.982736};

  std::vector<double> points_seq;
  for (size_t i = 0; i < points_all.size(); i++) {
    double x = l2leja.compute(i, i);
    points_seq.push_back(x);
    quad.setGridPoints(points_seq);

    BOOST_CHECK_CLOSE(x, points_all[i], 1e-3);
    BOOST_CHECK(quad.getAbsoluteWeightSum() < 1.1);
  }
}

#ifdef USE_DAKOTA

BOOST_AUTO_TEST_CASE(testWeightedL2LejaPoints) {
  sgpp::combigrid::CO2 co2Model;
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::BOUNDED_LOGNORMAL;
  config.polyParameters.logmean_ = co2Model.logmean;
  config.polyParameters.stddev_ = co2Model.stddev;
  config.polyParameters.lowerBound_ = co2Model.bounds[0];
  config.polyParameters.upperBound_ = co2Model.bounds[1];
  auto basisFunction = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  auto func = [basisFunction, &co2Model](double x_unit) {
    double x_prob = (co2Model.bounds[1] - co2Model.bounds[0]) * x_unit + co2Model.bounds[0];
    return basisFunction->pdf(x_prob);
  };

  sgpp::combigrid::SingleFunction weightFunction(func);
  sgpp::combigrid::L2LejaPointDistribution l2leja(weightFunction);
  sgpp::combigrid::PolynomialQuadratureEvaluator quad;

  std::vector<double> points_all{
      0.308723,    0.677028,  0.101446,   0.916106,  0.470405,    0.0274901,  0.981783,  0.199154,
      0.798541,    0.566859,  0.00617981, 0.99608,   0.381582,    0.146763,   0.860146,  0.0600752,
      0.736673,    0.953056,  0.253828,   0.617714,  0.00131865,  0.426464,   0.889892,  0.0800088,
      0.999166,    0.51922,   0.172939,   0.827948,  0.0162007,   0.968396,   0.34278,   0.706966,
      0.0430044,   0.935459,  0.22816,    0.64524,   0.123125,    0.989896,   0.494621,  0.769765,
      0.0105465,   0.282442,  0.875482,   0.403843,  0.000271479, 0.975404,   0.59115,   0.0699475,
      0.84305,     0.159872,  0.999822,   0.542737,  0.0349612,   0.753389,   0.325882,  0.94432,
      0.213659,    0.448352,  0.904206,   0.0912584, 0.661743,    0.00347681, 0.993176,  0.362353,
      0.812321,    0.0510128, 0.9611,     0.267645,  0.720919,    0.133934,   0.0213068, 0.925701,
      0.604162,    0.986099,  0.1862,     0.506768,  0.78435,     0.111688,   0.415073,  0.997838,
      0.00825609,  0.691735,  0.29556,    0.867974,  0.0311674,   0.554938,   0.897395,  0.240426,
      5.80238e-05, 0.972005,  0.631906,   0.0854532, 0.459364,    0.835485,   0.352627,  0.994696,
      0.0554528,   0.94869,   0.153321,   0.745215};

  std::vector<double> points_seq;
  for (size_t i = 0; i < points_all.size(); i++) {
    double x = l2leja.compute(i, i);
    points_seq.push_back(x);
    quad.setGridPoints(points_seq);
    BOOST_CHECK_CLOSE(x, points_all[i], 1e-3);
    BOOST_CHECK(quad.getAbsoluteWeightSum() < 1.4);
  }
}

#endif
