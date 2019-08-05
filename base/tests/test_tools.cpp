// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>

#include <string>
#include <vector>

using sgpp::base::Printer;
using sgpp::base::RandomNumberGenerator;

double calculateMean(std::vector<double>& x) {
  double mean = 0.0;

  for (size_t i = 0; i < x.size(); i++) {
    mean += x[i];
  }

  mean /= static_cast<double>(x.size());
  return mean;
}

double calculateVariance(std::vector<double>& x) {
  double mean = calculateMean(x);
  double var = 0.0;

  for (size_t i = 0; i < x.size(); i++) {
    var += (x[i] - mean) * (x[i] - mean);
  }

  var /= static_cast<double>(x.size());
  return var;
}

BOOST_AUTO_TEST_CASE(TestPrinter) {
  // Test sgpp::base::Printer.

  // redirect std::cout to not confuse the user
  const std::string fileName = "test_tools_output.tmp";
  std::ofstream outStream(fileName);
  Printer::getInstance().setStream(&outStream);

  Printer::getInstance().setVerbosity(2);
  Printer::getInstance().printStatusBegin("Testing status printing...");
  Printer::getInstance().printStatusUpdate("Test status update 1");
  Printer::getInstance().printStatusUpdate("Test status update 2");
  Printer::getInstance().printStatusEnd("Testing status printing ended.");

  Printer::getInstance().getMutex();

  const double duration = Printer::getInstance().getLastDurationSecs();
  BOOST_CHECK_GE(duration, 0.0);
  BOOST_CHECK_LE(duration, 0.01);

  sgpp::base::DataMatrix A(3, 3, 0.0);
  A(0, 1) = 12.3;
  A(1, 2) = 42.1337;
  sgpp::base::FullSLE sle(A);
  Printer::getInstance().printSLE(sle);

  // undo redirection
  outStream.close();
  std::remove(fileName.c_str());
  Printer::getInstance().setStream(&std::cout);
}

BOOST_AUTO_TEST_CASE(TestRandomNumberGenerator) {
  // Test sgpp::base::RandomNumberGenerator.
  const size_t seed = 42;
  const size_t N = 20000;
  std::vector<double> numbers(N);

  // set and test seed getting/setting
  RandomNumberGenerator::getInstance().setSeed();
  RandomNumberGenerator::getInstance().setSeed(seed);
  BOOST_CHECK_EQUAL(RandomNumberGenerator::getInstance().getSeed(), seed);

  // test continuous uniform random numbers
  {
    for (size_t i = 0; i < N; i++) {
      numbers[i] = RandomNumberGenerator::getInstance().getUniformRN();
      BOOST_CHECK_GE(numbers[i], 0.0);
      BOOST_CHECK_LE(numbers[i], 1.0);
    }

    BOOST_CHECK_SMALL(calculateMean(numbers) - 0.5, 1e-3);
    BOOST_CHECK_SMALL(calculateVariance(numbers) - 1.0 / 12.0, 1e-3);
  }

  // test Gaussian random numbers
  {
    std::vector<double> mus = {0.0, 12.3, -42.0, 13.37};
    std::vector<double> sigmas = {1.0, 2.6, 8.1, 0.3};

    for (size_t k = 0; k < mus.size(); k++) {
      for (size_t i = 0; i < N; i++) {
        numbers[i] = RandomNumberGenerator::getInstance().getGaussianRN(mus[k], sigmas[k]);
      }

      BOOST_CHECK_SMALL(calculateMean(numbers) - mus[k], 0.1 * sigmas[k]);
      BOOST_CHECK_SMALL(calculateVariance(numbers) - sigmas[k] * sigmas[k],
                        0.1 * sigmas[k] * sigmas[k]);
    }
  }

  // test discrete uniform random numbers
  for (size_t k = 1; k < 11; k++) {
    for (size_t i = 0; i < N; i++) {
      numbers[i] = static_cast<double>(RandomNumberGenerator::getInstance().getUniformIndexRN(k));
      BOOST_CHECK_EQUAL(numbers[i], static_cast<int>(numbers[i]));
      BOOST_CHECK_GE(numbers[i], 0);
      BOOST_CHECK_LE(numbers[i], static_cast<double>(k - 1));
    }

    double kDbl = static_cast<double>(k);
    BOOST_CHECK_SMALL(calculateMean(numbers) - (kDbl - 1.0) / 2.0, 0.01 * kDbl);
    BOOST_CHECK_SMALL(calculateVariance(numbers) - (kDbl * kDbl - 1.0) / 12.0, 0.01 * kDbl * kDbl);
  }
}
