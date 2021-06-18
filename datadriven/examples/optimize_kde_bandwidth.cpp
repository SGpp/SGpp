// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>
#include <sgpp/datadriven/application/KernelDensityEstimator.hpp>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <random>
#include <vector>

using sgpp::datadriven::KernelType;
using sgpp::datadriven::BandwidthOptimizationType;
using sgpp::datadriven::KernelDensityEstimator;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;

/**
 * Prints a separator line.
 */
void printLine() {
  std::cout << "----------------------------------------"
               "----------------------------------------\n";
}

void randn(DataVector& rvar, std::mt19937& generator) {
  std::normal_distribution<double> distribution(0.5, 0.02);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
  }
}

void randn(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();

  std::mt19937 generator(static_cast<std::mt19937::result_type>(seedValue));
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randn(sample, generator);
    rvar.setRow(i, sample);
  }
}

/**
 * Main method.
 *
 * @param argc ignored
 * @param argv ignored
 */
int main(int argc, const char* argv[]) {
  (void)argc;
  (void)argv;

  std::cout << "sgpp::optimization kde bandwidth optimization.\n\n";

  // dimension of domain
  const size_t d = 1;

  DataVector bandwidths(d);

  // estimate a kernel density
  DataMatrix samples(1000, d);
  randn(samples);

  KernelDensityEstimator kdeRot(samples, KernelType::GAUSSIAN,
                                BandwidthOptimizationType::SILVERMANSRULE);
  kdeRot.getBandwidths(bandwidths);
  std::cout << "h_rot = " << bandwidths.toString() << " -> " << kdeRot.crossEntropy(samples)
            << std::endl;

  KernelDensityEstimator kdeOpt(samples, KernelType::GAUSSIAN,
                                BandwidthOptimizationType::MAXIMUMLIKELIHOOD);
  kdeOpt.getBandwidths(bandwidths);
  std::cout << "h_opt = " << bandwidths.toString() << " -> " << kdeOpt.crossEntropy(samples)
            << std::endl;

  return 0;
}
