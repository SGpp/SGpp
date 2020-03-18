// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REGULARIZATIONCONFIGURATION_HPP_
#define REGULARIZATIONCONFIGURATION_HPP_

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

enum class RegularizationMetricType { mse, nll, accuracy, residual };

enum class RegularizationType {
  Identity,
  Laplace,
  Diagonal,
  Lasso,
  ElasticNet,
  GroupLasso
};

struct RegularizationConfiguration {
  RegularizationType type_ = RegularizationType::Identity;
  double lambda_ = 0.01;
  double l1Ratio_ = 0.0;
  double exponentBase_ = 1.0;
  double lamda_start_ = 0.01;
  double lambda_end_ = 0.01;
  double lambda_steps_ = 0;
  bool lambda_log_scale_ = false;
  bool optimizeLambda_ = false;
  double optimizerTolerance_ = 1e-15;
  double convergenceThreshold_ = 1e-5;
  double intervalA_ = 1e-15;
  double intervalB_ = 1.0;

  RegularizationMetricType regularizationMetric_ =
      RegularizationMetricType::residual;

  /*
  // Debug method to neatly print internal data
  void dumpToStream(std::ostream& stream_out = std::cout) const {
    stream_out << "type: \t\t\t" << RegularizationTypeParser::toString(type_)
               << std::endl;
    stream_out << "lambda: \t\t\t" << lambda_ << std::endl;
    stream_out << "l1Ratio: \t\t\t" << l1Ratio_ << std::endl;
    stream_out << "exponentBase: \t\t" << exponentBase_ << std::endl;
    stream_out << "lambda_start: \t\t" << lamda_start_ << std::endl;
    stream_out << "lambda_end: \t\t" << lambda_end_ << std::endl;
    stream_out << "lambda_steps: \t\t" << lambda_steps_ << std::endl;
    stream_out << "lambda_log_scale: \t" << std::boolalpha << lambda_log_scale_
               << std::endl;
    stream_out << "optimizeLambda: \t" << std::boolalpha << optimizeLambda_
               << std::endl;
    stream_out << "optimizerTolerance: \t" << optimizerTolerance_ << std::endl;
    stream_out << "convergenceThreshold: \t" << convergenceThreshold_
               << std::endl;
    stream_out << "intervalA: \t\t" << intervalA_ << std::endl;
    stream_out << "intervalB: \t\t" << intervalB_ << std::endl;
  }
  */
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* REGULARIZATIONCONFIGURATION_HPP_ */
