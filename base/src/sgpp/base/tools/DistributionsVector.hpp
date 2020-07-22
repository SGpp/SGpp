// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <vector>

namespace sgpp {
namespace base {

class DistributionsVector {
 public:
  DistributionsVector();
  explicit DistributionsVector(size_t dim);
  DistributionsVector(size_t dim, std::shared_ptr<sgpp::base::Distribution> pdf);
  DistributionsVector(const DistributionsVector& other);

  virtual ~DistributionsVector();

  std::vector<std::shared_ptr<sgpp::base::Distribution>> getDistributions();

  sgpp::base::DataVector sample() const;

  sgpp::base::DataMatrix getBounds() const;

  void push_back(std::shared_ptr<sgpp::base::Distribution> pdf);
  std::shared_ptr<sgpp::base::Distribution> get(size_t i);
  size_t getSize() const;
  void clear();

 private:
  std::vector<std::shared_ptr<sgpp::base::Distribution>> distributions;
};

}  // namespace base
} /* namespace sgpp */
