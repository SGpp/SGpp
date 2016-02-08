// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYCONDITIONALKDE_HPP_
#define OPERATIONDENSITYCONDITIONALKDE_HPP_

#include <sgpp/datadriven/application/GaussianKDE.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

class OperationDensityConditionalKDE {
 public:
  OperationDensityConditionalKDE(datadriven::GaussianKDE& kde);
  virtual ~OperationDensityConditionalKDE();

  /**
   * Conditional (Density) Functions
   *
   * @param mdim Marginalize in dimension mdim
   * @param xbar Point at which to conditionalize
   * @param conditionalizedKDE conditionalized kernel density
   */
  virtual void doConditional(size_t mdim, float_t xbar,
                             datadriven::GaussianKDE& conditionalizedKDE);

  /**
   * Conditional (Density) Functions
   *
   * @param mdims Conditionalize in dimensions mdims
   * @param xbar Point at which to conditionalize
   * @param conditionalizedKDE conditionalized kernel density
   */
  virtual void doConditional(std::vector<size_t>& mdims,
                             base::DataVector& xbar,
                             datadriven::GaussianKDE& conditionalizedKDE);

  /**
   * Conditional (Density) Functions to a 1d density where the remaining
   * dimension is mdim
   *
   * @param mdim conditionalize all dimensions but mdim
   * @param xbar point at which to conditionalize
   * @param conditionalizedKDE conditionalized kernel density
   */
  virtual void condToDimX(size_t mdim, base::DataVector& xbar,
                          datadriven::GaussianKDE& conditionalizedKDE);

  /**
   * Conditional (Density) Functions to a dd density where the remaining
   * dimensions are mdims
   *
   * @param mdims conditionalize all dimensions but the ones in mdims
   * @param xbar point at which to conditionalize
   * @param conditionalizedKDE conditionalized kernel density
   */
  virtual void condToDimXs(std::vector<size_t>& mdims, base::DataVector& xbar,
                           datadriven::GaussianKDE& conditionalizedKDE);

 private:
  datadriven::GaussianKDE* kde;
};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* OPERATIONDENSITYCONDITIONALKDE_HPP_ */
