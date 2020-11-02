// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONINVERSEROSENBLATTTRANSFORMATIONPOLYCLENSHAWCURTIS1D_HPP
#define OPERATIONINVERSEROSENBLATTTRANSFORMATIONPOLYCLENSHAWCURTIS1D_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTransformation1D.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/tools/Printer.hpp>

#include <sgpp/globaldef.hpp>
#include <map>
#include <vector>
#include <functional>

namespace sgpp {
namespace datadriven {

class OperationInverseRosenblattTransformation1DPolyClenshawCurtis
    : public OperationTransformation1D {
 private:
  base::GaussLegendreQuadRule1D gauss;
  base::DataVector weights;
  base::DataVector gauss_coordinates;
  std::unique_ptr<base::OperationEval> opEval;
  std::vector<double> patch_areas;
  std::vector<bool> is_negative_patch;
  std::vector<double> ordered_grid_points;
  std::vector<std::function<double(double)>> patch_functions;
  std::multimap<double, double> coord_cdf;
  double sum;
  size_t quadOrder;

  /**
   * this function computes the CDF (i.e. the patch areas)
   * and saves it into the vectors of this object
   */
  void init(base::DataVector* alpha1d);

  /**
   * this performs the actual sampling after the CDF has been computed
   * @param alpha1d coefficient vector in the current direction
   * @param coord1d point where to evaluate the CDF
   * @return the value of the CDF
   */
  double sample(base::DataVector* alpha1d, double coord1d);

 protected:
  base::Grid* grid;

 public:
  explicit OperationInverseRosenblattTransformation1DPolyClenshawCurtis(base::Grid* grid);
  virtual ~OperationInverseRosenblattTransformation1DPolyClenshawCurtis();

  /**
   * Inverse Rosenblatt Transformation 1D
   * @param alpha1d
   * @param coord1d
   * @return
   */
  double doTransformation1D(base::DataVector* alpha1d, double coord1d);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* OPERATIONINVERSEROSENBLATTTRANSFORMATIONPOLYCLENSHAWCURTIS1D_HPP */
