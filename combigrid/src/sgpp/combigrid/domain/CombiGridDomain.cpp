// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/domain/CombiGridDomain.hpp>

#include <vector>

combigrid::GridDomain::GridDomain(int dim, const std::vector<int>& levels,
                                  const std::vector<double>& min, const std::vector<double>& max,
                                  combigrid::AbstractStretchingMaker& stretchingMaker) {
  dim_ = dim;
  _stretching_type = UNKNOWN;
  _min = min;
  _max = max;

  // add for each dimension
  for (int d = 0; d < dim_; d++) {
    _axisDomains.push_back(combigrid::Domain1D(levels[d], _min[d], _max[d], stretchingMaker));
  }

  _stretching_type = stretchingMaker.getStretchingType();
}

combigrid::GridDomain::GridDomain(int dim, const std::vector<int>& levels,
                                  const std::vector<double>& min, const std::vector<double>& max,
                                  std::vector<AbstractStretchingMaker*> stretchingMaker) {
  dim_ = dim;
  _stretching_type = UNKNOWN;
  _min = min;
  _max = max;

  // add for each dimension
  for (int d = 0; d < dim_; d++) {
    _axisDomains.push_back(combigrid::Domain1D(levels[d], _min[d], _max[d], *stretchingMaker[d]));
  }

  _stretching_type = UNKNOWN;
}

/* copy constructor*/
combigrid::GridDomain::GridDomain(const GridDomain& domain) {
  dim_ = domain.getDim();
  _stretching_type = domain.getStretchingType();
  _max = domain.getMax();
  _min = domain.getMin();

  for (int d = 0; d < dim_; d++) {
    Domain1D newDomain(domain.get1DDomain(d));  // invoke the copy constructor of
    _axisDomains.push_back(newDomain);
  }
}

void combigrid::GridDomain::transformRealToUnit(std::vector<double>& coords,
                                                const std::vector<int>& levels_in,
                                                const std::vector<bool>& boundaryFlag) const {
  // for each dimension make the transformation
  // int verb = 6;
  double tmp = 0.0;

  // COMBIGRID_OUT_LEVEL3( verb , " combigrid::GridDomain::transformRealToUnit()
  // ");
  for (int d = 0; d < dim_; d++) {
    _axisDomains[d].transformRealToUnit(coords[d], tmp, levels_in[d], boundaryFlag[d]);
    coords[d] = tmp;
  }
}

void combigrid::GridDomain::printDomain() {
  std::cout << "------------------" << std::endl;

  for (int d = 0; d < dim_; d++) {
    std::cout << _axisDomains[d].getMinDomain() << "\t" << _axisDomains[d].getMaxDomain() << "\t"
              << _axisDomains[d].isAxisScaled() << std::endl;
  }

  std::cout << "------------------" << std::endl;
}
