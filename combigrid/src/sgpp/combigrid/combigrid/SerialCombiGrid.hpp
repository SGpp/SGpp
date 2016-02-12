// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SERIALCOMBIGRID_HPP_
#define SERIALCOMBIGRID_HPP_

#include <sgpp/combigrid/combigrid/CombiGrid.hpp>
#include <sgpp/combigrid/domain/CombiGridDomain.hpp>
#include <sgpp/combigrid/utils/combigrid_utils.hpp>

#include <vector>

namespace combigrid {
template <typename _Tp>
class SerialCombiGrid : public CombiGrid<_Tp> {
 public:
  SerialCombiGrid(int in_dim, const std::vector<bool>& in_hasBoundaryPts)
      : CombiGrid<_Tp>(in_dim, in_hasBoundaryPts) {
    /**QUE NADA*/
  }

  /** sets the domain of all the full grids, this is the correct way for
   * extrapolation */
  void setDomainAllFG(GridDomain* gridDomain) const;

  /** create the actual vector for the full grids. <br>
   * This is be different for serial and parallel implementations */
  void createFullGrids();
};
}  // namespace combigrid

#endif /* SERIALCOMBIGRID_HPP_ */
