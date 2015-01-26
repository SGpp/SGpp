// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef ADAPTIVESERIALCOMBIGRIDVARIABLECOEFFICIENTS_HPP_
#define ADAPTIVESERIALCOMBIGRIDVARIABLECOEFFICIENTS_HPP_

#include <sgpp/combigrid/combigrid/SerialCombiGrid.hpp>
#include <sgpp/combigrid/utils/combigrid_ultils.hpp>

namespace combigrid {

  class AdaptiveSerialCombiGridVariableCoefficients: public combigrid::SerialCombiGrid {
    public:
      AdaptiveSerialCombiGridVariableCoefficients(
        const CombiSchemeBasis* combischeme,
        const std::vector<bool>& hasBoundaryPts) :
        SerialCombiGrid(combischeme, hasBoundaryPts) {
        ;
      }
      //  virtual ~AdaptiveSerialCombiGridVariableCoefficients();

      void changeCoefficients(std::vector<double> newCoef);
      void changeCoefficients(int i, double newCoef);
  };

} /* namespace combigrid */
#endif /* ADAPTIVESERIALCOMBIGRIDVARIABLECOEFFICIENTS_HPP_ */