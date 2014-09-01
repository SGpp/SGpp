/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include <cstring>

#include "opt/operation/OpFactory.hpp"
#include "base/exception/factory_exception.hpp"

#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/LinearClenshawCurtisGrid.hpp"
#include "base/grid/type/ModLinearGrid.hpp"
#include "base/grid/type/BsplineGrid.hpp"
#include "base/grid/type/BsplineTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/BsplineClenshawCurtisGrid.hpp"
#include "base/grid/type/ModBsplineGrid.hpp"
#include "base/grid/type/WaveletGrid.hpp"
#include "base/grid/type/WaveletTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/ModWaveletGrid.hpp"

#include "opt/basis/bspline/boundary/operation/OperationMultipleHierarchisationBsplineBoundary.hpp"
#include "opt/basis/bspline/clenshawcurtis/operation/OperationMultipleHierarchisationBsplineClenshawCurtis.hpp"
#include "opt/basis/bspline/modified/operation/OperationMultipleHierarchisationModBspline.hpp"
#include "opt/basis/bspline/noboundary/operation/OperationMultipleHierarchisationBspline.hpp"
#include "opt/basis/linear/boundary/operation/OperationMultipleHierarchisationLinearBoundary.hpp"
#include "opt/basis/linear/clenshawcurtis/operation/OperationMultipleHierarchisationLinearClenshawCurtis.hpp"
#include "opt/basis/linear/modified/operation/OperationMultipleHierarchisationModLinear.hpp"
#include "opt/basis/linear/noboundary/operation/OperationMultipleHierarchisationLinear.hpp"
#include "opt/basis/wavelet/boundary/operation/OperationMultipleHierarchisationWaveletBoundary.hpp"
#include "opt/basis/wavelet/modified/operation/OperationMultipleHierarchisationModWavelet.hpp"
#include "opt/basis/wavelet/noboundary/operation/OperationMultipleHierarchisationWavelet.hpp"

namespace sg {
  namespace op_factory {

    opt::OperationMultipleHierarchisation* createOperationMultipleHierarchisation(base::Grid& grid) {
      if (std::strcmp(grid.getType(), "linear") == 0) {
        return new opt::OperationMultipleHierarchisationLinear(
                 dynamic_cast<base::LinearGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        return new opt::OperationMultipleHierarchisationLinearBoundary(
                 dynamic_cast<base::LinearTrapezoidBoundaryGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0) {
        return new opt::OperationMultipleHierarchisationLinearClenshawCurtis(
                 dynamic_cast<base::LinearClenshawCurtisGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "modlinear") == 0) {
        return new opt::OperationMultipleHierarchisationModLinear(
                 dynamic_cast<base::ModLinearGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "Bspline") == 0) {
        return new opt::OperationMultipleHierarchisationBspline(
                 dynamic_cast<base::BsplineGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0) {
        return new opt::OperationMultipleHierarchisationBsplineBoundary(
                 dynamic_cast<base::BsplineTrapezoidBoundaryGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0) {
        return new opt::OperationMultipleHierarchisationBsplineClenshawCurtis(
                 dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "modBspline") == 0) {
        return new opt::OperationMultipleHierarchisationModBspline(
                 dynamic_cast<base::ModBsplineGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "Wavelet") == 0) {
        return new opt::OperationMultipleHierarchisationWavelet(
                 dynamic_cast<base::WaveletGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0) {
        return new opt::OperationMultipleHierarchisationWaveletBoundary(
                 dynamic_cast<base::WaveletTrapezoidBoundaryGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "modWavelet") == 0) {
        return new opt::OperationMultipleHierarchisationModWavelet(
                 dynamic_cast<base::ModWaveletGrid&>(grid));
      } else {
        throw base::factory_exception(
          "OperationMultipleHierarchisation is not implemented for this grid type.");
      }
    }

  }
}
