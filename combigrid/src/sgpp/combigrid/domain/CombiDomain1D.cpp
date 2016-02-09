/* ****************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Janos Benk (benk@in.tum.de)
// @author Christoph Kowitz (kowitz@in.tum.de)
#include <sgpp/combigrid/domain/CombiDomain1D.hpp>
#include <vector>
#include <algorithm>

combigrid::Domain1D::Domain1D(int level, double min, double max,
                              combigrid::AbstractStretchingMaker& stretching) {
  if (level < 0) {
    std::cout << "Unsupported Operation Exception -> Cannot initialize a "
                 "Domain with less than 1 point!";
    throw UNSUPPORTED_OPERATION_EXCEPTION;
  }

  _isStretched = true;
  _level = level;
  stretching.get1DStretching(level, min, max, _stretching, _jacobian);
  _min = min;
  _max = max;
  _stretching_type = stretching.getStretchingType();
}

combigrid::Domain1D::Domain1D(const Domain1D& domain) {
  _min = domain.getMinDomain();
  _max = domain.getMaxDomain();
  _level = domain.getLevel();
  _isStretched = domain.isAxisScaled();
  _stretching_type = domain.getStretchingType();
  _stretching = domain.axisScaling();
  _jacobian = domain.axisJacobian();
}

void combigrid::Domain1D::transformRealToUnit(double coordReal,
                                              double& coordUnit, int level_in,
                                              bool noBoundary) const {
  // int verb = 6;
  if (_isStretched) {
    int startInd = 0, mid = 0;
    int endInd = static_cast<int>(_stretching.size()) - 1;
    double intersec = 0.0;
    int level_diff = (_level < level_in) ? 0 : _level - level_in;

    // stop when the difference is one, which means we found the cell
    while (endInd - startInd > combigrid::powerOfTwo[level_diff]) {
      mid = ((endInd + startInd) / 2);

      // make the bisection
      if (_stretching[mid] < coordReal) {
        startInd = mid;
      } else {
        endInd = mid;
      }
    }

    // startInd should be now at the beginning of the cell
    intersec = (_stretching[endInd] - coordReal) /
               (_stretching[endInd] - _stretching[startInd]);
    //   - make the transformation to the unit domain -> unitCoords
    coordUnit =
        (((double)endInd) * combigrid::oneOverPowOfTwo[level_diff] - intersec) /
        combigrid::powerOfTwo[_level - level_diff];

    if (noBoundary) {
      int offs = combigrid::powerOfTwo[level_diff];

      // for boundary cells make extrapolation
      if (startInd == 0) {
        double h = _stretching[endInd + offs] - _stretching[startInd + offs];
        intersec = (_stretching[endInd] - coordReal) / h;
        coordUnit = (((double)endInd) * combigrid::oneOverPowOfTwo[level_diff] -
                     intersec) /
                    combigrid::powerOfTwo[_level - level_diff];
        // COMBIGRID_OUT_LEVEL3( verb , " combigrid::Domain1D::tran 1");
      }

      if (endInd == (int)_stretching.size() - 1) {
        double h = _stretching[startInd] - _stretching[startInd - offs];
        intersec = (coordReal - _stretching[startInd]) / h;
        coordUnit = (((double)(endInd - offs)) *
                         combigrid::oneOverPowOfTwo[level_diff] +
                     intersec) /
                    combigrid::powerOfTwo[_level - level_diff];
        // COMBIGRID_OUT_LEVEL3( verb , " combigrid::Domain1D::tran 2 ,
        // intersec:" << intersec << " , endInd:"<<endInd);
        // COMBIGRID_OUT_LEVEL3( verb , " stretching_[startInd]:" <<
        // stretching_[startInd]<< " , stretching_[endInd]:" <<
        // stretching_[endInd]);
        // COMBIGRID_OUT_LEVEL3( verb , " h:" << h << " , h_old:" <<
        // stretching_[endInd] - stretching_[startInd] );
      }
    }

    // COMBIGRID_OUT_LEVEL3( verb , " combigrid::Domain1D::transformRealToUnit
    // coordReal:" << coordReal << " coordUnit:"
    //    << coordUnit << "  level_in:"<<level_in);
  } else {
    // no stretching , just simple scaling
    coordUnit = (coordReal - _min) / (_max - _min);
    // COMBIGRID_OUT_LEVEL3( verb , " combigrid::Domain1D::transformRealToUnit
    // NO STRETCH coordReal:" << coordReal << " coordUnit:"
    //    << coordUnit << "  level_in:"<<level_in);
  }
}

void combigrid::Domain1D::transformUnitToReal(int level, int index,
                                              double& realCoord) const {
  if (_isStretched) {
    // get the stretched index
    int level_diff = (_level < level) ? 0 : _level - level;
    realCoord = _stretching[index * combigrid::powerOfTwo[level_diff]];
  } else {
    // no stretching , just simple scaling
    realCoord = _min + (_max - _min) * ((double)index) * oneOverPowOfTwo[level];
  }
}

void combigrid::Domain1D::findEntry(double coordReal, int level_in,
                                    int& startIndex, double& intersect) const {
  if (_isStretched) {
    int startInd = 0, mid = 0;
    int endInd = static_cast<int>(_stretching.size()) - 1;
    int level_diff = (_level < level_in) ? 0 : _level - level_in;

    // stop when the difference is one, which means we found the cell
    while (endInd - startInd > combigrid::powerOfTwo[level_diff]) {
      mid = ((endInd + startInd) / 2);

      // make the bisection
      if (_stretching[mid] < coordReal) {
        startInd = mid;
      } else {
        endInd = mid;
      }
    }

    // startInd should be now at the beginning of the cell
    intersect = (coordReal - _stretching[startInd]) / _stretching[endInd] -
                (_stretching[endInd] - _stretching[startInd]);
    // this must be the start index of the grid, not from the stretching (the
    // levels could be different)
    startIndex = startInd / combigrid::powerOfTwo[level_diff];
  } else {
    // for the non-stretched case
    double unitC = (coordReal - _min) / (_max - _min);
    startIndex = static_cast<int>(
        ::floor((double)combigrid::powerOfTwo[level_in] * unitC));
    startIndex = (startIndex < 0) ? 0 : startIndex;
    startIndex = (startIndex >= combigrid::powerOfTwo[level_in] - 1)
                     ? combigrid::powerOfTwo[level_in] - 1
                     : startIndex;
    intersect = (unitC * (double)combigrid::powerOfTwo[level_in] -
                 (double)(startIndex));
  }
}

void combigrid::Domain1D::getMeshWidth(int index, int level_in, double& h0,
                                       double& h1) const {
  if (_isStretched) {
    int level_diff = (_level < level_in) ? 0 : _level - level_in;
    // index checking should be done in the debug mode
    h0 = _stretching[index * combigrid::powerOfTwo[level_diff]] -
         _stretching[(index - 1) * combigrid::powerOfTwo[level_diff]];
    h1 = _stretching[(index + 1) * combigrid::powerOfTwo[level_diff]] -
         _stretching[index * combigrid::powerOfTwo[level_diff]];
  } else {
    h1 = h0 = (_max - _min) * (1 / (double)combigrid::powerOfTwo[level_in]);
  }
}
