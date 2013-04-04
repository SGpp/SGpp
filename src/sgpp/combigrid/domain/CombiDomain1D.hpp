/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef COMBIDOMAIN1D_HPP_
#define COMBIDOMAIN1D_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/domain/AbstractStretchingMaker.hpp"

namespace combigrid {

  /** */
  class Domain1D {
    public:

      /** */
      Domain1D(double min, double max);

      /** */
      Domain1D(const std::vector<double>& inputStretching);

      /** */
      Domain1D(int level, double min, double max, const AbstractStretchingMaker& stretching);

      virtual ~Domain1D() {
        ;
      }

      /** return if the axis is scaled */
      inline bool isAxisScaled() const {
        return isStretched_;
      }

      /** return the minimum of the domain */
      inline double getMinDomain() const {
        return min_;
      }

      /** return the maximum of the domain */
      inline double getMaxDomain() const {
        return max_;
      }

      /** return if the axis is scaled */
      const std::vector<double>& axisScaling() const {
        return stretching_;
      }

      /** return the level of the domain , then the number of points are (2^L) + 1*/
      inline int getLevel() const {
        return level_;
      }

      /** transform the real coordinates to unit coordinates
       * @param coordReal real non-unit coordinates
       * @param coordUnit coordinates in the unit cube
       * @param level_in level of the resolution required
       * @param noBoundary make extrapolation for the boundary cells*/
      void transformRealToUnit( double coordReal, double& coordUnit ,
                                int level_in = 0 , bool noBoundary = false ) const;

      /** transform from unit index to real coordinates
       * @param level input level
       * @param index the index 0..2^level
       * @param realCoord the real coordinate */
      void transformUnitToReal( int level , int index , double& realCoord) const;

      /** flocated the point on one axis
       * @param coordReal the coord on real domain
       * @param level_in input level
       * @param startIndex the left index of the cell
       * @param intersect the intersection of the cell
       */
      void findEntry(double coordReal, int level_in ,
                     int& startIndex , double& intersect) const;

      /** returns the mesh width /scaling
       * @param index point index
       * @param level_in the level resolution
       * @param h0 the first mesh width
       * @param h1 the second mesh width */
      void getMeshWidth(int index , int level_in , double& h0 , double& h1) const;

    private:

      /** the level in the case of stretched */
      int level_;

      /** if stretching is needed*/
      bool isStretched_;

      /** minimum value of the axis */
      double min_;

      /** minimum value of the axis */
      double max_;

      /** if the axis is scaled then here we store the scaling */
      std::vector<double> stretching_;

  };

}

#endif /* COMBIDOMAIN1D_HPP_ */
