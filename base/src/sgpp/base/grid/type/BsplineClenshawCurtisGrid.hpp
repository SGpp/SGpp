// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BSPLINECLENSHAWCURTISGRID_HPP
#define BSPLINECLENSHAWCURTISGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Grid with Clenshaw-Curtis Bspline basis functions with boundaries, pentagon cut
     */
    class BsplineClenshawCurtisGrid : public Grid {
      protected:
        /**
         * This constructor creates a new GridStorage out of the stream.
         *
         * @param istr inputstream that contains the grid information
         */
        BsplineClenshawCurtisGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with Clenshaw-Curtis Bspline basis functions with boundaries, pentagon cut
         *
         * @param dim the dimension of the grid
         * @param degree the bspline's degree
         * @param boundaryLevel 1 + how much levels the boundary is coarser than
         *                      the main axes, 0 means one level finer,
         *                      1 means same level,
         *                      2 means one level coarser, etc.
         */
        BsplineClenshawCurtisGrid(size_t dim,
                                  size_t degree,
                                  level_t boundaryLevel = 1);

        /**
         * Destructor.
         */
        virtual ~BsplineClenshawCurtisGrid() override;

        /**
         * @return string that identifies the grid type uniquely
         */
        virtual SGPP::base::GridType getType() override;

        /**
         * @return B-spline basis
         */
        virtual const SBasis& getBasis() override;

        /**
         * @return pointer to a GridGenerator object
         */
        virtual GridGenerator* createGridGenerator() override;

        /**
         * reads a grid out of a string
         *
         * @param istr string that contains the grid information
         * @return grid
         */
        static Grid* unserialize(std::istream& istr);

        /**
         * Serializes the grid.
         *
         * @param ostr stream to which the grid is written
         */
        virtual void serialize(std::ostream& ostr) override;

        /**
         * @return B-spline degree
         */
        virtual size_t getDegree();

      protected:
        /// B-spline degree
        size_t degree;
        /// B-spline basis
        const SBsplineClenshawCurtisBase* basis_;
        /// 1 + how much levels the boundary is coarser than the main axes
        level_t boundaryLevel;
    };

  }
}

#endif /* BSPLINECLENSHAWCURTISGRID_HPP */
