// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POLYTRUNCATEDBOUNDARYGRID_HPP
#define POLYTRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * trapezoid boundary grid with polynomial base functions
     */
    class PolyBoundaryGrid : public Grid {
      protected:
        PolyBoundaryGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with polynomial base functions
         *
         * @param dim the dimension of the grid
         * @param degree the max. polynom's degree
         * @param boundaryLevel 1 + how much levels the boundary is coarser than
         *                      the main axes, 0 means one level finer,
         *                      1 means same level,
         *                      2 means one level coarser, etc.
         */
        PolyBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel = 1);

        /**
         * Destructor
         */
        virtual ~PolyBoundaryGrid() override;

        virtual const SBasis& getBasis() override;
        virtual SGPP::base::GridType getType() override;
        virtual void serialize(std::ostream& ostr) override;

        virtual GridGenerator* createGridGenerator() override;

        static Grid* unserialize(std::istream& istr);
        size_t getDegree() const;

      protected:
        /// max. polynom's degree
        size_t degree;
        /// polynomial basis
        const SPolyBoundaryBase* basis_;
        /// 1 + how much levels the boundary is coarser than the main axes
        level_t boundaryLevel;
    };

  }
}

#endif /* POLYTRUNCATEDBOUNDARYGRID_HPP */
