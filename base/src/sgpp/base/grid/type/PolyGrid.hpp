// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POLYGRID_HPP
#define POLYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with polynomial base functions
     */
    class PolyGrid : public Grid {
      protected:
        PolyGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with polynomial base functions
         *
         * @param dim the dimension of the grid
         * @param degree the max. polynom's degree
         */
        PolyGrid(size_t dim, size_t degree);

        /**
         * Destructor
         */
        virtual ~PolyGrid();

        virtual const char* getType();
        virtual const SBasis& getBasis();
        virtual void serialize(std::ostream& ostr);

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
        size_t getDegree() const;

      protected:
        /// max. polynom's degree
        size_t degree;
        const SPolyBase* basis_;
    };

  }
}

#endif /* POLYGRID_HPP */
