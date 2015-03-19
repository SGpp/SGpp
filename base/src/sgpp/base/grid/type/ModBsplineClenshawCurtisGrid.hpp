// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODBSPLINECLENSHAWCURTISGRID_HPP
#define MODBSPLINECLENSHAWCURTISGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Grid with modified Clenshaw-Curtis Bspline basis functions
     */
    class ModBsplineClenshawCurtisGrid : public Grid {
      protected:
        ModBsplineClenshawCurtisGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with modified Clenshaw-Curtis Bspline basis functions
         *
         * @param dim the dimension of the grid
           * @param degree the bspline's degree
         */
        ModBsplineClenshawCurtisGrid(size_t dim, size_t degree);

        /**
         * Destructor
         */
        virtual ~ModBsplineClenshawCurtisGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);

        virtual void serialize(std::ostream& ostr);
        virtual size_t getDegree();

      protected:
        // degree of Bspline
        size_t degree;
        const SBsplineModifiedClenshawCurtisBase* basis_;
    };

  }
}

#endif /* MODBSPLINECLENSHAWCURTISGRID_HPP */
