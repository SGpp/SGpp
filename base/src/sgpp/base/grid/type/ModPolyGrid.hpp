// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef MODPOLYGRID_HPP
#define MODPOLYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/basis/modpoly/ModifiedPolyBasis.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with modified polynomial base functions
     */
    class ModPolyGrid : public Grid {
      protected:
        ModPolyGrid(std::istream& istr);

      public:
        /**
         * Constructor of grid with modified polynomial base functions
         *
         * @param dim the dimension of the grid
         * @param degree the max. polynom's degree
         */
        ModPolyGrid(size_t dim, size_t degree);

        /**
         * Destructor
         */
        virtual ~ModPolyGrid();

        virtual const char* getType();
        virtual void serialize(std::ostream& ostr);

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);
        virtual size_t getDegree() const;

      protected:
        /// max. polynom's degree
        size_t degree;
        const SModPolyBase* basis_;
    };

  }
}

#endif /* MODPOLYGRID_HPP */