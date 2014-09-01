/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_BASE_GRID_TYPE_LINEARCLENSHAWCURTISGRID_HPP
#define SGPP_BASE_GRID_TYPE_LINEARCLENSHAWCURTISGRID_HPP

#include <iostream>

#include "base/grid/Grid.hpp"
#include "base/tools/CosineTable.hpp"

namespace sg {
  namespace base {

    /**
     * Clenshaw-Curtis grid with linear basis functions.
     */
    class LinearClenshawCurtisGrid : public Grid {
      public:
        /**
         * Constructor.
         *
         * @param dim       number of dimensions
         * @param cosine_table  lookup cosine table (optional)
         */
        LinearClenshawCurtisGrid(size_t dim, const CosineTable* cosine_table = NULL);

        /**
         * Destructor.
         */
        virtual ~LinearClenshawCurtisGrid();

        /**
         * @return  identifying grid type string
         */
        virtual const char* getType();

        /**
         * @return grid generator for this grid type
         */
        virtual GridGenerator* createGridGenerator();

        /**
         * @param istr  input stream containing the serialization
         * @return      pointer to newly generated deserialized grid
         */
        static Grid* unserialize(std::istream& istr);

        /**
         * @return  cosine table
         */
        virtual const CosineTable* getCosineTable() const;

        /**
         * @param cosine_table  new cosine table
         */
        virtual void setCosineTable(const CosineTable* cosine_table);

      protected:
        /// cosine table
        const CosineTable* cosine_table;

        /**
         * Deserialization constructor.
         *
         * @param istr  serialized grid
         */
        LinearClenshawCurtisGrid(std::istream& istr);
    };

  }
}

#endif
