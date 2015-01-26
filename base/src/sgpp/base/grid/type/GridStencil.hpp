/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Gerrit Buse (buse@in.tum.de)

#ifndef GRIDSTENCIL_HPP
#define GRIDSTENCIL_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * grid with linear base functions
     */
    class GridStencil : public Grid {
      protected:
        GridStencil(std::istream& istr);

      public:

        typedef std::vector<unsigned int> IndexStencil;

        typedef std::vector<float>      WeightStencil;


        /**
         * Constructor Linear Grid without boundaries
         *
         * @param dim the dimension of the grid
         */
        GridStencil(size_t dim);

        /**
         * Constructor Linear Grid
         *
         * @param BB the BoundingBox of the grid
         */
        GridStencil(BoundingBox& BB);

        /**
         * Destructor
         */
        virtual ~GridStencil();

        /**
         * Get the surplus stencil, in fact an array of unsigned ints.
         */
        virtual const IndexStencil&
        getSurplusStencil() const;

        /**
         * Get the neighbor stencil, in fact an array of unsigned ints.
         */
        virtual const IndexStencil&
        getNeighborStencil() const;

        /**
         * Get the weight stencil, in fact an array of floats.
         */
        virtual const WeightStencil&
        getWeightStencil() const;


      protected:

        IndexStencil surplusStencil;
        IndexStencil neighborStencil;
        WeightStencil weightStencil;
    };

  }
}

#endif /* GRIDSTENCIL_HPP */
