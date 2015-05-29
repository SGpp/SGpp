// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FUNDAMENTALSPLINEGRID_HPP
#define FUNDAMENTALSPLINEGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    class FundamentalSplineGrid : public Grid {
      protected:
        FundamentalSplineGrid(std::istream& istr);

      public:
        FundamentalSplineGrid(size_t dim, size_t degree);

        virtual ~FundamentalSplineGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);

        virtual void serialize(std::ostream& ostr);
        virtual size_t getDegree();

      protected:
        size_t degree;

        const SFundamentalSplineBase* basis_;
    };

  }
}

#endif /* FUNDAMENTALSPLINEGRID_HPP */
