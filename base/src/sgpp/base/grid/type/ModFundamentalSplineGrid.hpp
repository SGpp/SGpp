// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODFUNDAMENTALSPLINEGRID_HPP
#define MODFUNDAMENTALSPLINEGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    class ModFundamentalSplineGrid : public Grid {
      protected:
        ModFundamentalSplineGrid(std::istream& istr);

      public:
        ModFundamentalSplineGrid(size_t dim, size_t degree);

        virtual ~ModFundamentalSplineGrid();

        virtual const char* getType();

        virtual const SBasis& getBasis();

        virtual GridGenerator* createGridGenerator();

        static Grid* unserialize(std::istream& istr);

        virtual void serialize(std::ostream& ostr);
        virtual size_t getDegree();

      protected:
        size_t degree;

        const SFundamentalSplineModifiedBase* basis_;
    };

  }
}

#endif /* MODFUNDAMENTALSPLINEGRID_HPP */
