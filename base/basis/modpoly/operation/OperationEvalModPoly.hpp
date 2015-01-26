/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONEVALMODPOLY_HPP
#define OPERATIONEVALMODPOLY_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/modpoly/ModifiedPolyBasis.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace base {

    /**
     * This class implements OperationEval for a grids with mod poly basis ansatzfunctions with
     *
     * @version $HEAD$
     */
    class OperationEvalModPoly : public OperationEval {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         * @param degree the polynom's max. degree
         */
        OperationEvalModPoly(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

        /**
         * Destructor
         */
        virtual ~OperationEvalModPoly() {}

        virtual double eval(DataVector& alpha, std::vector<double>& point);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
        /// Mod Poly Basis object
        sg::base::SModPolyBase base;
    };

  }
}

#endif /* OPERATIONEVALMODPOLY_HPP */
