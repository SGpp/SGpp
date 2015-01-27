// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALMODBSPLINE_HPP
#define OPERATIONEVALMODBSPLINE_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/basis/bsplineModified/BsplineModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements OperationEval for a grids with modified Bspline basis functions with a certain degree
     *
     * @version $HEAD$
     */
    class OperationEvalModBspline : public OperationEval {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         * @param degree the polynom's max. degree
         */
        OperationEvalModBspline(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

        /**
         * Destructor
         */
        virtual ~OperationEvalModBspline() {}

        virtual double eval(DataVector& alpha, std::vector<double>& point);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
        /// Mod Bspline Basis object
        SBsplineModifiedBase base;
    };

  }
}

#endif /* OPERATIONEVALMODBSPLINE_HPP */